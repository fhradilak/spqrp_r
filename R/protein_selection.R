# Pairwise random-forest classifier + 3 selectable backends.
# Port of spqrp/spqrp/protein_selection.py.

# Split a vector of unique patients into train/test.
# Uses withr::with_seed so the user's global RNG is not mutated
# (required by CRAN policy).
split_by_patient <- function(unique_individuals, test_size, seed = 42L) {
  n <- length(unique_individuals)
  n_test <- max(1L, round(n * test_size))
  test_idx <- withr::with_seed(seed, sample.int(n, size = n_test))
  list(
    train = unique_individuals[-test_idx],
    test  = unique_individuals[test_idx]
  )
}

# Pivot a long-format intensity frame to wide and prefix columns.
pivot_df <- function(df,
                      index = c(PATIENT, SAMPLE),
                      columns = "Protein",
                      values = "Intensity") {
  out <- tidyr::pivot_wider(
    df,
    id_cols = dplyr::all_of(index),
    names_from = dplyr::all_of(columns),
    values_from = dplyr::all_of(values)
  )
  tibble::as_tibble(out)
}

# Element-wise differences between two numeric vectors.
compute_feature_differences <- function(s1, s2) s1 - s2

# NaN-aware Euclidean distance (sklearn nan_euclidean_distances formula).
compute_nan_euclidean_distance <- function(s1, s2) {
  s1 <- as.numeric(s1); s2 <- as.numeric(s2)
  n <- length(s1)
  mask <- !is.na(s1) & !is.na(s2)
  m <- sum(mask)
  if (m == 0L) return(NA_real_)
  sq <- (s1[mask] - s2[mask]) ^ 2
  sqrt(n / m * sum(sq))
}

# Generate all i<j pairs of rows of X, with labels y[i]==y[j].
create_pairs <- function(X_data, y_data) {
  n <- nrow(X_data)
  pair_indices <- if (n >= 2L) {
    t(utils::combn(seq_len(n), 2L))
  } else {
    matrix(integer(0L), ncol = 2L)
  }
  if (nrow(pair_indices) == 0L) {
    return(list(pairwise = list(), labels = integer(),
                 pair_indices = pair_indices))
  }
  X_mat <- as.matrix(X_data)
  pairwise <- lapply(seq_len(nrow(pair_indices)), function(k) {
    list(X_mat[pair_indices[k, 1L], ], X_mat[pair_indices[k, 2L], ])
  })
  labels <- as.integer(
    y_data[pair_indices[, 1L]] == y_data[pair_indices[, 2L]]
  )
  list(pairwise = pairwise, labels = labels, pair_indices = pair_indices)
}

# Threshold selection on predicted probabilities; mirrors Python get_threshold.
#' Pick a binary-classifier threshold from probabilities
#'
#' Supports `"ROC"` (Youden's J via `pROC`), `"F1"` (best F1 via grid),
#' `"J"` (alias for ROC), and `"MinFP"` (largest threshold with FPR <= max_fpr).
#'
#' @param y_test Binary 0/1 vector.
#' @param y_prob Predicted probability of class 1.
#' @param method One of `"ROC"`, `"F1"`, `"J"`, `"MinFP"`.
#' @param max_fpr Used only for `"MinFP"`.
#' @return List with `y_pred_adjusted` and `threshold`.
#' @keywords internal
get_threshold <- function(y_test, y_prob,
                            method = c("ROC", "F1", "J", "MinFP"),
                            max_fpr = 0.01) {
  method <- match.arg(method)
  y_test <- as.integer(y_test); y_prob <- as.numeric(y_prob)

  if (method %in% c("ROC", "J")) {
    roc <- pROC::roc(response = y_test, predictor = y_prob,
                       levels = c(0L, 1L), direction = "<", quiet = TRUE)
    coords <- pROC::coords(roc, "all", ret = c("threshold", "tpr", "fpr"))
    youden <- coords$tpr - coords$fpr
    threshold <- coords$threshold[which.max(youden)]
  } else if (method == "F1") {
    grid <- sort(unique(c(0, y_prob, 1)))
    f1_vals <- vapply(grid, function(thr) {
      pred <- as.integer(y_prob > thr)
      tp <- sum(pred == 1L & y_test == 1L)
      fp <- sum(pred == 1L & y_test == 0L)
      fn <- sum(pred == 0L & y_test == 1L)
      precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
      recall    <- if ((tp + fn) > 0) tp / (tp + fn) else 0
      if (precision + recall == 0) return(0)
      2 * precision * recall / (precision + recall)
    }, numeric(1L))
    threshold <- grid[which.max(f1_vals)]
  } else if (method == "MinFP") {
    roc <- pROC::roc(response = y_test, predictor = y_prob,
                       levels = c(0L, 1L), direction = "<", quiet = TRUE)
    coords <- pROC::coords(roc, "all", ret = c("threshold", "fpr"))
    valid <- coords[coords$fpr <= max_fpr, , drop = FALSE]
    threshold <- if (nrow(valid) > 0L) {
      utils::tail(valid$threshold, 1L)
    } else {
      coords$threshold[1L]
    }
  }
  list(y_pred_adjusted = as.integer(y_prob > threshold),
        threshold = threshold)
}

#' Pairwise balanced random-forest classifier
#'
#' Builds a pairwise design matrix (feature-wise differences of every
#' sample pair, optionally augmented with the Euclidean distance), labels
#' each pair 1 if the two samples share a patient ID, then trains a class-
#' balanced random forest. The classifier backend is selectable.
#'
#' @param X_train,X_test Sample x feature matrices.
#' @param y_train,y_test Patient labels (vectors with one entry per row).
#' @param df_pivot_test Wide test frame including `Sample_ID` column --
#'   used to label misclassified pairs by sample.
#' @param compute_euclid Add a NaN-aware Euclidean distance feature.
#' @param method Threshold selection (see [get_threshold()]).
#' @param classifier_backend `"randomForest"` (default -- closest behaviour
#'   to Python's `imblearn.BalancedRandomForestClassifier` via per-tree
#'   balanced bootstrap), `"ranger"` (faster; class-weighted impurity),
#'   or `"themis_smote"` (SMOTE oversampling). See
#'   <https://github.com/fhradilak/spqrp_r/blob/main/articles/numerical-divergence.md>
#'   for the tradeoffs. Importance
#'   values returned in the results are normalised to sum to 1.0 across
#'   features (matching sklearn's `clf.feature_importances_` convention)
#'   regardless of backend.
#' @param k Fold number for diagnostic printing.
#' @param plots_per_sample Per-sample probability plots.
#' @param sample_decision_curve If `TRUE`, draw ROC + PR + threshold plots.
#' @param absolute Take absolute value of feature differences before
#'   passing to the model. (Stored after training is complete.)
#' @param quiet If `TRUE` (default), suppress informational status messages.
#'   Set `FALSE` to print progress and per-call summaries (sample counts,
#'   chosen cutoff, etc.). Warnings about genuine data issues -- e.g.
#'   samples dropped from the analysis -- are emitted regardless.
#' @return Named list as described in the package docs.
#' @examples
#' \donttest{
#' df <- spqrp_example_data("input_cohort_df")
#' # In practice, call the high-level [train_with_normalise()] instead --
#' # it handles the train/test split, normalisation, and pivoting for you.:
#' res <- train_with_normalise(df, plate_corrected = FALSE,
#'                              outlier_removal = FALSE)
#' res$classifier_backend
#' }
#' @export
train_pairwise_balanced_rand_forest <- function(X_train, y_train,
                                                  X_test, y_test,
                                                  df_pivot_test,
                                                  compute_euclid = TRUE,
                                                  method = "F1",
                                                  classifier_backend = c("randomForest", "ranger", "themis_smote"),
                                                  k = 0L,
                                                  plots_per_sample = FALSE,
                                                  sample_decision_curve = FALSE,
                                                  absolute = FALSE,
                                                  quiet = TRUE) {
  classifier_backend <- match.arg(classifier_backend)

  results_dict <- list(
    feature_importances = list(),
    feature_importances_test = list(),
    misclassified_pairs = list()
  )

  train_pairs <- create_pairs(X_train, y_train)
  test_pairs  <- create_pairs(X_test,  y_test)

  original_feature_names <- colnames(X_train)
  feature_names <- paste0("diff_", original_feature_names)
  if (compute_euclid) feature_names <- c(feature_names, "euclidean_distance")

  build_design <- function(pairwise, feature_names, compute_euclid) {
    if (length(pairwise) == 0L) {
      return(matrix(numeric(0L), ncol = length(feature_names),
                    dimnames = list(NULL, feature_names)))
    }
    rows <- lapply(pairwise, function(pair) {
      feat_diff <- compute_feature_differences(pair[[1L]], pair[[2L]])
      if (compute_euclid) {
        c(feat_diff,
          compute_nan_euclidean_distance(pair[[1L]], pair[[2L]]))
      } else {
        feat_diff
      }
    })
    mat <- do.call(rbind, rows)
    colnames(mat) <- feature_names
    mat
  }

  X_train_pairs <- build_design(train_pairs$pairwise, feature_names, compute_euclid)
  X_test_pairs  <- build_design(test_pairs$pairwise,  feature_names, compute_euclid)
  y_train_pairs <- train_pairs$labels
  y_test_pairs  <- test_pairs$labels

  # NaN/NA/Inf -> 0 to keep classifiers happy (randomForest, RANN, etc.)
  X_train_pairs[!is.finite(X_train_pairs)] <- 0
  X_test_pairs[!is.finite(X_test_pairs)]   <- 0

  X_train_df <- as.data.frame(X_train_pairs, check.names = FALSE)
  X_test_df  <- as.data.frame(X_test_pairs,  check.names = FALSE)
  y_train_fac <- factor(y_train_pairs, levels = c(0L, 1L))

  clf <- switch(
    classifier_backend,
    "ranger" = {
      class_counts <- table(y_train_fac)
      cw <- if (any(class_counts == 0L)) c(1, 1) else as.numeric(1 / class_counts)
      cw <- cw / sum(cw)
      ranger::ranger(
        x = X_train_df, y = y_train_fac,
        num.trees = 100L, probability = TRUE,
        class.weights = cw, seed = 42L,
        importance = "impurity"
      )
    },
    "randomForest" = {
      min_n <- min(table(y_train_fac))
      randomForest::randomForest(
        x = X_train_df, y = y_train_fac,
        ntree = 100L,
        sampsize = c(min_n, min_n),
        strata = y_train_fac,
        importance = TRUE
      )
    },
    "themis_smote" = {
      rlang::check_installed(c("themis", "recipes"))
      # `recipes` rejects column names starting with "." in formulas
      # (treats them as in-line functions), so use a plain name.
      train_df <- X_train_df
      train_df$spqrp_target <- y_train_fac
      rec <- recipes::recipe(spqrp_target ~ ., data = train_df) |>
        themis::step_smote(recipes::all_outcomes(), seed = 42L) |>
        recipes::prep()
      baked <- recipes::bake(rec, new_data = NULL)
      ranger::ranger(
        x = baked[, setdiff(names(baked), "spqrp_target"), drop = FALSE],
        y = baked$spqrp_target,
        num.trees = 100L, probability = TRUE, seed = 42L,
        importance = "impurity"
      )
    }
  )

  y_prob <- switch(
    classifier_backend,
    "ranger" = stats::predict(clf, data = X_test_df)$predictions[, "1"],
    "randomForest" = stats::predict(clf, newdata = X_test_df, type = "prob")[, "1"],
    "themis_smote" = stats::predict(clf, data = X_test_df)$predictions[, "1"]
  )

  thr <- get_threshold(y_test_pairs, y_prob, method = method)
  y_pred_adjusted <- thr$y_pred_adjusted
  optimal_threshold <- thr$threshold

  k <- k + 1L
  if (!quiet) cli::cli_h2("Fold {k}")
  evaluate_model(y_test_pairs, y_pred_adjusted, y_prob,
                   sample_decision_curve = sample_decision_curve,
                   optimal_threshold = optimal_threshold,
                   quiet = quiet)

  if (absolute) {
    X_train_pairs <- abs(X_train_pairs)
    X_test_pairs  <- abs(X_test_pairs)
  }

  if (plots_per_sample) {
    for (i in seq_len(nrow(df_pivot_test))) {
      plot_pairwise_probabilities_for_sample(
        sample_index = i,
        y_prob = y_prob,
        pair_indices_test = test_pairs$pair_indices,
        df_pivot = df_pivot_test,
        optimal_threshold = optimal_threshold,
        quiet = quiet
      )
    }
  }

  misclassified_indices <- which(y_pred_adjusted != y_test_pairs)
  results_dict <- get_misclassified_samples(
    df_pivot_test, misclassified_indices,
    test_pairs$pair_indices,
    y_pairs = y_test_pairs,
    y_pred_adjusted = y_pred_adjusted,
    results_dict = results_dict,
    quiet = quiet
  )

  importances <- switch(
    classifier_backend,
    "ranger" = clf$variable.importance,
    "randomForest" = clf$importance[, "MeanDecreaseGini"],
    "themis_smote" = clf$variable.importance
  )
  # Normalise to sum to 1.0 so the values match sklearn's
  # `clf.feature_importances_` semantics. Both `ranger$variable.importance`
  # and `randomForest::importance(type = 2)` return raw impurity decreases;
  # sklearn normalises across features. Rank order is preserved.
  total <- sum(importances, na.rm = TRUE)
  if (is.finite(total) && total > 0) {
    importances <- importances / total
  }
  results_dict <- print_and_add_feature_importance(
    importances = importances,
    feature_names = names(importances),
    results_dict = results_dict,
    quiet = quiet
  )

  train_pairs_names <- if (nrow(train_pairs$pair_indices) > 0L) {
    lapply(seq_len(nrow(train_pairs$pair_indices)), function(i) {
      c(rownames(X_train)[train_pairs$pair_indices[i, 1L]],
        rownames(X_train)[train_pairs$pair_indices[i, 2L]])
    })
  } else list()
  test_pairs_names <- if (nrow(test_pairs$pair_indices) > 0L) {
    lapply(seq_len(nrow(test_pairs$pair_indices)), function(i) {
      c(rownames(X_test)[test_pairs$pair_indices[i, 1L]],
        rownames(X_test)[test_pairs$pair_indices[i, 2L]])
    })
  } else list()

  out <- list(
    results_dict      = results_dict,
    train_pairs_names = train_pairs_names,
    test_pairs_names  = test_pairs_names,
    clf               = clf,
    classifier_backend = classifier_backend,
    X_train_pairs     = tibble::as_tibble(X_train_df)
  )
  class(out) <- c("spqrp_train", class(out))
  if (quiet) invisible(out) else out
}

#' Print a one-line summary of an spqrp_train object
#'
#' Displays the classifier backend, the number of training/test pairs,
#' and the feature count for the pairwise random-forest model returned
#' by [train_with_normalise()] and
#' [train_pairwise_balanced_rand_forest()].
#'
#' @param x A `spqrp_train` object.
#' @param ... Unused; present for S3 generic compatibility.
#' @return `x`, invisibly.
#' @export
print.spqrp_train <- function(x, ...) {
  cli::cli_h2("spqrp pairwise classifier ({x$classifier_backend})")
  cli::cli_text("Pairs (train): {length(x$train_pairs_names)}")
  cli::cli_text("Pairs (test): {length(x$test_pairs_names)}")
  cli::cli_text("Features: {ncol(x$X_train_pairs)}")
  invisible(x)
}

# Misclassified-pair bookkeeping (used by train_pairwise_balanced_rand_forest)
get_misclassified_samples <- function(df_pivot, misclassified_indices,
                                       pair_indices,
                                       y_pairs, y_pred_adjusted,
                                       results_dict,
                                       quiet = TRUE) {
  seen_pairs <- character()
  misclassified <- list()
  for (idx in misclassified_indices) {
    i <- pair_indices[idx, 1L]
    j <- pair_indices[idx, 2L]
    name_1 <- df_pivot$Sample_ID[i]
    name_2 <- df_pivot$Sample_ID[j]
    key <- paste(sort(c(name_1, name_2)), collapse = "||")
    if (key %in% seen_pairs) next
    seen_pairs <- c(seen_pairs, key)
    misclassified[[length(misclassified) + 1L]] <- list(
      `Sample 1` = name_1,
      `Sample 2` = name_2,
      `True Label` = y_pairs[idx],
      `Predicted Label` = y_pred_adjusted[idx]
    )
  }
  results_dict$misclassified_pairs <- misclassified
  if (!quiet) {
    for (p in misclassified) {
      cli::cli_text(
        "Sample 1: {p$`Sample 1`}, Sample 2: {p$`Sample 2`}, True:{p$`True Label`}, Predicted:{p$`Predicted Label`}"
      )
    }
  }
  results_dict
}

print_and_add_feature_importance <- function(importances, feature_names,
                                              results_dict, test = FALSE,
                                              quiet = TRUE) {
  ord <- order(-importances)
  sorted <- importances[ord]
  names_sorted <- feature_names[ord]
  key <- if (test) "feature_importances_test" else "feature_importances"
  n_show <- min(10L, length(sorted))
  if (!quiet) {
    cli::cli_h3("Top {n_show} important features")
    for (i in seq_len(n_show)) {
      cli::cli_text("{names_sorted[i]}: {round(sorted[i], 4)}")
    }
  }
  results_dict[[key]] <- lapply(seq_along(sorted), function(i) {
    list(feature = names_sorted[i], importance = unname(sorted[i]))
  })
  results_dict
}

# Plot eval curves (ROC / PR / threshold-vs-accuracy / probability split)
plot_evaluation <- function(y_test, y_probs, optimal_threshold, quiet = TRUE) {
  roc <- pROC::roc(y_test, y_probs, levels = c(0L, 1L), direction = "<",
                     quiet = TRUE)
  coords <- pROC::coords(roc, "all", ret = c("threshold", "tpr", "fpr"))
  auc_val <- as.numeric(pROC::auc(roc))

  roc_plot <- ggplot2::ggplot(coords,
                                ggplot2::aes(x = .data$fpr, y = .data$tpr)) +
    ggplot2::geom_line(colour = "blue") +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "dashed") +
    ggplot2::labs(
      x = "False Positive Rate", y = "True Positive Rate",
      title = sprintf("ROC Curve (AUC = %.2f)", auc_val)
    ) +
    ggplot2::theme_minimal()
  if (!quiet) print(roc_plot)

  # Precision-Recall curve via threshold sweep
  thresholds <- sort(unique(c(0, y_probs, 1)))
  pr_df <- do.call(rbind, lapply(thresholds, function(thr) {
    pred <- as.integer(y_probs > thr)
    tp <- sum(pred == 1L & y_test == 1L)
    fp <- sum(pred == 1L & y_test == 0L)
    fn <- sum(pred == 0L & y_test == 1L)
    prec <- if ((tp + fp) > 0) tp / (tp + fp) else 1
    rec  <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    data.frame(threshold = thr, precision = prec, recall = rec)
  }))
  pr_plot <- ggplot2::ggplot(pr_df,
                               ggplot2::aes(x = .data$recall, y = .data$precision)) +
    ggplot2::geom_line(colour = "blue") +
    ggplot2::labs(title = "Precision-Recall Curve",
                  x = "Recall", y = "Precision") +
    ggplot2::theme_minimal()
  if (!quiet) print(pr_plot)

  thr_plot <- ggplot2::ggplot(
    data = data.frame(prob = y_probs, label = factor(y_test)),
    mapping = ggplot2::aes(x = .data$prob, fill = .data$label)
  ) +
    ggplot2::geom_histogram(bins = 30L, alpha = 0.6, position = "identity") +
    ggplot2::geom_vline(xintercept = optimal_threshold, colour = "black",
                          linetype = "dashed") +
    ggplot2::scale_fill_manual(values = c("0" = "blue", "1" = "red")) +
    ggplot2::labs(title = "Predicted probability distribution",
                  x = "Predicted probability", y = "Count") +
    ggplot2::theme_minimal()
  if (!quiet) print(thr_plot)
}

evaluate_model <- function(y_test_pairs, y_pred_adjusted, y_prob,
                            sample_decision_curve = TRUE,
                            optimal_threshold = 0.5,
                            quiet = TRUE) {
  auc_val <- as.numeric(pROC::auc(pROC::roc(
    y_test_pairs, y_prob, levels = c(0L, 1L), direction = "<", quiet = TRUE
  )))
  accuracy <- mean(y_pred_adjusted == y_test_pairs)
  TP <- sum(y_pred_adjusted == 1L & y_test_pairs == 1L)
  TN <- sum(y_pred_adjusted == 0L & y_test_pairs == 0L)
  FP <- sum(y_pred_adjusted == 1L & y_test_pairs == 0L)
  FN <- sum(y_pred_adjusted == 0L & y_test_pairs == 1L)
  if (!quiet) {
    cli::cli_text("AUC-ROC: {round(auc_val, 4)}")
    cli::cli_text("Accuracy: {round(accuracy, 4)}")
    cli::cli_text("TP: {TP}, FP: {FP}, FN: {FN}, TN: {TN}")
  }
  if (sample_decision_curve) {
    plot_evaluation(y_test_pairs, y_prob, optimal_threshold, quiet = quiet)
  }
  invisible(NULL)
}

plot_pairwise_probabilities_for_sample <- function(sample_index, y_prob,
                                                     pair_indices_test, df_pivot,
                                                     optimal_threshold,
                                                     quiet = TRUE) {
  sample_name <- df_pivot$Sample_ID[sample_index]
  rows <- which(pair_indices_test[, 1L] == sample_index |
                  pair_indices_test[, 2L] == sample_index)
  if (length(rows) == 0L) return(invisible(NULL))
  other_idx <- ifelse(pair_indices_test[rows, 1L] == sample_index,
                       pair_indices_test[rows, 2L],
                       pair_indices_test[rows, 1L])
  plot_df <- data.frame(
    other = df_pivot$Sample_ID[other_idx],
    prob = y_prob[rows]
  )
  p <- ggplot2::ggplot(plot_df,
                        ggplot2::aes(x = .data$other, y = .data$prob)) +
    ggplot2::geom_point(colour = "blue", alpha = 0.7) +
    ggplot2::geom_hline(yintercept = optimal_threshold,
                         colour = "red", linetype = "dashed") +
    ggplot2::labs(title = paste0("Pairwise probabilities for sample: ", sample_name),
                  x = "Compared sample", y = "P(class = 1)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, size = 7))
  if (!quiet) print(p)
}

#' Convert classifier output to a Protein / Importance ranking
#'
#' Strips the `diff_` prefix the pairwise model adds to feature names.
#' Importance values are normalised to sum to 1.0 across features at
#' training time (matching sklearn's `clf.feature_importances_`
#' convention), so the numbers in the returned tibble are directly
#' comparable to Python output. Rank order is preserved across the
#' normalisation.
#'
#' @param results Output of [train_with_normalise()].
#' @return Tibble with `Protein` and `Importance` columns; `Importance`
#'   sums to ~1.0.
#' @examples
#' \donttest{
#' df <- spqrp_example_data("input_cohort_df")
#' results <- train_with_normalise(df, plate_corrected = FALSE,
#'                                  outlier_removal = FALSE)
#' retrieve_ranking(results)
#' }
#' @export
retrieve_ranking <- function(results) {
  fi <- results$results_dict$feature_importances
  if (length(fi) == 0L) return(tibble::tibble(Protein = character(), Importance = numeric()))
  df <- do.call(rbind, lapply(fi, function(x) {
    data.frame(feature = x$feature, importance = x$importance,
               stringsAsFactors = FALSE)
  }))
  feats <- df$feature
  feats <- ifelse(startsWith(feats, "diff_"), substring(feats, 6L), feats)
  tibble::tibble(Protein = feats, Importance = df$importance)
}

#' End-to-end ranking pipeline: filter, normalise, optionally plate-correct, train RF
#'
#' Mirrors `protein_selection.train_with_normalise` from the Python package
#' but exposes `classifier_backend` so users can compare three RF variants
#' (`"ranger"`, `"randomForest"`, `"themis_smote"`). See
#' <https://github.com/fhradilak/spqrp_r/blob/main/articles/numerical-divergence.md>
#' for the tradeoffs.
#'
#' @param df Long-format cohort data frame.
#' @param threshold Occurrence-filter threshold.
#' @param test_size Patient-level test fraction.
#' @param plate_corrected If `TRUE`, run plate-effect residualisation.
#' @param individual Patient column.
#' @param sample Sample column.
#' @param compute_euclid Add NaN-aware Euclidean distance feature.
#' @param method Threshold-selection strategy.
#' @param outlier_removal Run [by_isolation_forest()] on each split.
#' @param train_individuals,test_individuals Explicit split overrides.
#' @param sample_decision_curve Draw ROC/PR curves.
#' @param classifier_backend `"randomForest"` (default -- closest behaviour
#'   to Python's `imblearn.BalancedRandomForestClassifier`), `"ranger"`
#'   (faster), or `"themis_smote"`. The default was changed from
#'   `"ranger"` to `"randomForest"` to bring R rankings closer to the
#'   Python port. See
#'   <https://github.com/fhradilak/spqrp_r/blob/main/articles/numerical-divergence.md>.
#' @param importance_method Unused placeholder (kept for API parity).
#' @param plot_per_sample Per-sample probability plots.
#' @param absolute Use absolute pairwise differences.
#' @param quiet If `TRUE` (default), suppress informational status messages
#'   (train/test split listing, "Proteins only in test set", outliers
#'   removed, fold headers, per-fold metrics, top-importance list, and
#'   per-misclassified-pair prints) and skip auto-rendering of the ROC
#'   / PR / probability plots. Set `FALSE` to print everything.
#'   Warnings about genuine data issues are emitted regardless.
#' @return `spqrp_train` S3 object (a named list with classifier, pair
#'   indices, feature importances, misclassified pairs).
#' @examples
#' \donttest{
#' df <- spqrp_example_data("input_cohort_df")
#' res <- train_with_normalise(df, plate_corrected = FALSE,
#'                              outlier_removal = FALSE)
#' retrieve_ranking(res)
#' }
#' @export
train_with_normalise <- function(df,
                                   threshold = 0.7,
                                   test_size = 0.3,
                                   plate_corrected = TRUE,
                                   individual = PATIENT,
                                   sample = SAMPLE,
                                   compute_euclid = FALSE,
                                   method = "F1",
                                   outlier_removal = TRUE,
                                   train_individuals = NULL,
                                   test_individuals = NULL,
                                   sample_decision_curve = FALSE,
                                   classifier_backend = c("randomForest", "ranger", "themis_smote"),
                                   importance_method = "impurity",
                                   plot_per_sample = FALSE,
                                   absolute = FALSE,
                                   quiet = TRUE) {
  classifier_backend <- match.arg(classifier_backend)
  check_input_data_format(df)

  # Reindex to full (Protein x Sample) grid filling missing intensities
  proteins <- unique(df$Protein)
  samples <- unique(df[[sample]])
  full <- tidyr::expand_grid(Protein = proteins,
                              !!sample := samples)
  df_complete <- dplyr::left_join(full, df, by = c("Protein", sample))

  df_clean <- filter_by_occurrence(df_complete, threshold)
  df_clean <- log_transform(df_clean)

  if (is.null(train_individuals) && is.null(test_individuals)) {
    unique_individuals <- unique(df_clean[[individual]])
    split_res <- split_by_patient(unique_individuals, test_size, seed = 42L)
    train_individuals <- split_res$train
    test_individuals  <- split_res$test
    if (!quiet) {
      cli::cli_text("train_individuals: {.val {train_individuals}}")
      cli::cli_text("test_individuals: {.val {test_individuals}}")
    }
  }

  df_train <- df_clean[df_clean[[individual]] %in% train_individuals, , drop = FALSE]
  df_test  <- df_clean[df_clean[[individual]] %in% test_individuals,  , drop = FALSE]

  extra_in_test <- setdiff(unique(df_test$Protein), unique(df_train$Protein))
  if (!quiet) {
    cli::cli_text("Proteins only in test set: {.val {extra_in_test}}")
  }

  norm_train <- normalize_medianintensity(df_train, plot = FALSE)
  norm_test  <- normalize_medianintensity(df_test,  plot = FALSE)
  df_train_norm <- norm_train$data
  df_test_norm  <- norm_test$data

  df_train_normalized <- df_train_norm
  df_test_normalized  <- df_test_norm
  if (plate_corrected) {
    df_train_normalized <- plate_correct_residuals_by_protein(df_train_normalized,
                                                                 individual = individual,
                                                                 sample = sample)
    df_test_normalized  <- plate_correct_residuals_by_protein(df_test_normalized,
                                                                 individual = individual,
                                                                 sample = sample)
  }

  df_train_pivot <- pivot_df(df_train_normalized,
                              index = c(individual, sample))
  df_test_pivot  <- pivot_df(df_test_normalized,
                              index = c(individual, sample))
  # Align columns
  protein_cols <- setdiff(names(df_train_pivot), c(individual, sample))
  missing_in_test <- setdiff(protein_cols, names(df_test_pivot))
  for (col in missing_in_test) df_test_pivot[[col]] <- NA_real_
  df_test_pivot <- df_test_pivot[, c(individual, sample, protein_cols), drop = FALSE]

  if (outlier_removal) {
    out_train <- by_isolation_forest(df_train_normalized, impute_median = TRUE,
                                       quiet = quiet)$outlier_list
    out_test  <- by_isolation_forest(df_test_normalized,  impute_median = TRUE,
                                       quiet = quiet)$outlier_list
    if (!quiet) {
      cli::cli_text("Outliers (train): {.val {out_train}}")
      cli::cli_text("Outliers (test): {.val {out_test}}")
    }
    df_train_clean <- df_train_pivot[!df_train_pivot[[sample]] %in% out_train, , drop = FALSE]
    df_test_clean  <- df_test_pivot[!df_test_pivot[[sample]]  %in% out_test,  , drop = FALSE]
  } else {
    df_train_clean <- df_train_pivot
    df_test_clean  <- df_test_pivot
  }

  X_train <- as.data.frame(df_train_clean[, protein_cols, drop = FALSE])
  rownames(X_train) <- df_train_clean[[sample]]
  X_train[is.na(X_train)] <- 0
  y_train <- df_train_clean[[individual]]

  X_test <- as.data.frame(df_test_clean[, protein_cols, drop = FALSE])
  rownames(X_test) <- df_test_clean[[sample]]
  X_test[is.na(X_test)] <- 0
  y_test <- df_test_clean[[individual]]

  train_pairwise_balanced_rand_forest(
    X_train = X_train, y_train = y_train,
    X_test  = X_test,  y_test  = y_test,
    df_pivot_test = df_test_clean,
    compute_euclid = compute_euclid,
    method = method,
    classifier_backend = classifier_backend,
    sample_decision_curve = sample_decision_curve,
    plots_per_sample = plot_per_sample,
    absolute = absolute,
    quiet = quiet
  )
}
