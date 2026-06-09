# Distance computation, evaluation metrics, optimization comparators.
# Port of spqrp/spqrp/helpers.py.

#' Percentile cutoff (numpy.percentile-equivalent)
#'
#' Returns the `percentile`-th percentile of `distances` using linear
#' interpolation (`stats::quantile` type 7), matching `numpy.percentile`'s
#' default for the values used in the package.
#'
#' @param distances Numeric vector of pairwise distances.
#' @param percentile Percentile in `[0, 100]`.
#' @return A single numeric.
#' @keywords internal
percentile_cutoff <- function(distances, percentile = 25) {
  unname(stats::quantile(
    distances, probs = percentile / 100,
    type = 7, na.rm = TRUE, names = FALSE
  ))
}

#' Pairwise distance matrix on a long-format intensity table
#'
#' Pivots the long-format `df_dist` to a (sample x protein) wide matrix
#' (missing values filled with 0) and computes pairwise distances using
#' the requested metric.
#'
#' @param df_dist Long-format data frame with `Sample_ID`, `Protein`, and
#'   the `intensity` column.
#' @param metric One of `"correlation"`, `"euclidean"`, `"manhattan"`,
#'   `"minkowski"`, or `"fractional"` (Minkowski with `p = fractional_p`).
#' @param intensity Name of the intensity column. Defaults to `"Intensity"`.
#' @param fractional_p Exponent for the fractional / Minkowski metric.
#' @param index Column used as the sample identifier. Defaults to
#'   `"Sample_ID"`.
#' @return A list with `distance_matrix` (numeric matrix with row/col names
#'   set to sample IDs) and `df_pivot` (the wide tibble used to compute it).
#' @keywords internal
get_distances <- function(df_dist,
                          metric = "correlation",
                          intensity = "Intensity",
                          fractional_p = 0.5,
                          index = SAMPLE) {
  # Build wide pivot, averaging duplicate (sample, protein) cells (matches
  # Python's pivot_table(..., aggfunc="mean")).
  df_pivot <- tidyr::pivot_wider(
    df_dist,
    id_cols = dplyr::all_of(index),
    names_from = "Protein",
    values_from = dplyr::all_of(intensity),
    values_fn = function(x) mean(x, na.rm = TRUE),
    values_fill = 0
  )

  # Replace NaN (from all-NA groups) with 0 to match Python's fillna(0)
  mat_full <- as.matrix(df_pivot[, setdiff(names(df_pivot), index), drop = FALSE])
  mat_full[is.nan(mat_full) | is.na(mat_full)] <- 0
  rownames(mat_full) <- as.character(df_pivot[[index]])

  distance_matrix <- pairwise_distance_matrix(mat_full, metric, fractional_p)

  list(distance_matrix = distance_matrix, df_pivot = df_pivot)
}

# Internal distance dispatcher: takes a numeric matrix (samples in rows)
# and returns a symmetric distance matrix with row/col names preserved.
pairwise_distance_matrix <- function(mat, metric, fractional_p = 0.5) {
  n <- nrow(mat)
  if (metric == "correlation") {
    # 1 - Pearson correlation between rows
    suppressWarnings(corr <- stats::cor(t(mat)))
    corr[is.na(corr)] <- 0
    d <- 1 - corr
    d[d < 0] <- 0
  } else if (metric == "fractional") {
    d <- as.matrix(stats::dist(mat, method = "minkowski", p = fractional_p))
  } else if (metric %in% c("euclidean", "manhattan", "minkowski")) {
    d <- as.matrix(stats::dist(mat, method = metric, p = fractional_p))
  } else {
    cli::cli_abort(c(
      "Unsupported metric {.val {metric}}.",
      "i" = "Supported: {.val {c('correlation','euclidean','manhattan','minkowski','fractional')}}."
    ))
  }
  diag(d) <- 0
  rownames(d) <- rownames(mat)
  colnames(d) <- rownames(mat)
  d
}

#' k nearest neighbours of every sample
#'
#' @param df_pivot Wide tibble (rows = samples, columns = proteins).
#' @param distance_matrix Symmetric distance matrix in the same row order as
#'   `df_pivot`.
#' @param k Number of neighbours per sample.
#' @return Tibble with columns `Neighbor_1..k` and `Distance_1..k`,
#'   row-aligned with `df_pivot`.
#' @keywords internal
get_nearest_neighbours <- function(df_pivot, distance_matrix, k = 4) {
  n_samples <- nrow(df_pivot)
  k <- min(k, n_samples - 1L)
  if (k < 1L) {
    return(tibble::tibble(.rows = n_samples))
  }
  sample_ids <- rownames(distance_matrix) %||% as.character(df_pivot[[SAMPLE]])

  neighbour_names <- matrix(NA_character_, n_samples, k)
  neighbour_dists <- matrix(NA_real_, n_samples, k)
  for (i in seq_len(n_samples)) {
    ord <- order(distance_matrix[i, ])
    # First entry is self (distance 0)
    keep <- ord[2:(k + 1L)]
    neighbour_names[i, ] <- sample_ids[keep]
    neighbour_dists[i, ] <- round(distance_matrix[i, keep], 2)
  }

  out <- tibble::tibble(!!SAMPLE := sample_ids)
  for (j in seq_len(k)) {
    out[[paste0("Neighbor_", j)]]  <- neighbour_names[, j]
    out[[paste0("Distance_", j)]]  <- neighbour_dists[, j]
  }
  out
}

#' Split pairwise distances into belonging / not-belonging by cutoff
#'
#' Pairs with `distance <= cutoff` go into `belonging`, others into
#' `not_belonging`. Patient IDs are looked up via `sample_patient_mapping`.
#'
#' @param distance_matrix Symmetric numeric matrix.
#' @param cutoff Numeric cutoff.
#' @param sample_patient_mapping Named character vector (names = `Sample_ID`,
#'   values = `Patient_ID`).
#' @param sample_order Row/col order of `distance_matrix` as a character
#'   vector of sample IDs.
#' @param quiet If `TRUE` (default), suppress informational status messages.
#'   Set `FALSE` to print progress and per-call summaries (sample counts,
#'   chosen cutoff, etc.). Warnings about genuine data issues -- e.g.
#'   samples dropped from the analysis -- are emitted regardless.
#' @return List with `belonging` and `not_belonging` tibbles.
#' @keywords internal
get_sample_relations_by_cutoff <- function(distance_matrix,
                                            cutoff,
                                            sample_patient_mapping,
                                            sample_order,
                                            quiet = TRUE) {
  n <- length(sample_order)
  triu <- which(upper.tri(distance_matrix), arr.ind = TRUE)
  i <- triu[, 1L]
  j <- triu[, 2L]

  s1 <- sample_order[i]
  s2 <- sample_order[j]
  d  <- distance_matrix[cbind(i, j)]

  df <- tibble::tibble(
    sample1      = s1,
    sample2      = s2,
    patient_id_1 = unname(sample_patient_mapping[s1]),
    patient_id_2 = unname(sample_patient_mapping[s2]),
    distance     = d
  )

  belonging     <- df[df$distance <= cutoff, , drop = FALSE]
  not_belonging <- df[df$distance >  cutoff, , drop = FALSE]
  rownames(belonging) <- NULL
  rownames(not_belonging) <- NULL

  if (!quiet) {
    cli::cli_h2("Belonging pairs")
    print(belonging)
    cli::cli_h2("Not-belonging pairs")
    print(not_belonging)
  }

  list(belonging = belonging, not_belonging = not_belonging)
}

#' Evaluate pairwise belonging/not-belonging classification
#'
#' Computes TP/FP/FN/TN, precision, recall (sensitivity), F1, accuracy,
#' balanced accuracy from the result of [get_sample_relations_by_cutoff()].
#'
#' @param belonging Tibble of pairs flagged as belonging (same cluster).
#' @param not_belonging Tibble of pairs flagged as not belonging.
#' @param quiet If `TRUE` (default), suppress informational status messages.
#'   Set `FALSE` to print progress and per-call summaries (sample counts,
#'   chosen cutoff, etc.). Warnings about genuine data issues -- e.g.
#'   samples dropped from the analysis -- are emitted regardless.
#' @return Named list of metrics; preserves the keys used by the Python
#'   port (e.g. `False_Negative_Pairs`, `False_Positive_Distances`).
#' @keywords internal
get_evaluation_metrics <- function(belonging, not_belonging, quiet = TRUE) {
  tp_mask <- belonging$patient_id_1 == belonging$patient_id_2
  true_positives        <- Map(c, belonging$sample1[tp_mask],  belonging$sample2[tp_mask])
  true_positive_dists   <- belonging$distance[tp_mask]

  fp_mask <- !tp_mask
  false_positives       <- Map(c, belonging$sample1[fp_mask],  belonging$sample2[fp_mask])
  false_positive_dists  <- belonging$distance[fp_mask]

  tn_mask <- not_belonging$patient_id_1 != not_belonging$patient_id_2
  true_negatives        <- Map(c, not_belonging$sample1[tn_mask],  not_belonging$sample2[tn_mask])
  true_negative_dists   <- not_belonging$distance[tn_mask]

  fn_mask <- !tn_mask
  false_negatives       <- Map(c, not_belonging$sample1[fn_mask],  not_belonging$sample2[fn_mask])
  false_negative_dists  <- not_belonging$distance[fn_mask]

  TP <- length(true_positives)
  FP <- length(false_positives)
  TN <- length(true_negatives)
  FN <- length(false_negatives)
  n_all <- TP + FP + TN + FN

  accuracy <- if (n_all > 0) (TP + TN) / n_all else 0
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else 0
  sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else 0
  balanced_accuracy <- if ((TP + FN) > 0 && (TN + FP) > 0) {
    0.5 * (TP / (TP + FN) + TN / (TN + FP))
  } else {
    0
  }
  f1 <- if ((precision + sensitivity) > 0) {
    2 * precision * sensitivity / (precision + sensitivity)
  } else {
    0
  }

  if (!quiet) {
    cli::cli_h2("Pairwise classification metrics: Threshold-based")
    cli::cli_text("FP + FN: {FP + FN}")
    cli::cli_text("TP: {TP}, FP: {FP}, FN: {FN}, TN: {TN}")
    cli::cli_text("Accuracy: {round(accuracy, 4)}")
    cli::cli_text("Balanced accuracy: {round(balanced_accuracy, 4)}")
    cli::cli_text("Precision: {round(precision, 4)}")
    cli::cli_text("Sensitivity: {round(sensitivity, 4)}")
    cli::cli_text("F1: {round(f1, 4)}")
  }

  list(
    TP = TP, FP = FP, FN = FN, TN = TN,
    Accuracy          = accuracy,
    Balanced_Accuracy = balanced_accuracy,
    Precision         = precision,
    Sensitivity       = sensitivity,
    F1                = f1,
    False_Negative_Pairs     = false_negatives,
    False_Negative_Distances = false_negative_dists,
    False_Positive_Pairs     = false_positives,
    False_Positive_Distances = false_positive_dists,
    True_Negative_Pairs      = true_negatives,
    True_Negative_Distances  = true_negative_dists,
    True_Positive_Pairs      = true_positives,
    True_Positive_Distances  = true_positive_dists
  )
}

# Comparator used by optimize_parameters() to decide whether a new result
# replaces the best-so-far given the chosen strategy.
is_better_result <- function(result, strategy, thresholds) {
  fp <- result$FP
  fn <- result$FN
  precision <- result$Precision
  sensitivity <- result$Sensitivity
  f1 <- result$F1 %||% 0

  switch(
    strategy,
    "fp+fn" = {
      curr_sum <- fp + fn
      curr_sum < thresholds[["fp+fn"]] ||
        (curr_sum == thresholds[["fp+fn"]] && precision > thresholds[["f1"]])
    },
    "fp" = {
      fp < thresholds[["fp"]] ||
        (fp == thresholds[["fp"]] && (fp + fn) < thresholds[["fp+fn"]])
    },
    "fn" = {
      fn < thresholds[["fn"]] ||
        (fn == thresholds[["fn"]] && (fp + fn) < thresholds[["fp+fn"]])
    },
    "F1"          = f1 > thresholds[["f1"]],
    "sensitivity" = sensitivity > thresholds[["sensitivity"]],
    "precision"   = precision > thresholds[["precision"]],
    FALSE
  )
}
