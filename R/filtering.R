# Isolation-forest outlier detection (port of spqrp/spqrp/filtering.py).

#' Long-to-wide pivot (samples in rows, proteins in columns)
#'
#' Pivots a long-format intensity table into a (sample x protein) wide
#' tibble. Rows are sorted by `Sample_ID` and columns by protein name
#' (both via codepoint / radix sort) so the matrix layout matches
#' pandas' `pivot_table` output in the Python port -- a prerequisite for
#' reproducible IsolationForest outputs across the two implementations.
#'
#' @param intensity_df Long-format data frame with `Sample_ID`, `Protein`,
#'   and intensity values.
#' @param value_name Optional name of the intensity column.
#' @return Wide tibble.
#' @keywords internal
long_to_wide <- function(intensity_df, value_name = NULL) {
  values_name <- value_name %||% "Intensity"
  out <- tidyr::pivot_wider(
    intensity_df,
    id_cols = dplyr::all_of(SAMPLE),
    names_from = "Protein",
    values_from = dplyr::all_of(values_name)
  )
  # Row order seeds the isolation-forest subsampling RNG; column order
  # seeds feature-index choices. Force both to the same codepoint sort
  # pandas uses so subsamples line up across ports.
  out <- out[order(out[[SAMPLE]], method = "radix"), , drop = FALSE]
  non_id_cols <- setdiff(names(out), SAMPLE)
  out[, c(SAMPLE, sort(non_id_cols, method = "radix")), drop = FALSE]
}

#' Isolation Forest outlier detection
#'
#' Pivots to a (sample x protein) matrix, runs an Isolation Forest via
#' the `solitude` package (a pure-R port of the same Liu et al. 2008
#' algorithm that scikit-learn's `IsolationForest` uses), and returns
#' the data frame with outlier rows removed plus a tibble of per-sample
#' anomaly scores.
#'
#' Two ways to decide which samples are outliers, mirroring sklearn's
#' `IsolationForest` API:
#'
#' * `contamination = "auto"` (default) -- flag every sample whose
#'   anomaly score exceeds `outlier_threshold`.
#' * `contamination` set to a numeric in `[0, 1]` -- flag exactly the top
#'   `contamination * 100`% of samples by anomaly score, ignoring
#'   `outlier_threshold`. Mirrors sklearn's
#'   `IsolationForest(contamination = 0.1)`.
#'
#' On the **sklearn score scale**, `contamination = "auto"` corresponds
#' to a threshold of `0.5`. solitude's scores, however, are
#' systematically shifted upward because `solitude` (via `ranger`) uses
#' `mtry = ncol - 1` and `extratrees` split bounds drawn from the full
#' dataset rather than from the per-tree subsample. The result is that
#' inlier scores typically sit between `0.55` and `0.60` even on clean
#' data, so the sklearn-calibrated `0.5` cutoff would flag everything.
#' The default `outlier_threshold = 0.6` below is calibrated empirically
#' for solitude's distribution and reproduces sklearn's
#' "few-to-zero outliers on clean data" behaviour. Lower it
#' (e.g. `0.55`) for more aggressive flagging, or use `contamination`
#' for a percentile-based rule.
#'
#' @param protein_df Long-format intensity data frame.
#' @param peptide_df Optional peptide-level data frame; subset alongside
#'   `protein_df` using the same outlier list.
#' @param n_estimators Number of trees.
#' @param impute_zero Replace NA intensities with 0 before fitting.
#' @param impute_median Replace NA intensities with column-wise median.
#' @param outlier_threshold Used when `contamination = "auto"`. Anomaly
#'   score above which a sample is flagged. Default `0.6` (calibrated
#'   for solitude's score scale; on sklearn's scale this would be `0.5`).
#' @param contamination Either `"auto"` (default; use `outlier_threshold`)
#'   or a numeric in `[0, 1]` specifying the fraction of the data to
#'   flag as outliers (top-by-score). Mirrors sklearn's
#'   `IsolationForest(contamination=...)` API.
#' @param quiet If `TRUE` (default), suppress informational status messages.
#'   Set `FALSE` to print progress and per-call summaries (sample counts,
#'   chosen cutoff, etc.). Warnings about genuine data issues -- e.g.
#'   samples dropped from the analysis -- are emitted regardless.
#' @return Invisibly returns a named list with `protein_df`, `peptide_df`,
#'   `outlier_list`, `anomaly_df`, and possibly `messages` on failure.
#'   `invisible()` keeps the REPL silent on unassigned calls; assign the
#'   result to a name and inspect with `result$protein_df` etc.
#' @examples
#' \donttest{
#' df <- spqrp_example_data("input_cohort_df")
#' res <- by_isolation_forest(df, impute_median = TRUE)
#' res$outlier_list
#' }
#' @export
by_isolation_forest <- function(protein_df,
                                 peptide_df = NULL,
                                 n_estimators = 100L,
                                 impute_zero = FALSE,
                                 impute_median = FALSE,
                                 outlier_threshold = 0.6,
                                 contamination = "auto",
                                 quiet = TRUE) {

  if (!(identical(contamination, "auto") ||
        (is.numeric(contamination) && length(contamination) == 1L &&
         contamination >= 0 && contamination <= 1))) {
    cli::cli_abort(
      "{.arg contamination} must be {.val auto} or a number in [0, 1]."
    )
  }

  transformed_df <- long_to_wide(protein_df)
  sample_ids <- transformed_df[[SAMPLE]]
  mat <- as.matrix(transformed_df[, setdiff(names(transformed_df), SAMPLE), drop = FALSE])

  if (impute_zero) {
    mat[is.na(mat)] <- 0
  } else if (impute_median) {
    # Per-column median, matching pandas' fillna(df.median(skipna=True)).
    # All-NA columns stay NA -- sklearn errors on them, and we let the
    # anyNA check below produce the same diagnostic Python emits.
    for (j in seq_len(ncol(mat))) {
      col <- mat[, j]
      if (anyNA(col)) {
        col[is.na(col)] <- stats::median(col, na.rm = TRUE)
        mat[, j] <- col
      }
    }
  }

  if (anyNA(mat)) {
    msg <- paste(
      "Outlier detection by IsolationForest does not accept missing values.",
      "Set impute_zero=TRUE or impute_median=TRUE or pre-fill NAs."
    )
    return(invisible(list(
      protein_df = protein_df,
      peptide_df = peptide_df,
      outlier_list = NULL,
      anomaly_df = NULL,
      messages = list(list(level = "ERROR", msg = msg))
    )))
  }

  # Match Python sklearn IsolationForest parameters:
  #   max_samples  = len(transformed_df) // 2  -> sample_size
  #   n_estimators = 100                       -> num_trees
  #   n_jobs       = 1                         -> nproc (single-threaded
  #     is required for run-to-run reproducibility; Python comparisons
  #     must also set n_jobs=1 for parity).
  # Seed note: sklearn's RNG path is incompatible with ranger's anyway
  # (different libraries, different streams), so `random_state` and
  # `seed` cannot share a value. We use `1L` because ranger treats
  # `seed = 0` as "ignore the seed and pick randomly" -- that would make
  # solitude non-deterministic across calls.
  sample_size <- nrow(mat) %/% 2L
  mat_df <- as.data.frame(mat, stringsAsFactors = FALSE)
  iforest <- solitude::isolationForest$new(
    sample_size = sample_size,
    num_trees   = n_estimators,
    seed        = 1L,
    nproc       = 1L
  )
  # solitude logs INFO lines via `lgr::lgr` (the root lgr logger) on
  # fit/predict; raise its threshold to "warn" so they don't clutter
  # user output, then restore.
  old_threshold <- lgr::lgr$threshold
  lgr::lgr$set_threshold("warn")
  on.exit(lgr::lgr$set_threshold(old_threshold), add = TRUE)
  iforest$fit(mat_df)
  pred <- iforest$predict(mat_df)
  scores <- as.numeric(pred$anomaly_score)
  
  outlier_flag <- if (identical(contamination, "auto")) {
    scores > outlier_threshold
  } else if (contamination == 0) {
    rep(FALSE, length(scores))
  } else {
    # Flag the top `contamination * 100`% of samples by score.
    n_outliers <- max(1L, round(length(scores) * contamination))
    dynamic_cutoff <- sort(scores, decreasing = TRUE)[n_outliers]
    scores >= dynamic_cutoff
  }

  anomaly_df <- tibble::tibble(
    !!SAMPLE := sample_ids,
    `Anomaly Score` = scores,
    Outlier = outlier_flag
  )
  outlier_list <- sample_ids[anomaly_df$Outlier]

  if (!quiet) {
    n_flagged <- length(outlier_list)
    if (n_flagged > 0L) {
      cli::cli_inform(c(
        "Isolation Forest flagged {n_flagged} of {length(sample_ids)} sample{?s} as outlier{?s}:",
        "i" = "{.val {outlier_list}}"
      ))
    } else {
      cli::cli_inform(
        "Isolation Forest flagged no outliers (of {length(sample_ids)} samples)."
      )
    }
  }

  protein_df <- protein_df[!protein_df[[SAMPLE]] %in% outlier_list, , drop = FALSE]
  if (!is.null(peptide_df)) {
    peptide_df <- peptide_df[!peptide_df[[SAMPLE]] %in% outlier_list, , drop = FALSE]
  }

  invisible(list(
    protein_df = tibble::as_tibble(protein_df),
    peptide_df = peptide_df,
    outlier_list = outlier_list,
    anomaly_df = anomaly_df
  ))
}

PLOT_COLOR_SEQUENCE <- c("#4A536A", "#CE5A5A", "#87A8B9", "#8E3325", "#E2A46D")
PLOT_PRIMARY_COLOR   <- PLOT_COLOR_SEQUENCE[1L]
PLOT_SECONDARY_COLOR <- PLOT_COLOR_SEQUENCE[2L]

create_anomaly_score_bar_plot <- function(anomaly_df,
                                          colour_outlier = PLOT_SECONDARY_COLOR,
                                          colour_non_outlier = PLOT_PRIMARY_COLOR,
                                          title = "") {
  rlang::check_installed("plotly", reason = "to render the anomaly score plot")
  anomaly_df$Outlier <- as.character(anomaly_df$Outlier)
  plotly::plot_ly(
    data = anomaly_df,
    x = stats::as.formula(paste0("~`", SAMPLE, "`")),
    y = ~`Anomaly Score`,
    type = "bar",
    color = ~Outlier,
    colors = stats::setNames(
      c(colour_non_outlier, colour_outlier),
      c("FALSE", "TRUE")
    ),
    hoverinfo = "text",
    text = ~paste0(
      "Sample: ", .data[[SAMPLE]], "<br>",
      "Anomaly Score: ", round(`Anomaly Score`, 4), "<br>",
      "Outlier: ", Outlier
    )
  ) |>
    plotly::layout(
      title = title,
      xaxis = list(visible = FALSE, showticklabels = FALSE),
      yaxis = list(title = "Anomaly Score", visible = TRUE)
    )
}

#' Plot per-sample anomaly scores from the isolation forest
#'
#' @param output_anomaly_df Tibble returned in `anomaly_df` from
#'   [by_isolation_forest()].
#' @param title Plot title.
#' @return A `plotly` figure (printed when invoked at top level).
#' @keywords internal
by_isolation_forest_plot <- function(output_anomaly_df, title = "") {
  fig <- create_anomaly_score_bar_plot(output_anomaly_df, title = title)
  fig
}

#' Remove samples flagged as outliers by Isolation Forest
#'
#' Convenience wrapper around [by_isolation_forest()] with median
#' imputation. Removes samples (not proteins) whose intensity profile
#' looks anomalous compared to the rest of the cohort.
#'
#' Pass `contamination = 0.1` (or any fraction) to mimic sklearn's
#' `IsolationForest(contamination = 0.1)` behaviour, or keep the default
#' `contamination = "auto"` to use the conservative absolute threshold.
#'
#' The returned list includes `anomaly_plot`, a `plotly` bar chart of
#' per-sample anomaly scores coloured by outlier flag. Printing the
#' object at the R REPL (or `print(result$anomaly_plot)` inside a
#' script) renders the chart -- mirroring the Python wrapper's
#' auto-shown bar plot, but without surprising side effects when the
#' function is called non-interactively.
#'
#' @param df Long-format intensity data frame.
#' @param sample Sample column (defaults to `"Sample_ID"`).
#' @param contamination `"auto"` (default) or a numeric in `[0, 1]`.
#'   See [by_isolation_forest()] for details.
#' @param outlier_threshold Anomaly-score cutoff used when
#'   `contamination = "auto"`. Default `0.6`, calibrated empirically
#'   for solitude's anomaly-score distribution. See
#'   [by_isolation_forest()] for the rationale.
#' @param quiet If `TRUE` (default), suppress informational status messages.
#'   Set `FALSE` to print progress and per-call summaries (sample counts,
#'   chosen cutoff, etc.). Warnings about genuine data issues -- e.g.
#'   samples dropped from the analysis -- are emitted regardless.
#' @return Invisibly returns a named list with components:
#'   * `df` -- filtered tibble (same shape as `df`, fewer rows)
#'   * `anomaly_df` -- per-sample tibble of `Sample_ID`, `Anomaly Score`,
#'     `Outlier`
#'   * `outlier_list` -- character vector of flagged `Sample_ID`s
#'   * `anomaly_plot` -- a `plotly` figure; `print(result$anomaly_plot)`
#'     to view the bar chart. `NULL` if the optional `plotly` package is
#'     not installed (a message explains how to enable it).
#'
#'   The return is wrapped in `invisible()` so unassigned REPL calls stay
#'   silent (matches `quiet = TRUE`). Assign to a name to inspect.
#' @examples
#' \donttest{
#' df <- spqrp_example_data("input_cohort_df")
#' filtered <- remove_outlier_samples(df, contamination = "auto")
#' filtered$outlier_list
#' head(filtered$df)
#' }
#' @export
remove_outlier_samples <- function(df, sample = SAMPLE,
                                    contamination = "auto",
                                    outlier_threshold = 0.6,
                                    quiet = TRUE) {
  forest_dict <- by_isolation_forest(df, impute_median = TRUE,
                                       contamination = contamination,
                                       outlier_threshold = outlier_threshold,
                                       quiet = quiet)
  df_filtered <- df[!df[[sample]] %in% forest_dict$outlier_list, , drop = FALSE]
  # `plotly` is an optional dependency; only build the chart when it is
  # available so the core filtering workflow does not require it. Let the
  # user know how to get the plot if they expected it.
  if (rlang::is_installed("plotly")) {
    anomaly_plot <- by_isolation_forest_plot(forest_dict$anomaly_df)
  } else {
    anomaly_plot <- NULL
    cli::cli_inform(c(
      "!" = "{.field anomaly_plot} is {.code NULL} because the optional package
             {.pkg plotly} is not installed.",
      "i" = "Install it with {.run install.packages(\"plotly\")} to get the
             per-sample anomaly score chart."
    ))
  }
  invisible(list(
    df           = tibble::as_tibble(df_filtered),
    anomaly_df   = forest_dict$anomaly_df,
    outlier_list = forest_dict$outlier_list,
    anomaly_plot = anomaly_plot
  ))
}

