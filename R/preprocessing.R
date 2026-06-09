# Preprocessing: protein-occurrence filter, log transform,
# per-sample median normalization, plate-effect residualisation.
# Port of spqrp/spqrp/preprocessing.py.

#' Keep proteins present in at least a given fraction of samples
#'
#' @param df Long-format intensity data frame.
#' @param cutoff Fraction in `[0, 1]`. A protein is kept when its non-NA
#'   intensity covers at least `cutoff` of the samples in `df`.
#' @return Filtered tibble in the same shape as `df`.
#' @examples
#' df <- spqrp_example_data("input_cohort_df")
#' kept <- filter_by_occurrence(df, cutoff = 0.7)
#' length(unique(kept$Protein))
#' @export
filter_by_occurrence <- function(df, cutoff = 0.7) {
  if (!is.numeric(cutoff) || length(cutoff) != 1L || cutoff < 0 || cutoff > 1) {
    cli::cli_abort("{.arg cutoff} must be a single number in [0, 1].")
  }
  n_samples <- length(unique(df[[SAMPLE]]))
  counts <- tapply(df$Intensity, df$Protein, function(x) sum(!is.na(x)))
  fractions <- counts / n_samples
  proteins_to_keep <- names(fractions)[fractions >= cutoff]
  tibble::as_tibble(df[df$Protein %in% proteins_to_keep, , drop = FALSE])
}

#' Log2-transform the intensity column in place (long format)
#'
#' @param df Long-format intensity data frame.
#' @return `df` with `Intensity = log2(Intensity)`.
#' @examples
#' df <- spqrp_example_data("input_cohort_df")
#' head(log_transform(df))
#' @export
log_transform <- function(df) {
  df$Intensity <- log2(df$Intensity)
  df
}

#' Inverse of log2-transform: raise intensities to the power of 2
#'
#' @inheritParams log_transform
#' @return `df` with `Intensity = 2^Intensity`.
#' @keywords internal
revert_log_transform <- function(df) {
  df$Intensity <- 2 ^ df$Intensity
  df
}

#' Per-sample median normalisation (log-space subtraction)
#'
#' Subtracts each sample's median intensity from its intensities, then
#' re-centers on the dataset's overall median. By default `df` is assumed
#' to already be in log space. If `revert_log = TRUE`, the function reverts
#' the log transform first, then divides by the per-sample median ratio.
#'
#' Returns a list with `data` (the normalized tibble) and `plot` (a
#' ggplot showing before/after boxplots).
#'
#' @param dataset Long-format intensity data frame.
#' @param string_of_pool If non-empty, samples whose ID contains this
#'   substring are excluded from normalization (kept out of the
#'   post-normalization data).
#' @param revert_log If `TRUE`, run [revert_log_transform()] first.
#' @param sample Column to group by (defaults to `"Sample_ID"`).
#' @param plot If `TRUE`, attach a before/after boxplot.
#' @return List with `data` and (optionally) `plot`.
#' @examples
#' df <- spqrp_example_data("input_cohort_df")
#' norm <- normalize_medianintensity(log_transform(df), plot = FALSE)
#' head(norm$data)
#' @export
normalize_medianintensity <- function(dataset,
                                       string_of_pool = "",
                                       revert_log = FALSE,
                                       sample = SAMPLE,
                                       plot = TRUE) {
  if (revert_log) {
    dataset <- revert_log_transform(dataset)
  }

  per_sample_median <- tapply(dataset$Intensity, dataset[[sample]], stats::median, na.rm = TRUE)
  median_of_summed_intensities <- stats::median(per_sample_median, na.rm = TRUE)
  dataset_before <- dataset

  if (nzchar(string_of_pool)) {
    dataset <- dataset[!grepl(string_of_pool, dataset[[sample]], fixed = TRUE), , drop = FALSE]
  }

  # Per-sample median scaling
  global_median <- stats::median(dataset$Intensity, na.rm = TRUE)
  if (revert_log) {
    dataset <- dplyr::group_by(dataset, .data[[sample]])
    dataset <- dplyr::mutate(
      dataset,
      NormalizeFactor = median_of_summed_intensities / stats::median(.data$Intensity, na.rm = TRUE),
      Intensity = .data$Intensity * (median_of_summed_intensities / stats::median(.data$Intensity, na.rm = TRUE))
    )
    dataset <- dplyr::ungroup(dataset)
  } else {
    dataset <- dplyr::group_by(dataset, .data[[sample]])
    dataset <- dplyr::mutate(
      dataset,
      NormalizeFactor = stats::median(.data$Intensity, na.rm = TRUE) - global_median,
      Intensity = .data$Intensity - (stats::median(.data$Intensity, na.rm = TRUE) - global_median)
    )
    dataset <- dplyr::ungroup(dataset)
  }

  p <- NULL
  if (plot) {
    plot_df <- dplyr::bind_rows(
      dplyr::mutate(dataset_before, stage = "Before normalization"),
      dplyr::mutate(dataset, stage = "After normalization")
    )
    plot_df$stage <- factor(plot_df$stage,
                            levels = c("Before normalization", "After normalization"))
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = .data[[sample]], y = .data$Intensity)
    ) +
      ggplot2::geom_boxplot(outlier.size = 0.5) +
      ggplot2::facet_wrap(~ stage, ncol = 2L) +
      ggplot2::labs(x = sample, y = "Intensity (log2)") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, size = 6))
  }

  list(data = tibble::as_tibble(dataset), plot = p)
}

#' Regress intensity on plate and replace it with OLS residuals
#'
#' Identifies a `plate` (or `Plate`) column, encodes it as integers, fits
#' `lm(Intensity ~ plate)`, and replaces `Intensity` with the model's
#' residuals. If no plate column is present, returns the input unchanged
#' with a message.
#'
#' @param group_data Long-format intensity data frame.
#' @param individual Patient identifier column (default `"Patient_ID"`).
#' @param sample Sample identifier column (default `"Sample_ID"`).
#' @param impute If `TRUE`, impute missing intensities by patient/protein
#'   median before regression; otherwise drop NA rows.
#' @param verbose If `TRUE`, also build before/after boxplots.
#' @return Tibble with corrected `Intensity`. Attribute `"plot"` carries
#'   the diagnostic ggplot when `verbose = TRUE`.
#' @examples
#' \donttest{
#' df <- spqrp_example_data("input_cohort_df")
#' df$plate <- rep(c("A", "B"), length.out = nrow(df))
#' corrected <- plate_correct_residuals_by_protein(df)
#' head(corrected)
#' }
#' @export
plate_correct_residuals_by_protein <- function(group_data,
                                                 individual = PATIENT,
                                                 sample = SAMPLE,
                                                 impute = FALSE,
                                                 verbose = FALSE) {
  if (impute) {
    group_data <- dplyr::group_by(group_data, .data$Protein, .data[[individual]])
    group_data <- dplyr::mutate(
      group_data,
      Intensity = ifelse(is.na(.data$Intensity),
                         stats::median(.data$Intensity, na.rm = TRUE),
                         .data$Intensity)
    )
    group_data <- dplyr::ungroup(group_data)
  } else {
    group_data <- group_data[!is.na(group_data$Intensity), , drop = FALSE]
  }

  plate_col <- intersect(c("plate", "Plate"), names(group_data))[1L]
  if (is.na(plate_col)) {
    cli::cli_warn(
      "No 'plate' or 'Plate' column found. Skipping plate effect correction."
    )
    return(tibble::as_tibble(group_data))
  }

  group_data[[plate_col]] <- as.integer(factor(as.character(group_data[[plate_col]])))
  original_intensity <- group_data$Intensity

  fit <- tryCatch(
    stats::lm(stats::as.formula(paste("Intensity ~", plate_col)), data = group_data),
    error = function(e) {
      cli::cli_warn("OLS model fitting failed: {e$message}")
      NULL
    }
  )
  if (is.null(fit)) return(tibble::as_tibble(group_data))

  coefs <- summary(fit)$coefficients
  plate_p_value <- if (plate_col %in% rownames(coefs)) {
    coefs[plate_col, "Pr(>|t|)"]
  } else {
    NA_real_
  }
  if (!is.na(plate_p_value) && plate_p_value >= 0.05) {
    cli::cli_warn(
      c("Plate effect is not statistically significant (p = {round(plate_p_value, 4)}).",
        "i" = "Consider whether plate correction is appropriate.")
    )
  }

  group_data$Intensity <- as.numeric(stats::residuals(fit))

  if (verbose) {
    diag_df <- dplyr::bind_rows(
      dplyr::mutate(group_data, Intensity = original_intensity, stage = "Before"),
      dplyr::mutate(group_data, stage = "After")
    )
    diag_df$stage <- factor(diag_df$stage, levels = c("Before", "After"))
    p_protein <- ggplot2::ggplot(diag_df,
                                  ggplot2::aes(x = .data$Protein, y = .data$Intensity)) +
      ggplot2::geom_boxplot(outlier.size = 0.4) +
      ggplot2::facet_wrap(~ stage, ncol = 2L) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, size = 7)) +
      ggplot2::labs(title = "Plate correction (per protein)")
    attr(group_data, "plot") <- p_protein
  }

  tibble::as_tibble(group_data)
}
