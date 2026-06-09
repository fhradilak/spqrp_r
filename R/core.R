# Top-level user-facing pipeline. Port of spqrp/spqrp/core.py.

#' Pairwise distances on the top-n ranked proteins
#'
#' Sub-procedure used by [perform_distance_evaluation_on_ranked_proteins()]
#' and [run_clustering()]. Selects the top-`n` proteins from `top_importance`
#' (after optionally dropping `remove_list` and after restricting the
#' ranking to proteins actually present in `df`), pivots `df` to a wide
#' matrix, and computes pairwise distances.
#'
#' The restriction to proteins present in `df` happens **before** the
#' top-n cut, so a ranking whose highest-importance entries are absent
#' from `df` still yields `n` usable proteins (the next-best ones).
#' This matches [optimize_parameters()]'s behaviour and avoids producing
#' near-empty distance matrices when the ranking and `df` were built
#' from different protein universes.
#'
#' @param top_importance Data frame with `Protein` and `Importance` columns.
#' @param n Number of top-ranked proteins to keep.
#' @param df Long-format cohort data frame.
#' @param metric See [get_distances()].
#' @param fractional_p Fractional/Minkowski exponent.
#' @param remove_list Optional character vector of proteins to drop.
#' @param number_display_neighbours Number of nearest neighbours to return.
#' @param quiet If `TRUE` (default), suppress informational status messages.
#'   Set `FALSE` to print progress and per-call summaries (sample counts,
#'   chosen cutoff, etc.). Warnings about genuine data issues -- e.g.
#'   samples dropped from the analysis -- are emitted regardless.
#' @return List with `sample_order`, `distance_matrix`, `df_dist`, and
#'   `nearest_neighbours`.
#' @keywords internal
calculate_pairwise_distances <- function(top_importance,
                                          n,
                                          df,
                                          metric = "correlation",
                                          fractional_p = 0.5,
                                          remove_list = NULL,
                                          number_display_neighbours = 1L,
                                          quiet = TRUE) {
  top_importance <- top_importance[order(-top_importance$Importance), , drop = FALSE]
  if (!is.null(remove_list)) {
    top_importance <- top_importance[!top_importance$Protein %in% remove_list, , drop = FALSE]
  }
  # Restrict the ranking to proteins actually present in df BEFORE taking
  # the top-n. Otherwise a ranking whose top-n proteins are largely
  # absent from df produces a near-empty df_dist, which silently drops
  # most samples (and then breaks UMAP downstream with
  # "n_vertices > 0 is not TRUE"). This mirrors what optimize_parameters
  # already does.
  in_df <- top_importance$Protein %in% df$Protein
  if (any(!in_df) && !quiet) {
    n_missing <- sum(!in_df)
    cli::cli_inform(c(
      "i" = "{n_missing} of {nrow(top_importance)} ranked protein{?s} {?is/are} not present in {.arg df} and will be skipped before selecting the top {n}."
    ))
  }
  top_importance <- top_importance[in_df, , drop = FALSE]
  top_importance <- utils::head(top_importance, as.integer(n))

  check_input_data_format(df, top_importance)

  df_dist <- df[df$Protein %in% top_importance$Protein, , drop = FALSE]

  res <- get_distances(df_dist, metric = metric, fractional_p = fractional_p)
  distance_matrix <- res$distance_matrix
  df_pivot <- res$df_pivot
  sample_order <- as.character(df_pivot[[SAMPLE]])

  available <- length(unique(df_dist$Sample_ID)) - 1L
  k <- min(as.integer(number_display_neighbours), available)
  nearest_neighbours <- get_nearest_neighbours(df_pivot, distance_matrix, k = k)

  invisible(list(
    sample_order = sample_order,
    distance_matrix = distance_matrix,
    df_dist = tibble::as_tibble(df_dist),
    nearest_neighbours = nearest_neighbours
  ))
}

#' Threshold-based pairwise distance evaluation
#'
#' Computes pairwise distances on the top-`n` proteins, splits sample pairs
#' by a percentile cutoff (`p`) on the distance distribution, and computes
#' classification metrics against the patient ID ground truth.
#'
#' @param df Long-format cohort data frame.
#' @param top_importance_path Optional path to a CSV with `Protein` and
#'   `Importance`. Used only when `top_importance_df` is `NULL`.
#' @param top_importance_df Optional pre-loaded ranking data frame. If
#'   supplied, `top_importance_path` is ignored. When both are `NULL`
#'   (the default) the bundled [cohort_a_ranking] dataset is used.
#' @param n Number of top-ranked proteins.
#' @param p Percentile (0-100) for the distance cutoff.
#' @param remove_list Proteins to exclude from the ranking.
#' @param metric Distance metric (see [get_distances()]).
#' @param fractional_p Fractional/Minkowski exponent.
#' @param threshold_based If `FALSE`, only return distances and skip
#'   classification.
#' @param quiet If `TRUE` (default), suppress informational status messages.
#'   Set `FALSE` to print progress and per-call summaries (sample counts,
#'   chosen cutoff, etc.). Warnings about genuine data issues -- e.g.
#'   samples dropped from the analysis -- are emitted regardless.
#' @param number_display_neighbours Number of nearest neighbours to report.
#' @param name Plot title suffix; appended to "Distribution of Pairwise
#'   Distances". Set this to a cohort label (e.g. `name = "Cohort A"`)
#'   so saved plots are self-documenting.
#' @param plot If `TRUE`, draw the distance histogram with FN/FP overlays
#'   and a legend matching the Python figure.
#' @param save_path Where to save a high-resolution render of the
#'   distance-distribution plot. Accepts `NULL` (default, don't save),
#'   `TRUE` (auto-save to a timestamped file in `tempdir()`), or a
#'   character path (e.g. `"distances.png"`). Same semantics as
#'   [run_clustering()]'s `save_path`. Only used when `plot = TRUE`.
#' @return Invisibly returns a list with `top_importance`,
#'   `nearest_neighbours`, `cutoff`, `belonging`, `not_belonging`,
#'   `eval_metrics`, `distance_matrix`, and `plot` (the ggplot built
#'   when `plot = TRUE`; `NULL` otherwise). `invisible()` keeps the
#'   REPL silent on unassigned calls. Assign to a name and use
#'   `result$plot`, `result$eval_metrics`, etc. To render the
#'   distance-distribution histogram on demand: `print(result$plot)`.
#' @examples
#' \donttest{
#' df      <- spqrp_example_data("input_cohort_df")
#' ranking <- spqrp_example_data("protein_ranking")
#' result <- perform_distance_evaluation_on_ranked_proteins(
#'   df = df, top_importance_df = ranking,
#'   metric = "manhattan", p = 0.989, n = 4L
#' )
#' result$eval_metrics[c("TP", "FP", "FN", "TN", "F1")]
#' result$plot 
#' }
#' @export
perform_distance_evaluation_on_ranked_proteins <- function(df,
                                                             top_importance_path = NULL,
                                                             top_importance_df = NULL,
                                                             n = 10L,
                                                             p = 0.5,
                                                             remove_list = NULL,
                                                             metric = "correlation",
                                                             fractional_p = 0.5,
                                                             threshold_based = TRUE,
                                                             quiet = TRUE,
                                                             number_display_neighbours = 4L,
                                                             name = "",
                                                             plot = TRUE,
                                                             save_path = NULL) {
  top_importance <- if (!is.null(top_importance_df)) {
    top_importance_df
  } else if (!is.null(top_importance_path)) {
    utils::read.csv(top_importance_path, stringsAsFactors = FALSE)
  } else {
    default_ranking()
  }

  res <- calculate_pairwise_distances(
    top_importance = top_importance,
    n = n,
    df = df,
    metric = metric,
    fractional_p = fractional_p,
    remove_list = remove_list,
    number_display_neighbours = number_display_neighbours,
    quiet = quiet
  )
  distance_matrix <- res$distance_matrix
  df_dist <- res$df_dist
  sample_order <- res$sample_order
  nearest_neighbours <- res$nearest_neighbours

  cutoff <- NULL
  belonging <- not_belonging <- NULL
  eval_metrics <- NULL
  distance_plot <- NULL

  if (threshold_based) {
    distances <- distance_matrix[upper.tri(distance_matrix)]
    cutoff <- percentile_cutoff(distances, percentile = p)

    sample_patient_mapping <- df_dist[!duplicated(df_dist$Sample_ID), c("Sample_ID", "Patient_ID")]
    spm <- stats::setNames(sample_patient_mapping$Patient_ID, sample_patient_mapping$Sample_ID)

    if (!quiet) cli::cli_text("Real number of proteins: {length(unique(df_dist$Protein))}")

    rel <- get_sample_relations_by_cutoff(
      distance_matrix, cutoff = cutoff,
      sample_patient_mapping = spm,
      sample_order = sample_order
    )
    belonging <- rel$belonging
    not_belonging <- rel$not_belonging
    eval_metrics <- get_evaluation_metrics(belonging, not_belonging, quiet = quiet)

    if (plot) {
      distance_plot <- plot_distribution_with_highlights(
        distances,
        fn_distances = eval_metrics$False_Negative_Distances,
        fp_distances = eval_metrics$False_Positive_Distances,
        percentiles = p,
        name = paste0(name, " ", n, " proteins"),
        quiet = quiet,
        save_path = save_path
      )
    }
  }

  invisible(list(
    top_importance     = tibble::as_tibble(top_importance),
    nearest_neighbours = nearest_neighbours,
    cutoff             = cutoff,
    belonging          = belonging,
    not_belonging      = not_belonging,
    eval_metrics       = eval_metrics,
    distance_matrix    = distance_matrix,
    plot               = distance_plot
  ))
}

# Default percentile grid used by optimize_parameters (matches Python).
default_percentile_grid <- function() {
  unique(c(
    25, 10,
    seq(5.0, 1.1, by = -0.1),
    seq(0.1, 0.99, by = 0.01)
  ))
}

# Default strategies for optimize_parameters
default_strategies <- function() {
  c("fp+fn", "fp", "fn", "F1", "precision", "sensitivity")
}

#' Grid-search the cutoff that optimises a chosen performance metric
#'
#' For each value of `n` (the number of top-ranked proteins) and each
#' fractional-p (only used when `metric = "fractional"`), sweeps a fixed
#' grid of percentile cutoffs and records the parameters that optimize
#' `optimization_strategy`.
#'
#' @param df Long-format cohort data frame.
#' @param metric Distance metric. `"fractional"` enables a sweep over
#'   `fractional_p_values`.
#' @param log_file Optional path; if non-NULL the optimization log is
#'   written there. Default `NULL` (no log).
#' @param top_importance_path Optional CSV path with `Protein`, `Importance`.
#'   Used only when `top_importance_df` is `NULL`.
#' @param top_importance_df Optional pre-loaded ranking. When both this and
#'   `top_importance_path` are `NULL` (the default) the bundled
#'   [cohort_a_ranking] dataset is used.
#' @param range Integer vector of `n` values to evaluate.
#' @param optimization_strategy One of `"fp+fn"`, `"fp"`, `"fn"`, `"F1"`,
#'   `"precision"`, `"sensitivity"`.
#' Optimizes for the lowest false negative (fn)
#' or false positive (fp) scores or for the highest F1, precision, sensitivity.
#' @param remove_list Proteins to drop from the ranking.
#' @param quiet If `TRUE` (default), suppress informational status messages.
#'   Set `FALSE` to print progress and per-call summaries (sample counts,
#'   chosen cutoff, etc.). Warnings about genuine data issues -- e.g.
#'   samples dropped from the analysis -- are emitted regardless.
#' @return Tibble of one row per `n`, listing the best parameters and
#'   their classification metrics.
#' @examples
#' \donttest{
#' df      <- spqrp_example_data("input_cohort_df")
#' ranking <- spqrp_example_data("protein_ranking")
#' best <- optimize_parameters(
#'   df = df, top_importance_df = ranking,
#'   metric = "manhattan", range = 2:4
#' )
#' best
#' }
#' @export
optimize_parameters <- function(df,
                                  metric = "correlation",
                                  log_file = NULL,
                                  top_importance_path = NULL,
                                  top_importance_df = NULL,
                                  range = 2:49,
                                  optimization_strategy = default_strategies(),
                                  remove_list = character(),
                                  quiet = TRUE) {
  optimization_strategy <- match.arg(optimization_strategy)
  fractional <- metric == "fractional"

  n_values <- as.integer(range)
  fractional_p_values <- if (fractional) seq(0.1, 0.99, by = 0.01) else NA_real_
  percentile_values <- default_percentile_grid()

  log_con <- if (!is.null(log_file)) {
    con <- file(log_file, open = "w")
    on.exit(close(con), add = TRUE)
    con
  } else {
    NULL
  }
  log_write <- function(text) {
    if (!is.null(log_con)) writeLines(text, log_con)
  }
  log_write("Starting parameter optimization...")

  classifier_top_proteins <- if (!is.null(top_importance_df)) {
    top_importance_df
  } else if (!is.null(top_importance_path)) {
    utils::read.csv(top_importance_path, stringsAsFactors = FALSE)
  } else {
    default_ranking()
  }

  best_results <- list()

  for (n in n_values) {
    best_result_for_n <- NULL
    best_params_for_n <- NULL
    thresholds <- list(
      `fp+fn` = Inf, fp = Inf, fn = Inf,
      f1 = 0, precision = 0, sensitivity = 0
    )

    top_importance <- classifier_top_proteins
    top_importance <- top_importance[!top_importance$Protein %in% remove_list, , drop = FALSE]
    top_prot <- top_importance[top_importance$Protein %in% df$Protein, , drop = FALSE]
    top_prot <- top_prot[order(-top_prot$Importance), , drop = FALSE]
    top_prot <- utils::head(top_prot, n)

    check_input_data_format(df, top_importance)
    df_dist <- df[df$Protein %in% top_prot$Protein, , drop = FALSE]
    spm <- df_dist[!duplicated(df_dist$Sample_ID), c("Sample_ID", "Patient_ID")]
    sample_patient_mapping <- stats::setNames(spm$Patient_ID, spm$Sample_ID)

    grid <- expand.grid(
      fractional_p = fractional_p_values,
      percentile = percentile_values,
      stringsAsFactors = FALSE
    )

    for (gi in seq_len(nrow(grid))) {
      fp_val <- grid$fractional_p[gi]
      perc <- grid$percentile[gi]

      dres <- get_distances(df_dist, metric = metric,
                              fractional_p = if (is.na(fp_val)) 0.5 else fp_val)
      distance_matrix <- dres$distance_matrix
      sample_order <- as.character(dres$df_pivot[[SAMPLE]])

      log_write(sprintf("Testing n=%d, fractional_p=%s, percentile=%s",
                          n, format(fp_val), format(perc)))

      distances <- distance_matrix[upper.tri(distance_matrix)]
      cutoff <- percentile_cutoff(distances, percentile = perc)

      rel <- get_sample_relations_by_cutoff(
        distance_matrix,
        cutoff = cutoff,
        sample_patient_mapping = sample_patient_mapping,
        sample_order = sample_order
      )
      eval_res <- get_evaluation_metrics(rel$belonging, rel$not_belonging, quiet = TRUE)
      log_write(sprintf(
        "Precision: %s, Sensitivity: %s, TP:%d, FP:%d, TN:%d, FN:%d",
        format(eval_res$Precision), format(eval_res$Sensitivity),
        eval_res$TP, eval_res$FP, eval_res$TN, eval_res$FN
      ))

      if (is_better_result(eval_res, optimization_strategy, thresholds)) {
        thresholds$`fp+fn`     <- eval_res$FP + eval_res$FN
        thresholds$fp          <- eval_res$FP
        thresholds$fn          <- eval_res$FN
        thresholds$f1          <- eval_res$F1 %||% 0
        thresholds$precision   <- eval_res$Precision
        thresholds$sensitivity <- eval_res$Sensitivity

        best_result_for_n <- eval_res
        best_params_for_n <- list(n = n, fractional_p = fp_val, percentile = perc)
      }

      if (!quiet) {
        cli::cli_text(
          "Testing n={n}, fractional_p={format(fp_val)}, percentile={format(perc)}"
        )
        cli::cli_text("Real Number Proteins: {length(unique(df_dist$Protein))}")
        cli::cli_text("Cutoff: {round(cutoff, 5)}")
      }
    }

    if (!is.null(best_result_for_n)) {
      best_results[[length(best_results) + 1L]] <- tibble::tibble(
        n            = best_params_for_n$n,
        fractional_p = best_params_for_n$fractional_p,
        percentile   = best_params_for_n$percentile,
        FP           = best_result_for_n$FP,
        FN           = best_result_for_n$FN,
        TP           = best_result_for_n$TP,
        TN           = best_result_for_n$TN,
        Precision    = best_result_for_n$Precision,
        Sensitivity  = best_result_for_n$Sensitivity,
        F1           = best_result_for_n$F1,
        Protein      = if (n - 1L < nrow(top_importance)) top_importance$Protein[n] else NA_character_
      )
    }
  }

  if (length(best_results) == 0L) {
    return(tibble::tibble())
  }
  dplyr::bind_rows(best_results)
}

# Clustering pipeline ---------------------------------------------------

#' End-to-end clustering pipeline
#'
#' Computes pairwise distances on the top-`n` ranked proteins, builds a
#' k-nearest-neighbour graph in a 2D embedding (default UMAP), iteratively
#' splits big components by max-weight edge, and visualises the result.
#'
#' @param df Long-format cohort data frame.
#' @param ranking Data frame with `Protein` and `Importance`.
#' @param n_neighbors Number of nearest-neighbour edges per sample.
#' @param max_component_size Maximum allowed connected component size.
#' @param metric Distance metric.
#' @param n Number of top-ranked proteins to use.
#' @param fractional_p Fractional/Minkowski exponent.
#' @param plot_name Plot title.
#' @param method Dimensionality reduction method (`"UMAP"`, `"PCA"`, `"MDS"`).
#' @param figsize Numeric vector of length 2: width and height in inches.
#'   Used both for `ggsave` (when `save_path` is set) and to auto-scale
#'   point sizes, line widths, and text on the plot. Larger values
#'   produce more readable plots. Default `c(16, 16)`.
#' @param dpi Resolution (dots per inch) for the saved file. Default
#'   `200` (matches Python matplotlib's default-ish output; bump to 300
#'   for print).
#' @param save_path Where to save a high-resolution PNG/SVG/PDF render.
#'   Accepts:
#'   * `NULL` (default) -- don't save; only return the ggplot object.
#'     The function still prints a hint about how to download the plot.
#'   * a character path (e.g. `"out.png"` or `"figs/cluster.svg"`) --
#'     save there via `ggsave()`. Extension chooses the format.
#' @param quiet If `TRUE` (default), suppress informational status messages.
#'   Set `FALSE` to print progress and per-call summaries (sample counts,
#'   chosen cutoff, etc.). Warnings about genuine data issues -- e.g.
#'   samples dropped from the analysis -- are emitted regardless.
#' @return Invisibly returns a list with `result_filtered`, `G` (the
#'   igraph object), `cluster_assignments`, `transitive_results`,
#'   `uncertain_samples`, `error_candidate_samples`, `plot`, and
#'   `saved_path` (the path passed in via `save_path`, or `NULL`).
#'   `invisible()` keeps the REPL silent on unassigned calls. Assign
#'   to a name to inspect; render the cluster plot on demand via
#'   `print(result$plot)`.
#' @examples
#' \donttest{
#' df      <- spqrp_example_data("input_cohort_df")
#' ranking <- spqrp_example_data("protein_ranking")
#' res <- run_clustering(
#'   df = df, ranking = ranking,
#'   n_neighbors = 1L, max_component_size = 2L,
#'   metric = "manhattan", method = "PCA"
#' )
#' head(res$cluster_assignments)
#' res$transitive_results
#' }
#' @export
run_clustering <- function(df,
                            ranking,
                            n_neighbors,
                            max_component_size,
                            metric = "manhattan",
                            n = 20L,
                            fractional_p = 0.98,
                            plot_name = "DF_Ranking_X on DF_Y",
                            method = "UMAP",
                            figsize = c(16, 16),
                            dpi = 200L,
                            save_path = NULL,
                            quiet = TRUE) {
  result_filtered <- calculate_pairwise_distances(
    top_importance = ranking,
    n = n,
    df = df,
    metric = metric,
    fractional_p = fractional_p,
    quiet = quiet
  )

  clust <- cluster_samples_iteratively(
    result = result_filtered,
    df = df,
    method = method,
    n_neighbors = n_neighbors,
    max_component_size = max_component_size,
    quiet = quiet
  )

  out <- plot_distances_neighbours_with_coloring_hue(
    df = df,
    G = clust$G,
    coords_2d = clust$coords_2d,
    method = method,
    figsize = figsize,
    dpi = dpi,
    save_path = save_path,
    df_name = plot_name,
    quiet = quiet
  )

  invisible(list(
    result_filtered          = result_filtered,
    G                        = out$G,
    cluster_assignments      = out$cluster_assignments,
    transitive_results       = out$transitive_results,
    uncertain_samples        = out$uncertain_nodes,
    error_candidate_samples  = out$error_candidates,
    saved_path               = out$saved_path,
    plot                     = out$plot
  ))
}

#' Iterative-clustering primitive: build a kNN graph + 2D coords
#'
#' Builds a 2D embedding via PCA/UMAP/MDS, connects each sample to its
#' `n_neighbors` nearest neighbours in distance space, and optionally
#' splits components larger than `max_component_size` by repeatedly
#' removing the largest-weight edge.
#'
#' @param result Output of [calculate_pairwise_distances()] (uses
#'   `distance_matrix`).
#' @param df Long-format cohort data frame.
#' @param method `"UMAP"`, `"PCA"`, or `"MDS"`.
#' @param random_state Seed for the dimensionality reduction.
#' @param n_neighbors Number of nearest-neighbour edges per sample.
#' @param max_component_size If non-NULL, iteratively split clusters
#'   above this size.
#' @param n_umap_neighbors UMAP's `n_neighbors` parameter.
#' @param precomputed_graph Optional precomputed igraph object.
#' @param mds_backend `"cmdscale"` (default) or `"smacof"` (Suggests).
#' @param quiet If `TRUE` (default), suppress informational status messages.
#'   Set `FALSE` to print progress and per-call summaries (sample counts,
#'   chosen cutoff, etc.). Warnings about genuine data issues -- e.g.
#'   samples dropped from the analysis -- are emitted regardless.
#' @return List with `G` (igraph) and `coords_2d` (matrix).
#' @keywords internal
cluster_samples_iteratively <- function(result,
                                          df,
                                          method = "UMAP",
                                          random_state = 42L,
                                          n_neighbors = 1L,
                                          max_component_size = NULL,
                                          n_umap_neighbors = 15L,
                                          precomputed_graph = NULL,
                                          mds_backend = c("cmdscale", "smacof"),
                                          quiet = TRUE) {
  mds_backend <- match.arg(mds_backend)
  n_neighbors <- as.integer(n_neighbors)
  dist_matrix <- result$distance_matrix
  # The distance matrix is the authoritative source of sample names: a
  # sample that survives `df` but has no measurements for any of the
  # top-N ranked proteins is silently dropped by the pivot inside
  # `get_distances`, so it won't appear in `dist_matrix`. Indexing
  # `dist_matrix` by such a name would error with "subscript out of
  # bounds". Pull names from the matrix itself and warn the caller if
  # any df samples got dropped.
  sample_names <- sort(rownames(dist_matrix))
  dropped <- setdiff(unique(df$Sample_ID), sample_names)
  if (length(dropped) > 0L) {
    cli::cli_warn(c(
      "{length(dropped)} sample{?s} in {.arg df} had no measurements for any of the top-ranked proteins and {?was/were} dropped from the clustering:",
      "i" = "{.val {dropped}}"
    ))
  }
  dist_matrix <- dist_matrix[sample_names, sample_names, drop = FALSE]

  graph_built <- create_graph_based_on_reduction_method(
    method = method,
    dist_matrix = dist_matrix,
    sample_names = sample_names,
    n_umap_neighbors = n_umap_neighbors,
    random_state = random_state,
    precomputed_graph = precomputed_graph,
    n_neighbors = n_neighbors,
    mds_backend = mds_backend
  )
  G <- graph_built$G
  coords_2d <- graph_built$coords_2d
  rownames(coords_2d) <- sample_names

  if (!is.null(max_component_size)) {
    G <- split_big_component_edges_by_weight(G, as.integer(max_component_size))
  }
  invisible(list(G = G, coords_2d = coords_2d))
}

#' Heavy clustering visualisation (TP hulls, FP edges, singleton markers)
#'
#' Builds the canonical SPQRP cluster plot: convex hulls around same-patient
#' true-positive clusters, dotted edges for cross-patient false-positive
#' edges, blue square markers for true-positive singletons, pink circles
#' for uncertain (isolated but should-be-connected) samples. Returns the
#' ggplot plus cluster bookkeeping.
#'
#' @param df Long-format cohort data frame.
#' @param G igraph object from [cluster_samples_iteratively()].
#' @param coords_2d 2D coordinates from [cluster_samples_iteratively()].
#' @param method Reduction method (`"UMAP"`, `"PCA"`, `"MDS"`).
#' @param subset_samples Optional sample subset to visualise.
#' @param highlight_singletons Mark same-patient-singletons with blue squares.
#' @param highlight_single_samples_missing_connections Mark uncertain
#'   samples with pink circles.
#' @param figsize Numeric vector of length 2: width and height in inches.
#'   Drives both the `ggsave` output dimensions (when `save_path` is set)
#'   and the auto-scaling of point sizes, line widths, fonts, and theme
#'   `base_size`. Default `c(14, 14)`. Use `c(20, 20)` or larger for
#'   publication-quality renders.
#' @param dpi Resolution (dots per inch) for the saved file. Default
#'   `150` matches the Python package's matplotlib default. Use `300`
#'   for print-quality.
#' @param label_patient_only Label nodes by patient instead of sample.
#' @param label_offset_x,label_offset_y Label nudge offsets. Values
#'   <= 0.05 are interpreted as a fraction of the coord range (auto-
#'   scaled to the data); larger values are absolute.
#' @param label_font Override the auto-scaled label font size. `NULL`
#'   (default) lets the function pick a size from `figsize`.
#' @param df_name Plot title.
#' @param save_path If non-NULL, save plot to this file via `ggsave`.
#' @param print If `TRUE`, print the plot.
#' @param quiet If `TRUE` (default), suppress informational status messages
#'   (save-path hints, cluster summaries, transitive performance metrics).
#'   Set `FALSE` to print them. Warnings about genuine data issues -- e.g.
#'   samples dropped from the analysis -- are emitted regardless.
#' @return List with `plot`, `G`, `cluster_assignments`,
#'   `transitive_results`, `uncertain_nodes`, `error_candidates`.
#' @keywords internal
plot_distances_neighbours_with_coloring_hue <- function(df,
                                                          G,
                                                          coords_2d,
                                                          method = "UMAP",
                                                          subset_samples = NULL,
                                                          highlight_singletons = TRUE,
                                                          highlight_single_samples_missing_connections = TRUE,
                                                          figsize = c(14, 14),
                                                          dpi = 150L,
                                                          label_patient_only = FALSE,
                                                          label_offset_x = 0.01,
                                                          label_offset_y = 0.01,
                                                          label_font = NULL,
                                                          df_name = "DF_NAME",
                                                          save_path = NULL,
                                                          print = TRUE,
                                                          quiet = TRUE) {
  # Auto-scale visual elements with figsize. Reference size = 14 inches;
  # all multipliers are clamped to (0.5, 2.5) so tiny or huge figsizes
  # don't produce absurd sizes.
  scale_factor <- max(0.5, min(2.5, mean(figsize) / 14))
  if (is.null(label_font)) {
    label_font <- 3.2 * scale_factor
  }
  base_size <- 12 * scale_factor
  node_size <- 2.4 * scale_factor
  fp_linewidth <- 1.6 * scale_factor
  hull_edge_linewidth <- 1.2 * scale_factor
  hull_outline_linewidth <- 0.4 * scale_factor
  highlight_size_square <- 4 * scale_factor
  highlight_size_circle <- 5 * scale_factor
  fp_color <- "#DC267F"
  tp_color <- "#1AFF1A"
  error_candidate_color <- "#FFBBF6"

  sample_to_patient <- stats::setNames(df$Patient_ID, df$Sample_ID)
  sample_to_patient <- sample_to_patient[!duplicated(names(sample_to_patient))]
  samples_by_patient <- split(names(sample_to_patient), sample_to_patient)
  # Intersect df's samples with the graph's vertices: a df sample without
  # measurements for any of the top ranked proteins was dropped during
  # distance computation and is therefore absent from G / coords_2d.
  # Indexing coords_2d by such a name would error with "subscript out of
  # bounds" further down in node_df construction.
  sample_names <- sort(intersect(unique(df$Sample_ID), igraph::V(G)$name))
  sample_index <- stats::setNames(seq_along(sample_names), sample_names)

  drawn_pairs_mat <- igraph::as_edgelist(G)
  drawn_pairs <- if (nrow(drawn_pairs_mat) > 0L) {
    lapply(seq_len(nrow(drawn_pairs_mat)), function(i) sort(drawn_pairs_mat[i, ]))
  } else {
    list()
  }
  connected_samples <- if (is.null(subset_samples)) {
    sample_names
  } else {
    intersect(sample_names, subset_samples)
  }

  classifications <- identify_clusters_singletons(
    G = G,
    sample_to_patient = sample_to_patient,
    samples_by_patient = samples_by_patient,
    drawn_pairs = drawn_pairs,
    sample_names = sample_names
  )

  # Build the layered ggplot
  node_df <- data.frame(
    sample = sample_names,
    x = coords_2d[sample_names, 1L],
    y = coords_2d[sample_names, 2L],
    patient = unname(sample_to_patient[sample_names]),
    stringsAsFactors = FALSE
  )
  node_df$color <- ifelse(
    node_df$sample %in% classifications$singleton_nodes |
      node_df$sample %in% classifications$nodes_in_tp_clusters,
    tp_color,
    ifelse(node_df$sample %in% classifications$nodes_in_fp_cluster, fp_color,
            ifelse(node_df$sample %in% classifications$isolated_nodes,
                    error_candidate_color, tp_color))
  )
  node_df <- node_df[node_df$sample %in% connected_samples, , drop = FALSE]

  hulls_edges <- tp_cluster_hull_data(
    samples_by_patient = samples_by_patient,
    connected_samples = connected_samples,
    drawn_pairs = drawn_pairs,
    sample_index = sample_index,
    coords_2d = coords_2d
  )

  fp_edges <- list()
  for (pair in drawn_pairs) {
    if (!any(pair %in% connected_samples)) next
    if (sample_to_patient[[pair[1L]]] != sample_to_patient[[pair[2L]]]) {
      fp_edges[[length(fp_edges) + 1L]] <- data.frame(
        x = coords_2d[pair[1L], 1L], y = coords_2d[pair[1L], 2L],
        xend = coords_2d[pair[2L], 1L], yend = coords_2d[pair[2L], 2L]
      )
    }
  }
  fp_edges_df <- if (length(fp_edges) > 0L) do.call(rbind, fp_edges) else NULL

  p <- ggplot2::ggplot()

  # Convex hulls (TP clusters of 3+ same-patient samples)
  for (hull in hulls_edges$hulls) {
    p <- p + ggplot2::geom_polygon(
      data = hull,
      mapping = ggplot2::aes(x = .data$x, y = .data$y),
      fill = tp_color, alpha = 0.45, colour = "#006400",
      linewidth = hull_outline_linewidth
    )
  }
  for (edge in hulls_edges$edges) {
    p <- p + ggplot2::geom_segment(
      data = edge,
      mapping = ggplot2::aes(x = .data$x, y = .data$y,
                              xend = .data$xend, yend = .data$yend),
      colour = tp_color, linewidth = hull_edge_linewidth
    )
  }

  if (!is.null(fp_edges_df) && nrow(fp_edges_df) > 0L) {
    p <- p + ggplot2::geom_segment(
      data = fp_edges_df,
      mapping = ggplot2::aes(x = .data$x, y = .data$y,
                              xend = .data$xend, yend = .data$yend),
      colour = fp_color, linetype = "dotted", linewidth = fp_linewidth
    )
  }

  p <- p +
    ggplot2::geom_point(
      data = node_df,
      mapping = ggplot2::aes(x = .data$x, y = .data$y),
      colour = node_df$color, size = node_size
    )

  # Scale label offset by the actual coord range so user's defaults
  # (0.01) don't disappear on PCA/MDS axes with span = 100+ units.
  x_range <- diff(range(coords_2d[, 1L]))
  y_range <- diff(range(coords_2d[, 2L]))
  ox <- if (label_offset_x <= 0.05) label_offset_x * x_range else label_offset_x
  oy <- if (label_offset_y <= 0.05) label_offset_y * y_range else label_offset_y

  label_df <- node_df
  label_df$label <- if (label_patient_only) label_df$patient else label_df$sample
  p <- p + ggplot2::geom_text(
    data = label_df,
    mapping = ggplot2::aes(x = .data$x + ox,
                            y = .data$y + oy,
                            label = .data$label),
    size = label_font, hjust = 0
  )

  if (highlight_singletons && length(classifications$singleton_nodes) > 0L) {
    sing_df <- node_df[node_df$sample %in% classifications$singleton_nodes, , drop = FALSE]
    if (nrow(sing_df) > 0L) {
      p <- p + ggplot2::geom_point(
        data = sing_df,
        mapping = ggplot2::aes(x = .data$x, y = .data$y),
        shape = 0, size = highlight_size_square, colour = "blue"
      )
    }
  }
  if (highlight_single_samples_missing_connections &&
      length(classifications$isolated_nodes) > 0L) {
    iso_df <- node_df[node_df$sample %in% classifications$isolated_nodes, , drop = FALSE]
    if (nrow(iso_df) > 0L) {
      p <- p + ggplot2::geom_point(
        data = iso_df,
        mapping = ggplot2::aes(x = .data$x, y = .data$y),
        shape = 1, size = highlight_size_circle,
        colour = error_candidate_color, stroke = 1.2 * scale_factor
      )
    }
  }

  # ----- Legend construction --------------------------------------------
  # All actual rendering above uses explicit (non-mapped) colours so the
  # ggplot legend machinery skips them. Add invisible dummy layers below
  # that ARE aesthetic-mapped so ggplot generates a legend matching the
  # Python original:
  #   * True positive cluster   -> filled green square (shape 22)
  #   * True positive singleton -> empty blue square   (shape 22, fill NA)
  #   * Uncertain sample        -> empty pink circle   (shape 21, fill NA)
  #   * Error candidate         -> magenta dotted line ("3 dots" key)
  #
  # Approach: drive everything off ONE colour scale per geom, then use
  # `override.aes` in guides() to set per-entry shape, fill, and stroke.
  # That avoids ggplot's quirky behaviour around merging fill+shape legends.
  node_levels <- c("True positive cluster",
                    "True positive singleton",
                    "Uncertain sample")
  node_legend_df <- data.frame(
    x = rep(NA_real_, length(node_levels)),
    y = rep(NA_real_, length(node_levels)),
    category = factor(node_levels, levels = node_levels),
    stringsAsFactors = FALSE
  )
  edge_legend_df <- data.frame(
    x = NA_real_, y = NA_real_, xend = NA_real_, yend = NA_real_,
    edge_category = factor("Error candidate", levels = "Error candidate"),
    stringsAsFactors = FALSE
  )

  node_stroke_colours <- c(
    "True positive cluster"   = "#006400",
    "True positive singleton" = "blue",
    "Uncertain sample"        = error_candidate_color
  )

  p <- p +
    ggplot2::geom_point(
      data = node_legend_df,
      mapping = ggplot2::aes(x = .data$x, y = .data$y,
                              colour = .data$category),
      size = node_size, na.rm = TRUE
    ) +
    ggplot2::geom_segment(
      data = edge_legend_df,
      mapping = ggplot2::aes(x = .data$x, y = .data$y,
                              xend = .data$xend, yend = .data$yend,
                              linetype = .data$edge_category),
      colour = fp_color, linewidth = fp_linewidth, na.rm = TRUE
    ) +
    ggplot2::scale_colour_manual(
      name = NULL,
      values = node_stroke_colours,
      breaks = node_levels
    ) +
    ggplot2::scale_linetype_manual(
      name = NULL,
      values = c("Error candidate" = "dotted")
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        order = 1L,
        override.aes = list(
          shape  = c(22L, 22L, 21L),
          fill   = c(tp_color, NA, NA),
          colour = c("#006400", "blue", error_candidate_color),
          stroke = c(1.0, 1.4, 1.4) * scale_factor,
          size   = c(6, 6, 6) * scale_factor
        )
      ),
      linetype = ggplot2::guide_legend(
        order = 2L,
        override.aes = list(colour    = fp_color,
                             linewidth = fp_linewidth)
      )
    )

  axis_labels <- switch(
    toupper(method),
    "PCA"  = list("PCA 1", "PCA 2"),
    "MDS"  = list("MDS Dimension 1", "MDS Dimension 2"),
    list(paste0(toupper(method), " Dim 1"), paste0(toupper(method), " Dim 2"))
  )
  p <- p +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.08, 0.18))
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.08, 0.12))
    ) +
    ggplot2::labs(
      title = paste0(toupper(method), " projection of ",
                      df_name, " sample Clustering"),
      x = axis_labels[[1L]], y = axis_labels[[2L]]
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = base_size * 1.2),
      plot.margin = ggplot2::margin(12, 30, 12, 12),
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.spacing.y = ggplot2::unit(8, "pt"),
      legend.text = ggplot2::element_text(size = base_size * 0.95),
      legend.key.size = ggplot2::unit(1.6, "lines"),
      legend.key.height = ggplot2::unit(1.4, "lines")
    )

  if (is.null(save_path)) {

    if (!quiet) {
      cli::cli_inform(c(
        "i" = paste0(
          "Tip: pass {.code save_path = \"cluster.png\"} to save a high-",
          "resolution version of this plot. RStudio's Plot pane re-",
          "renders on resize, which can distort proportions; a saved ",
          "PNG/SVG does not."
        )
      ))
    }

  } else {

    ggplot2::ggsave(
      filename = save_path,
      plot = p,
      width = figsize[1L],
      height = figsize[2L],
      dpi = dpi
    )

    if (!quiet) {
      px_w <- round(figsize[1L] * dpi)
      px_h <- round(figsize[2L] * dpi)

      cli::cli_inform(c(
        "v" = "Saved plot to {.path {save_path}}",
        "i" = paste0(
          "Dimensions: {figsize[1L]}x{figsize[2L]} in @ {dpi} dpi ",
          "({px_w}x{px_h} px). Open in an image viewer for a fixed, ",
          "undistorted render."
        )
      ))
    }
  }

  if (print && !quiet) {
    # Open a sized PDF in tempdir() for non-interactive runs (default
    # would dump 7x7" Rplots.pdf in cwd -- CRAN-noncompliant).
    if (!interactive() && grDevices::dev.cur() == 1L) {
      tmp_pdf <- tempfile(pattern = "spqrp_cluster_", fileext = ".pdf")
      grDevices::pdf(file = tmp_pdf,
                     width = figsize[1L], height = figsize[2L])
      on.exit(grDevices::dev.off(), add = TRUE)
      cli::cli_inform(c(
        "i" = "Auto-rendered cluster plot to {.path {tmp_pdf}} (no display device available)."
      ))
    }
    print(p)
    grDevices::dev.flush()
  }

  comps <- igraph::components(G)
  cluster_assignments <- stats::setNames(
    as.integer(comps$membership),
    igraph::V(G)$name
  )

  if (!print && !quiet) {
    # Emit cluster summary to stdout to mirror Python's print() behaviour.
    vnames <- igraph::V(G)$name
    for (cid in seq_len(comps$no)) {
      mem <- sort(vnames[comps$membership == cid])
      cli::cli_text(
        "Cluster {cid} ({length(mem)} samples): {.val {mem}}"
      )
    }
  }

  transitive_results <- transitive_performance(
    sample_names = sample_names,
    drawn_pairs = drawn_pairs,
    sample_to_patient = sample_to_patient,
    quiet = quiet
  )

  invisible(list(
    plot                = p,
    G                   = G,
    cluster_assignments = cluster_assignments,
    transitive_results  = transitive_results,
    uncertain_nodes     = classifications$isolated_nodes,
    error_candidates    = classifications$nodes_in_fp_cluster,
    saved_path          = save_path
  ))
}
