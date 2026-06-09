# Static plot helpers. Ports of:
# - plot_distribution_of_pairwise_dist
# - plot_distribution_with_highlights
# The heavy clustering plot (plot_distances_neighbours_with_coloring_hue)
# lives in core.R since it forms the user-facing output.

#' Histogram of pairwise distances with optional percentile lines
#'
#' @param distances Numeric vector.
#' @param percentiles Vector of percentile values (0-100) to mark with
#'   vertical lines.
#' @param print If `TRUE`, auto-render the plot (gated by `quiet`).
#' @param quiet If `TRUE` (default), suppress informational status messages
#'   and skip auto-rendering of the returned ggplot. Set `FALSE` to
#'   render. Warnings about genuine data issues are emitted regardless.
#' @return Invisible ggplot.
#' @keywords internal
plot_distribution_of_pairwise_dist <- function(distances,
                                                 percentiles = c(1, 2, 5, 10),
                                                 print = TRUE,
                                                 quiet = TRUE) {
  p <- ggplot2::ggplot(
    data = data.frame(distance = distances),
    mapping = ggplot2::aes(x = .data$distance)
  ) +
    ggplot2::geom_histogram(bins = 50L, fill = "grey60", colour = "black", alpha = 0.7) +
    ggplot2::labs(x = "Pairwise Distance", y = "Frequency",
                  title = "Distribution of Pairwise Distances")

  for (q in percentiles) {
    p <- p + ggplot2::geom_vline(
      xintercept = stats::quantile(distances, q / 100, type = 7L, na.rm = TRUE),
      colour = "#DC267F", linetype = "dashed", linewidth = 0.8
    )
  }
  if (print && !quiet) print(p)
  invisible(p)
}

#' Histogram with FN/FP/percentile overlays and a legend
#'
#' Mirrors Python's `plot_distribution_with_highlights` (helpers.py): a grey
#' histogram of all pairwise distances, overlaid with vertical lines for
#' false-negative pairs (blue), false-positive pairs (orange), and percentile
#' cutoffs (magenta). The legend names each category, matching the Python
#' figure key.
#'
#' @param distances Numeric vector of all pairwise distances.
#' @param fn_distances Distances of false-negative pairs. Pass `numeric(0)` or
#'   `NULL` to omit the FN legend entry.
#' @param fp_distances Distances of false-positive pairs. Pass `numeric(0)` or
#'   `NULL` to omit the FP legend entry.
#' @param percentiles Percentiles (0-100) to draw as vertical lines. Each
#'   percentile becomes its own legend entry.
#' @param name Title suffix appended to "Distribution of Pairwise Distances".
#'   Set this to a cohort label so saved plots are self-documenting.
#' @param print If `TRUE`, print the plot.
#' @param figsize Numeric vector of length 2: width and height in inches for
#'   `ggsave` output. Default `c(8, 5)` matches Python's `plt.figure(figsize=(8, 5))`.
#' @param dpi Resolution (dots per inch) for the saved file. Default `150`.
#' @param save_path Where to save a high-resolution PNG/SVG/PDF render.
#'   Accepts:
#'   * `NULL` (default) -- don't save; only return the ggplot object. The
#'     function prints a hint about how to save.
#'   * a character path (e.g. `"distances.png"`) -- save there via
#'     `ggsave()`. Extension chooses the format.
#' @param quiet If `TRUE` (default), suppress the informational
#'   `save_path` hints. Warnings about genuine data issues are emitted
#'   regardless.
#' @return Invisible ggplot.
#' @keywords internal
plot_distribution_with_highlights <- function(distances,
                                                fn_distances,
                                                fp_distances,
                                                percentiles = c(1, 2, 5, 10),
                                                name = "",
                                                print = TRUE,
                                                figsize = c(8, 5),
                                                dpi = 150L,
                                                save_path = NULL,
                                                quiet = TRUE) {
  hist_color <- "grey60"
  fn_color   <- "#648FFF"
  fp_color   <- "#FFB000"
  pct_color  <- "#DC267F"

  # Build the colour-scale level set dynamically -- only include categories
  # we actually have data for (matches Python's first-iteration-only label
  # injection in helpers.py).
  has_fn  <- length(fn_distances) > 0L
  has_fp  <- length(fp_distances) > 0L
  has_pct <- length(percentiles) > 0L
  percentile_labels <- if (has_pct) sprintf("%sth Percentile", percentiles) else character()
  levels_all <- c("All Distances",
                   if (has_fn)  "False Negatives" else NULL,
                   if (has_fp)  "False Positives" else NULL,
                   percentile_labels)

  colour_values <- c(
    "All Distances"   = hist_color,
    "False Negatives" = fn_color,
    "False Positives" = fp_color,
    stats::setNames(rep(pct_color, length(percentile_labels)), percentile_labels)
  )

  p <- ggplot2::ggplot(
    data = data.frame(distance = distances),
    mapping = ggplot2::aes(x = .data$distance)
  ) +
    ggplot2::geom_histogram(bins = 40L, fill = hist_color, colour = NA, alpha = 0.5) +
    ggplot2::labs(
      x = "Distance", y = "Frequency",
      title = paste0("Distribution of Pairwise Distances ", name)
    )

  # FN vertical lines -- aesthetic-mapped so the legend picks them up.
  if (has_fn) {
    fn_df <- data.frame(
      xintercept = as.numeric(fn_distances),
      category   = factor("False Negatives", levels = levels_all)
    )
    p <- p + ggplot2::geom_vline(
      data = fn_df,
      mapping = ggplot2::aes(xintercept = .data$xintercept, colour = .data$category),
      linetype = "solid", linewidth = 0.4, alpha = 0.2, show.legend = TRUE
    )
  }
  # FP vertical lines.
  if (has_fp) {
    fp_df <- data.frame(
      xintercept = as.numeric(fp_distances),
      category   = factor("False Positives", levels = levels_all)
    )
    p <- p + ggplot2::geom_vline(
      data = fp_df,
      mapping = ggplot2::aes(xintercept = .data$xintercept, colour = .data$category),
      linetype = "solid", linewidth = 0.4, alpha = 0.2, show.legend = TRUE
    )
  }
  # Percentile cutoff lines (one factor level each, all share the magenta colour).
  if (has_pct) {
    pct_df <- data.frame(
      xintercept = as.numeric(stats::quantile(distances, percentiles / 100,
                                                type = 7L, na.rm = TRUE)),
      category   = factor(percentile_labels, levels = levels_all)
    )
    p <- p + ggplot2::geom_vline(
      data = pct_df,
      mapping = ggplot2::aes(xintercept = .data$xintercept, colour = .data$category),
      linewidth = 1.0, show.legend = TRUE
    )
  }

  # Dummy point layer to inject "All Distances" as the first legend entry --
  # the histogram itself doesn't map to an aesthetic, so without this it
  # would be missing from the key.
  dummy_hist_df <- data.frame(
    x = NA_real_, y = NA_real_,
    category = factor("All Distances", levels = levels_all)
  )
  p <- p +
    ggplot2::geom_point(
      data = dummy_hist_df,
      mapping = ggplot2::aes(x = .data$x, y = .data$y, colour = .data$category),
      size = 5, na.rm = TRUE
    ) +
    ggplot2::scale_colour_manual(
      name = NULL,
      values = colour_values,
      breaks = levels_all,
      drop = FALSE
    )

  # Override the legend keys: filled square for "All Distances", thin solid
  # lines for FN/FP, thick solid lines for percentile cutoffs.
  n_pct <- length(percentile_labels)
  override_shape     <- c(15L,        # All Distances -> filled square (alpha 0.5 already on geom)
                           if (has_fn)  NA else NULL,
                           if (has_fp)  NA else NULL,
                           rep(NA, n_pct))
  override_linetype  <- c("blank",
                           if (has_fn)  "solid" else NULL,
                           if (has_fp)  "solid" else NULL,
                           rep("solid", n_pct))
  override_linewidth <- c(0,
                           if (has_fn)  0.9 else NULL,
                           if (has_fp)  0.9 else NULL,
                           rep(1.6, n_pct))
  override_alpha     <- c(0.5,
                           if (has_fn)  1 else NULL,
                           if (has_fp)  1 else NULL,
                           rep(1, n_pct))
  p <- p + ggplot2::guides(
    colour = ggplot2::guide_legend(
      override.aes = list(
        shape     = override_shape,
        linetype  = override_linetype,
        linewidth = override_linewidth,
        alpha     = override_alpha
      )
    )
  )

  # Resolve save_path:
  #   NULL      -> don't save (but optionally hint)
  #   character -> save there
  if (is.null(save_path)) {
    if (!quiet) {
      cli::cli_inform(c(
        "i" = paste0(
          "Tip: pass {.code save_path = \"distances.png\"} to save a ",
          "high-resolution version of this plot."
        )
      ))
    }
  } else {
    ggplot2::ggsave(save_path, plot = p,
                    width = figsize[1L], height = figsize[2L],
                    dpi = dpi)
    if (!quiet) {
      px_w <- round(figsize[1L] * dpi)
      px_h <- round(figsize[2L] * dpi)
      cli::cli_inform(c(
        "v" = "Saved distance-distribution plot to {.path {save_path}}",
        "i" = paste0(
          "Dimensions: {figsize[1L]}x{figsize[2L]} in @ {dpi} dpi ",
          "({px_w}x{px_h} px)."
        )
      ))
    }
  }

  # `quiet` is the master verbosity switch: when TRUE (default) the plot
  # is built and returned but not auto-rendered. Users who want to see
  # it without changing quiet can call `print(result_of_perform_…$plot)`.
  if (print && !quiet) print(p)
  invisible(p)
}
