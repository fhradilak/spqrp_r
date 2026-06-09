# Constants & internal utilities

SAMPLE  <- "Sample_ID"
PATIENT <- "Patient_ID"

REQUIRED_COLUMNS_DF      <- c("Sample_ID", "Patient_ID", "Intensity", "Protein")
REQUIRED_COLUMNS_RANKING <- c("Protein", "Importance")

# Assertions ------------------------------------------------------------

assert_required_columns <- function(df, required, type,
                                    arg = rlang::caller_arg(df),
                                    call = rlang::caller_env()) {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "{.arg {arg}} ({type}) is missing required column{?s}: {.val {missing}}.",
        "i" = "Required columns for {type}: {.val {required}}.",
        "i" = "Rename the columns of your dataset to match the requirements."
      ),
      call = call
    )
  }
  invisible(TRUE)
}

#' Validate required columns of a cohort and (optionally) a ranking
#'
#' Throws an informative error when required columns are missing.
#'
#' @param df Cohort data frame; must contain `Sample_ID`, `Patient_ID`,
#'   `Protein`, `Intensity`.
#' @param importance_ranking Optional ranking data frame; must contain
#'   `Protein`, `Importance` if supplied.
#' @return Invisible `TRUE`.
#' @examples
#' df <- spqrp_example_data("input_cohort_df")
#' ranking <- spqrp_example_data("protein_ranking")
#' check_input_data_format(df, ranking)
#' @export
check_input_data_format <- function(df, importance_ranking = NULL) {
  assert_required_columns(df, REQUIRED_COLUMNS_DF, "protein matrix")
  if (!is.null(importance_ranking)) {
    assert_required_columns(
      importance_ranking, REQUIRED_COLUMNS_RANKING, "importance ranking"
    )
  }
  invisible(TRUE)
}

# Path helpers ----------------------------------------------------------

#' Filesystem path to a bundled example CSV
#'
#' @param which One of `"input_cohort_df"`, `"protein_ranking"`.
#' @return Absolute character path inside `inst/extdata/`.
#' @examples
#' spqrp_example_path("input_cohort_df")
#' @export
spqrp_example_path <- function(which = c("input_cohort_df",
                                         "protein_ranking")) {
  which <- match.arg(which)
  file <- switch(
    which,
    input_cohort_df  = "example_input_cohort_df.csv",
    protein_ranking  = "example_protein_ranking.csv"
  )
  system.file("extdata", file, package = "spqrp", mustWork = TRUE)
}

#' Load a bundled example data file as a tibble
#'
#' @details
#' The package ships two example CSV files in `inst/extdata/`, both describing
#' a small synthetic cohort intended only for runnable examples and tests:
#'
#' * `example_input_cohort_df.csv` -- mock cohort (30 patients x 2 samples x
#'   5 proteins) in long format with the required columns `Sample_ID`,
#'   `Patient_ID`, `Protein`, `Intensity`.
#' * `example_protein_ranking.csv` -- protein importance ranking aligned with
#'   the mock cohort.
#'
#' The real-cohort protein-importance ranking is provided separately as the
#' lazy-loaded [cohort_a_ranking] dataset: a tibble of `Protein` / `Importance`
#' computed by the pairwise balanced random-forest classifier on plasma cohort
#' "A". It is the built-in default ranking for
#' [perform_distance_evaluation_on_ranked_proteins()] and
#' [optimize_parameters()], and is accessed with `data(cohort_a_ranking)` or
#' `spqrp::cohort_a_ranking` rather than through this function.
#'
#' Use [spqrp_example_path()] if you need the file path instead of the
#' loaded data.
#'
#' @inheritParams spqrp_example_path
#' @return A tibble.
#' @examples
#' spqrp_example_data("input_cohort_df")
#' @export
spqrp_example_data <- function(which = c("input_cohort_df",
                                         "protein_ranking")) {
  tibble::as_tibble(utils::read.csv(spqrp_example_path(which),
                                    stringsAsFactors = FALSE))
}

# Internal: the built-in default protein-importance ranking, fetched
# explicitly from the package namespace (the lazy-loaded `cohort_a_ranking`
# dataset). Used as the fallback ranking by the threshold-based functions
# when the caller supplies neither a data frame nor a CSV path.
default_ranking <- function() {
  get("cohort_a_ranking", envir = asNamespace("spqrp"))
}

# Misc ------------------------------------------------------------------

# Null-coalescing operator: prefer x, fall back to y if x is NULL
`%||%` <- function(x, y) if (is.null(x)) y else x
