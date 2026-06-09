#' spqrp: Sample Provenance Quality Resolver in Proteomics
#'
#' Detects sample-provenance inconsistencies in MS-based plasma-proteome
#' cohorts via pairwise distance, threshold-based classification, iterative
#' clustering, and a pairwise random-forest classifier for protein importance
#' ranking. Native R port of the Python package of the same name.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang .data %||% :=
#' @importFrom stats cmdscale cor dist lm median predict quantile residuals sd setNames
#' @importFrom utils combn read.csv
## usethis namespace: end
NULL
