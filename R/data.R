#' Protein-importance ranking for plasma cohort "A"
#'
#' A pre-computed protein-importance ranking produced by the pairwise balanced
#' random-forest classifier ([train_pairwise_balanced_rand_forest()]) on a real
#' mass-spectrometry plasma-proteome cohort. It serves as the built-in default
#' ranking for [perform_distance_evaluation_on_ranked_proteins()] and
#' [optimize_parameters()] when the caller supplies neither
#' `top_importance_df` nor `top_importance_path`.
#'
#' @format A [tibble][tibble::tibble] with one row per protein and two columns:
#' \describe{
#'   \item{Protein}{Character. Protein identifier (UniProt accession with gene
#'     suffix, e.g. `"P01861_IGHG4"`).}
#'   \item{Importance}{Numeric. Random-forest importance score; higher means
#'     more discriminative. Rows are ordered from most to least important.}
#' }
#'
#' @source Pairwise balanced random-forest importances computed on plasma
#'   cohort "A", derived from mass-spectrometry plasma-proteome measurements.
#'
#' @seealso [perform_distance_evaluation_on_ranked_proteins()],
#'   [optimize_parameters()], [retrieve_ranking()]
#'
#' @examples
#' head(cohort_a_ranking)
"cohort_a_ranking"
