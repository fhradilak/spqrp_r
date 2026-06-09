# Provenance for the `cohort_a_ranking` dataset.
#
# Source: pre-computed pairwise random-forest protein-importance ranking from
# real plasma-proteome cohort "A". Originally shipped as
# inst/extdata/ranked_classification_importance_cohort_a.csv; promoted to a
# documented, lazy-loaded package dataset so it is a citable internal resource
# rather than an anonymous example file.
#
# This script and the raw CSV live in data-raw/, which is listed in
# .Rbuildignore and therefore NOT included in the built/installed package.
#
# Re-run with:  Rscript data-raw/cohort_a_ranking.R   (from the package root)

csv <- "data-raw/ranked_classification_importance_cohort_a.csv"

cohort_a_ranking <- tibble::as_tibble(
  utils::read.csv(csv, stringsAsFactors = FALSE)
)

stopifnot(
  all(c("Protein", "Importance") %in% names(cohort_a_ranking)),
  is.numeric(cohort_a_ranking$Importance)
)

save(
  cohort_a_ranking,
  file = "data/cohort_a_ranking.rda",
  compress = "xz",
  version = 2
)
