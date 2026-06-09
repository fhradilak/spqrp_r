# Pin the silent-by-default contract introduced for CRAN-readiness.
#
# Two directions per user-facing entry point:
#   1. Default call (`quiet` unset): zero informational messages, zero
#      stdout chatter, and an invisibly-returned value (so the REPL
#      doesn't auto-print the list at the end either).
#   2. Explicit `quiet = FALSE`: at least one message line emitted.
#
# Legitimate `cli_warn` calls (data-loss / deprecation) are out of scope
# here -- they fire regardless of `quiet` and are tested separately.

# Helper: capture messages emitted by `expr`, return the count.
n_messages <- function(expr) {
  length(capture.output(force(expr), type = "message"))
}

# Helper: is the value invisible?
is_invisible <- function(expr) {
  !withVisible(expr)$visible
}

test_that("by_isolation_forest is silent by default", {
  df <- spqrp_example_data("input_cohort_df")
  expect_equal(
    n_messages(by_isolation_forest(df, impute_median = TRUE)),
    0L
  )
})

test_that("by_isolation_forest emits status when quiet = FALSE", {
  df <- spqrp_example_data("input_cohort_df")
  expect_gt(
    n_messages(by_isolation_forest(df, impute_median = TRUE, quiet = FALSE)),
    0L
  )
})

test_that("remove_outlier_samples is silent by default and returns invisibly", {
  df <- spqrp_example_data("input_cohort_df")
  expect_equal(n_messages(remove_outlier_samples(df)), 0L)
  expect_true(is_invisible(remove_outlier_samples(df)))
})

test_that("perform_distance_evaluation_on_ranked_proteins is silent by default", {
  df <- spqrp_example_data("input_cohort_df")
  ranking <- spqrp_example_data("protein_ranking")
  expect_equal(
    n_messages(perform_distance_evaluation_on_ranked_proteins(
      df = df, top_importance_df = ranking, n = 4L, p = 0.5, plot = FALSE
    )),
    0L
  )
})

test_that("perform_distance_evaluation_on_ranked_proteins emits when quiet=FALSE", {
  df <- spqrp_example_data("input_cohort_df")
  ranking <- spqrp_example_data("protein_ranking")
  expect_gt(
    n_messages(perform_distance_evaluation_on_ranked_proteins(
      df = df, top_importance_df = ranking, n = 4L, p = 0.5, plot = FALSE,
      quiet = FALSE
    )),
    0L
  )
})

test_that("run_clustering is silent by default", {
  df <- spqrp_example_data("input_cohort_df")
  ranking <- spqrp_example_data("protein_ranking")
  # `cluster_samples_iteratively` emits a legitimate `cli_warn` about
  # dropped samples -- intentionally unconditional. Suppress only that
  # warning channel; messages must still be zero.
  expect_equal(
    suppressWarnings(n_messages(
      run_clustering(df = df, ranking = ranking,
                      n_neighbors = 1L, max_component_size = 2L,
                      method = "PCA", n = 4L)
    )),
    0L
  )
})

test_that("train_with_normalise is silent by default", {
  df <- spqrp_example_data("input_cohort_df")
  expect_equal(
    n_messages(train_with_normalise(df, plate_corrected = FALSE,
                                      outlier_removal = FALSE)),
    0L
  )
})

test_that("train_with_normalise emits status when quiet = FALSE", {
  df <- spqrp_example_data("input_cohort_df")
  expect_gt(
    n_messages(train_with_normalise(df, plate_corrected = FALSE,
                                      outlier_removal = FALSE,
                                      quiet = FALSE)),
    0L
  )
})
