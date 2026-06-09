test_that("perform_distance_evaluation_on_ranked_proteins returns coherent outputs", {
  df <- spqrp_example_data("input_cohort_df")
  ranking <- spqrp_example_data("protein_ranking")
  # Exercise the verbose code path explicitly (default flipped to TRUE
  # when we moved the package to silent-by-default); suppressMessages
  # keeps the testthat output clean.
  res <- suppressMessages(perform_distance_evaluation_on_ranked_proteins(
    df = df, top_importance_df = ranking,
    metric = "manhattan", p = 0.989, n = 4L, plot = FALSE,
    quiet = FALSE
  ))
  expect_named(res, c("top_importance", "nearest_neighbours", "cutoff",
                       "belonging", "not_belonging", "eval_metrics",
                       "distance_matrix", "plot"))
  expect_true(is.numeric(res$cutoff))
  expect_true(is.matrix(res$distance_matrix))
  expect_equal(res$eval_metrics$TP + res$eval_metrics$FP,
                nrow(res$belonging))
  expect_equal(res$eval_metrics$TN + res$eval_metrics$FN,
                nrow(res$not_belonging))
  # The returned list must NOT have the old `not_be` key -- only `not_belonging`
  expect_false("not_be" %in% names(res))
  # `plot = FALSE` was passed, so $plot is NULL.
  expect_null(res$plot)
})

test_that("optimize_parameters returns a non-empty tibble", {
  df <- spqrp_example_data("input_cohort_df")
  ranking <- spqrp_example_data("protein_ranking")
  out <- suppressMessages(optimize_parameters(df, metric = "manhattan",
                                                top_importance_df = ranking,
                                                range = 3:4,
                                                quiet = FALSE))
  expect_s3_class(out, "tbl_df")
  expect_true(nrow(out) >= 1L)
  expect_true(all(c("n", "FP", "FN", "TP", "TN", "F1") %in% names(out)))
})
