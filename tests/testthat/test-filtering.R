make_outlier_df <- function() {
  set.seed(0)
  proteins <- paste0("P", 1:8)
  normal_samples <- paste0("S", 1:9)
  outlier_sample <- "Soutlier"
  rows <- list()
  for (s in normal_samples) {
    rows[[length(rows) + 1L]] <- data.frame(
      Sample_ID = s, Patient_ID = paste0("P", s),
      Protein = proteins,
      Intensity = stats::rnorm(length(proteins), mean = 0, sd = 0.1)
    )
  }
  rows[[length(rows) + 1L]] <- data.frame(
    Sample_ID = outlier_sample, Patient_ID = "POUT",
    Protein = proteins,
    Intensity = stats::rnorm(length(proteins), mean = 50, sd = 0.1)
  )
  list(df = do.call(rbind, rows), outlier_sample = outlier_sample)
}

# Tests exercise the verbose code path (`quiet = FALSE`) but wrap calls
# in `suppressMessages()` to keep testthat output clean. The package's
# default is `quiet = TRUE`; flipping to FALSE here ensures the cli_inform
# branches in the source are still hit at least once per CI run.

test_that("by_isolation_forest flags an obvious outlier", {
  fixture <- make_outlier_df()
  res <- suppressMessages(by_isolation_forest(fixture$df, impute_median = TRUE,
                                                outlier_threshold = 0.55,
                                                quiet = FALSE))
  expect_true(fixture$outlier_sample %in% res$outlier_list)
})

test_that("contamination = numeric flags exactly that fraction of samples", {
  fixture <- make_outlier_df()
  # 10 samples, contamination = 0.1 => exactly 1 outlier flagged
  res <- suppressMessages(by_isolation_forest(fixture$df, impute_median = TRUE,
                                                contamination = 0.1,
                                                quiet = FALSE))
  expect_equal(length(res$outlier_list), 1L)
  expect_true(fixture$outlier_sample %in% res$outlier_list)

  res2 <- suppressMessages(by_isolation_forest(fixture$df, impute_median = TRUE,
                                                 contamination = 0.3,
                                                 quiet = FALSE))
  expect_equal(length(res2$outlier_list), 3L)
})

test_that("contamination = 0 flags nothing", {
  fixture <- make_outlier_df()
  res <- suppressMessages(by_isolation_forest(fixture$df, impute_median = TRUE,
                                                contamination = 0,
                                                quiet = FALSE))
  expect_equal(length(res$outlier_list), 0L)
})

test_that("contamination rejects bad input", {
  fixture <- make_outlier_df()
  expect_error(
    by_isolation_forest(fixture$df, impute_median = TRUE,
                         contamination = 1.5, quiet = FALSE),
    "contamination"
  )
  expect_error(
    by_isolation_forest(fixture$df, impute_median = TRUE,
                         contamination = "wrong", quiet = FALSE),
    "contamination"
  )
})

test_that("remove_outlier_samples removes flagged samples", {
  fixture <- make_outlier_df()
  out <- suppressMessages(remove_outlier_samples(fixture$df,
                                                   contamination = 0.1,
                                                   quiet = FALSE))
  expect_named(out, c("df", "anomaly_df", "outlier_list", "anomaly_plot"),
               ignore.order = TRUE)
  expect_false(fixture$outlier_sample %in% out$df$Sample_ID)
  skip_if_not_installed("plotly")
  expect_s3_class(out$anomaly_plot, "plotly")
})

