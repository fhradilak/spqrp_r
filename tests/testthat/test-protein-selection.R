test_that("train_with_normalise runs end-to-end with the ranger backend", {
  df <- spqrp_example_data("input_cohort_df")
  res <- suppressMessages(train_with_normalise(
    df, plate_corrected = FALSE,
    outlier_removal = FALSE,
    classifier_backend = "ranger",
    method = "F1",
    quiet = FALSE   # exercise verbose path; suppressMessages keeps output clean
  ))
  expect_s3_class(res, "spqrp_train")
  expect_true(length(res$results_dict$feature_importances) > 0L)

  rk <- retrieve_ranking(res)
  expect_s3_class(rk, "tbl_df")
  expect_true(all(c("Protein", "Importance") %in% names(rk)))
  expect_false(any(startsWith(rk$Protein, "diff_")))
})

test_that("randomForest backend works (the default)", {
  df <- spqrp_example_data("input_cohort_df")
  res <- suppressMessages(train_with_normalise(
    df, plate_corrected = FALSE,
    outlier_removal = FALSE,
    classifier_backend = "randomForest",
    method = "F1",
    quiet = FALSE
  ))
  expect_s3_class(res, "spqrp_train")
  expect_equal(res$classifier_backend, "randomForest")
})

test_that("themis_smote backend works if packages are available", {
  skip_if_not_installed("themis")
  skip_if_not_installed("recipes")
  df <- spqrp_example_data("input_cohort_df")
  res <- suppressMessages(train_with_normalise(
    df, plate_corrected = FALSE,
    outlier_removal = FALSE,
    classifier_backend = "themis_smote",
    method = "F1",
    quiet = FALSE
  ))
  expect_s3_class(res, "spqrp_train")
  expect_equal(res$classifier_backend, "themis_smote")
})

test_that("get_threshold picks a threshold in (0, 1) for separable data", {
  y <- c(rep(0L, 50), rep(1L, 50))
  p <- c(stats::runif(50, 0, 0.3), stats::runif(50, 0.7, 1))
  thr <- get_threshold(y, p, method = "F1")
  expect_true(thr$threshold > 0 && thr$threshold < 1)
  expect_equal(length(thr$y_pred_adjusted), length(y))
})
