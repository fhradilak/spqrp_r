test_that("percentile_cutoff matches numpy.percentile (type-7 quantile)", {
  x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  # numpy.percentile(x, 25) = 3.25 (type 7 linear interpolation)
  expect_equal(percentile_cutoff(x, 25), 3.25, tolerance = 1e-12)
  expect_equal(percentile_cutoff(x, 50), 5.5,  tolerance = 1e-12)
  expect_equal(percentile_cutoff(x, 99), 9.91, tolerance = 1e-12)
})

test_that("get_distances euclidean matches stats::dist reference", {
  df <- data.frame(
    Sample_ID = rep(c("S1", "S2", "S3"), each = 3),
    Patient_ID = rep(c("P1", "P2", "P3"), each = 3),
    Protein = rep(c("Pa", "Pb", "Pc"), times = 3),
    Intensity = c(1, 2, 3,
                  4, 5, 6,
                  7, 8, 9)
  )
  res <- get_distances(df, metric = "euclidean")
  expect_equal(dim(res$distance_matrix), c(3, 3))
  expect_equal(unname(diag(res$distance_matrix)), c(0, 0, 0))
  # |S1 - S2| in euclidean = sqrt(3 * 9) = sqrt(27) (since each protein diff is 3)
  expect_equal(res$distance_matrix["S1", "S2"], sqrt(27), tolerance = 1e-10)
})

test_that("get_distances manhattan matches expected sum-of-abs", {
  df <- data.frame(
    Sample_ID = rep(c("S1", "S2"), each = 3),
    Patient_ID = rep(c("P1", "P2"), each = 3),
    Protein = rep(c("Pa", "Pb", "Pc"), times = 2),
    Intensity = c(1, 2, 3, 4, 5, 6)
  )
  res <- get_distances(df, metric = "manhattan")
  expect_equal(res$distance_matrix["S1", "S2"], 9, tolerance = 1e-10)
})

test_that("get_distances correlation is in [0, 2] and symmetric", {
  set.seed(1)
  df <- expand.grid(
    Sample_ID = c("S1", "S2", "S3", "S4"),
    Protein = c("Pa", "Pb", "Pc", "Pd")
  )
  df$Patient_ID <- substr(df$Sample_ID, 1, 2)
  df$Intensity <- stats::rnorm(nrow(df))
  res <- get_distances(df, metric = "correlation")
  expect_true(isSymmetric(res$distance_matrix, tol = 1e-10))
  expect_true(all(res$distance_matrix >= 0))
  expect_true(all(res$distance_matrix <= 2))
  expect_equal(unname(diag(res$distance_matrix)), rep(0, 4))
})

test_that("get_evaluation_metrics computes a clean confusion matrix", {
  belonging <- tibble::tibble(
    sample1 = c("a", "b"), sample2 = c("a2", "c"),
    patient_id_1 = c("P1", "P2"), patient_id_2 = c("P1", "P3"),
    distance = c(0.1, 0.2)
  )
  not_belonging <- tibble::tibble(
    sample1 = c("d", "e"), sample2 = c("d2", "f"),
    patient_id_1 = c("P4", "P5"), patient_id_2 = c("P4", "P6"),
    distance = c(0.9, 1.0)
  )
  # Exercise the verbose path; suppressMessages keeps testthat output clean.
  m <- suppressMessages(
    get_evaluation_metrics(belonging, not_belonging, quiet = FALSE)
  )
  expect_equal(m$TP, 1L)
  expect_equal(m$FP, 1L)
  expect_equal(m$FN, 1L)
  expect_equal(m$TN, 1L)
  expect_equal(m$Precision, 0.5)
  expect_equal(m$Sensitivity, 0.5)
  expect_equal(m$F1, 0.5)
})

test_that("check_input_data_format fails with informative error on missing column", {
  bad <- data.frame(Sample_ID = "S1", Intensity = 1, Protein = "P")  # missing Patient_ID
  expect_error(check_input_data_format(bad), "Patient_ID")
})

test_that("spqrp_example_path resolves bundled files", {
  p <- spqrp_example_path("input_cohort_df")
  expect_true(file.exists(p))
})

test_that("spqrp_example_data loads bundled cohort with required columns", {
  df <- spqrp_example_data("input_cohort_df")
  expect_true(all(c("Sample_ID", "Patient_ID", "Protein", "Intensity") %in% names(df)))
})
