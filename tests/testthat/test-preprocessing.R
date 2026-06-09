test_that("filter_by_occurrence removes infrequent proteins", {
  df <- data.frame(
    Sample_ID = c("S1", "S2", "S3", "S1", "S2", "S3"),
    Patient_ID = c("P1", "P2", "P3", "P1", "P2", "P3"),
    Protein   = c("A",  "A",  "A",  "B",  "B",  "B"),
    Intensity = c(1, 2, 3, NA, NA, 5)  # B present in 1/3 samples
  )
  out <- filter_by_occurrence(df, cutoff = 0.5)
  expect_true("A" %in% out$Protein)
  expect_false("B" %in% out$Protein)
})

test_that("log_transform and revert_log_transform invert", {
  df <- data.frame(
    Sample_ID = "S1", Patient_ID = "P1", Protein = "A",
    Intensity = c(1, 2, 4, 8)
  )
  back <- revert_log_transform(log_transform(df))
  expect_equal(back$Intensity, c(1, 2, 4, 8))
})

test_that("normalize_medianintensity centers per-sample median", {
  df <- data.frame(
    Sample_ID = rep(c("S1", "S2"), each = 5),
    Patient_ID = rep(c("P1", "P2"), each = 5),
    Protein = rep(LETTERS[1:5], times = 2),
    Intensity = c(1:5, 11:15) + 0.0
  )
  res <- normalize_medianintensity(df, plot = FALSE)
  per_sample_med <- tapply(res$data$Intensity, res$data$Sample_ID, stats::median)
  # After per-sample log-space normalization the per-sample medians should
  # cluster close together (within numerical noise of the global median).
  expect_lt(diff(range(per_sample_med)), 1e-9)
})
