test_that("public functions abort when required columns are missing", {
  bad <- data.frame(Sample = "S1", Patient = "P1", Protein = "P", Intensity = 1)
  expect_error(check_input_data_format(bad), "Sample_ID|Patient_ID")
  expect_error(
    calculate_pairwise_distances(
      top_importance = data.frame(Protein = "P", Importance = 1),
      n = 1, df = bad, metric = "euclidean"
    ),
    "Sample_ID|Patient_ID"
  )
})

test_that("filter_by_occurrence rejects out-of-range cutoffs", {
  df <- data.frame(Sample_ID = "S1", Patient_ID = "P1",
                    Protein = "A", Intensity = 1)
  expect_error(filter_by_occurrence(df, cutoff = -0.1), "in \\[0, 1\\]")
  expect_error(filter_by_occurrence(df, cutoff = 1.5),  "in \\[0, 1\\]")
})
