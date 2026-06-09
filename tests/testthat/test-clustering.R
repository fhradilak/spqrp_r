# Tests exercise the verbose code path (`quiet = FALSE`) but wrap calls
# in `suppressMessages()` and `suppressWarnings()` to keep testthat
# output clean. `cluster_samples_iteratively` emits a legitimate
# `cli_warn` about samples dropped from the kNN graph -- that's
# unconditional by design, so we suppress warnings here too.

test_that("cluster_samples_iteratively returns a 2D embedding and igraph", {
  df <- spqrp_example_data("input_cohort_df")
  ranking <- spqrp_example_data("protein_ranking")
  result <- suppressMessages(calculate_pairwise_distances(
    top_importance = ranking, n = 5L, df = df, metric = "manhattan",
    quiet = FALSE
  ))
  out <- suppressMessages(suppressWarnings(
    cluster_samples_iteratively(result, df, method = "PCA",
                                 n_neighbors = 1L,
                                 max_component_size = 2L,
                                 quiet = FALSE)
  ))
  expect_s3_class(out$G, "igraph")
  expect_equal(ncol(out$coords_2d), 2L)
  expect_equal(nrow(out$coords_2d), length(unique(df$Sample_ID)))

  comps <- igraph::components(out$G)
  expect_true(all(comps$csize <= 2L))
})

test_that("run_clustering with PCA returns the expected list shape", {
  df <- spqrp_example_data("input_cohort_df")
  ranking <- spqrp_example_data("protein_ranking")
  res <- suppressMessages(suppressWarnings(
    run_clustering(df = df, ranking = ranking,
                    n_neighbors = 1L, max_component_size = 2L,
                    plot_name = "test", method = "PCA",
                    quiet = FALSE)
  ))
  expect_named(res, c("result_filtered", "G", "cluster_assignments",
                       "transitive_results", "uncertain_samples",
                       "error_candidate_samples", "saved_path", "plot"))
  expect_true(length(res$cluster_assignments) == length(unique(df$Sample_ID)))
  expect_s3_class(res$plot, "ggplot")
})
