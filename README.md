# spqrp

[![License: GPL-3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)

**Sample Provenance Quality Resolver in Proteomics** — native R port of
the [Python `spqrp` package](https://github.com/fhradilak/spqrp).

Recent advancements in MS technology and lab methods opened the door for
large-scale proteomics but also led to a growing concern regarding sample
mix-ups. `spqrp` helps you evaluate whether sample data is safe for
further analysis by clustering samples and flagging probable mix-ups,
uncertain assignments, and outliers.

## Install

```r
# install.packages("remotes")
remotes::install_github("fhradilak/spqrp_r")
```

No Python install needed — this is a native R port.

## Input data format

A long-format data frame with these columns:

| Column      | Description                       |
|-------------|-----------------------------------|
| `Sample_ID` | Unique sample identifier          |
| `Patient_ID`| Patient identifier                |
| `Protein`   | Protein name/identifier           |
| `Intensity` | Numeric intensity value           |

Optionally a protein ranking with `Protein` and `Importance` columns.
If you don't supply one, the package uses a precomputed ranking from a
plasma cohort (`spqrp_example_data("ranking_cohort_a")`).

## Quick start

```r
library(spqrp)

df      <- spqrp_example_data("input_cohort_df")
ranking <- spqrp_example_data("protein_ranking")

# Clustering: build kNN graph, split big components, visualise
res <- run_clustering(
  df = df, ranking = ranking,
  n_neighbors = 1L,
  max_component_size = 2L,
  metric = "manhattan",
  method = "UMAP"   # or "PCA" / "MDS"
)

res$cluster_assignments       # sample -> cluster ID
res$uncertain_samples         # likely missing connections
res$error_candidate_samples   # likely sample mix-ups
res$plot                      # ggplot object
```

## Verbose output

All `spqrp` functions are silent by default --- no progress messages,
no per-call summaries. If you want progress and diagnostic prints
(which sample IDs got flagged, what cutoff was picked, how many
proteins overlapped between ranking and data, etc.) pass
`quiet = FALSE` to any function that emits status output:

```r
remove_outlier_samples(df, quiet = FALSE)         # prints flagged Sample_IDs
run_clustering(df, ranking, n_neighbors = 1,
                max_component_size = 3, quiet = FALSE)  # prints save-path hint,
                                                          # cluster listing,
                                                          # transitive metrics
perform_distance_evaluation_on_ranked_proteins(
  df, top_importance_df = ranking, quiet = FALSE
)                                                 # prints real-protein count
```

Warnings about genuine data issues --- e.g. samples dropped because
they lack measurements for any of the top-ranked proteins --- fire
**regardless** of `quiet`, because they signal a real problem you
need to see.

## Three random-forest backends for protein ranking

If you don't have a precomputed ranking, train one. The Python package
uses `imblearn.BalancedRandomForestClassifier`; this R port exposes
three substitute backends so you can pick the tradeoff that fits:

```r
results <- train_with_normalise(
  df,
  classifier_backend = "randomForest"  # default — closest to imblearn's BalancedRF
  # classifier_backend = "ranger"        # faster, class.weights on impurity
  # classifier_backend = "themis_smote"  # SMOTE rebalance + ranger
)

new_ranking <- retrieve_ranking(results)
```

See [`articles/numerical-divergence.md`](https://github.com/fhradilak/spqrp_r/blob/HEAD/articles/numerical-divergence.md)
for when to pick each.

## Threshold-based evaluation

```r
result <- perform_distance_evaluation_on_ranked_proteins(
  df = df,
  top_importance_df = ranking,
  metric = "manhattan",
  p = 0.989,
  n = 20L
)
result$cutoff
result$eval_metrics[c("TP", "FP", "FN", "TN", "Precision", "Sensitivity", "F1")]
```

`optimize_parameters()` sweeps `n` and the percentile cutoff to find
optimal values for your dataset.

## Preprocessing

Optional helpers that mirror the Python pipeline:

```r
df_pp <- df |>
  log_transform() |>
  filter_by_occurrence(cutoff = 0.7)

norm <- normalize_medianintensity(df_pp, plot = FALSE)
df_pp <- norm$data

# If your data has a `plate` column:
df_pp <- plate_correct_residuals_by_protein(df_pp)
```

## Function reference

| Function | Purpose |
|---|---|
| `run_clustering()` | End-to-end clustering pipeline |
| `cluster_samples_iteratively()` | Build kNN graph + 2D embedding |
| `plot_distances_neighbours_with_coloring_hue()` | Heavy clustering visualization |
| `perform_distance_evaluation_on_ranked_proteins()` | Threshold-based pairwise classification |
| `optimize_parameters()` | Grid-search optimal `n` and percentile |
| `calculate_pairwise_distances()` | Distance matrix on top-`n` proteins |
| `train_with_normalise()` | Full ranking pipeline (filter → normalize → RF) |
| `retrieve_ranking()` | Extract ranked proteins from a trained model |
| `train_pairwise_balanced_rand_forest()` | Pairwise RF (3 backends) |
| `get_threshold()` | ROC/F1/Youden/MinFP threshold selection |
| `get_distances()`, `get_nearest_neighbours()` | Distance + kNN helpers |
| `get_sample_relations_by_cutoff()`, `get_evaluation_metrics()` | Cutoff → metrics |
| `percentile_cutoff()` | numpy.percentile-equivalent |
| `filter_by_occurrence()`, `log_transform()`, `revert_log_transform()`, `normalize_medianintensity()`, `plate_correct_residuals_by_protein()` | Preprocessing |
| `by_isolation_forest()`, `by_isolation_forest_plot()`, `remove_outlier_samples()` | Outlier detection (Isolation Forest). `contamination = 0.1` for sklearn-like behaviour. |
| `spqrp_example_data()` | Access bundled example CSVs |
| `check_input_data_format()` | Validate required columns |

## Migrating from the Python version

The R API mirrors the Python one — function names are identical
snake_case. Outputs are R named lists (which work just like Python
dicts: `res$cluster_assignments`).

Because the underlying numerical libraries differ (`uwot` vs
`umap-learn`, `ranger` vs `imblearn`, `solitude` (wrapping `ranger`)
vs sklearn's IsolationForest), exact numbers can drift across runs
even with matched seeds. See
[`articles/numerical-divergence.md`](https://github.com/fhradilak/spqrp_r/blob/HEAD/articles/numerical-divergence.md)
for which outputs are bit-exact, which match up to rotation/reflection,
and which are only equivalent in expectation.

## License

GPL-3
