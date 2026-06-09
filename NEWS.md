# spqrp 0.1.0

Initial CRAN-targeted release. Native R port of the Python `spqrp`
package.

## Features

* **Clustering pipeline** (`run_clustering()`) -- kNN graph in a
  PCA/UMAP/MDS embedding, optional iterative split of large
  components, ggplot visualization with patient-hue colouring and
  legend.
* **Threshold-based evaluation** (`perform_distance_evaluation_on_ranked_proteins()`)
  -- pairwise sample classification from a percentile cutoff on
  pairwise distances, with FN/FP/percentile-overlay histogram.
* **Protein ranking** (`train_with_normalise()`) -- pairwise random-
  forest classifier with three selectable backends; `randomForest`
  is the default (closest behaviour to Python's
  `imblearn.BalancedRandomForestClassifier`). Importance values are
  normalised to sum to 1.0, matching sklearn's
  `clf.feature_importances_` convention.
* **Isolation-forest outlier filtering** (`remove_outlier_samples()`)
  -- pure-R via the `solitude` package; default `outlier_threshold`
  calibrated empirically for solitude's anomaly-score scale.

## Verbosity

All functions are silent by default. Pass `quiet = FALSE` to any
function that emits status output to see progress messages, per-call
summaries, save-path hints, and cluster listings. Warnings about
genuine data issues -- e.g. samples dropped from analysis -- fire
regardless of `quiet`.

## Documentation

* See [`articles/numerical-divergence.md`](https://github.com/fhradilak/spqrp_r/blob/main/articles/numerical-divergence.md) for known cross-language
  divergences (UMAP, random-forest backends, isolation-forest
  scales, MDS solvers) and recommendations for cross-language
  comparison.
* See `vignette("spqrp-mock-data")` for a worked example on a small
  bundled cohort.
