# Numerical divergence: spqrp R vs Python

This package is a native-R port of the Python package `spqrp`. Most of
the algorithms are 1:1 ports, but a few rely on libraries that differ
across languages. This document records **where** the outputs diverge,
**why**, **how much**, and **what to do** if you need parity.

## TL;DR

* **Identical to ≤ 1e-12**
  * Distance metrics (`euclidean`, `manhattan`, `minkowski` /
    `fractional`, `correlation`)
  * Percentile cutoffs (`percentile_cutoff` uses type-7 quantile, the
    `numpy.percentile` default)
  * Pairwise evaluation metrics (TP/FP/FN/TN, precision, recall, F1)
    once the cutoff is fixed
  * OLS plate-correction residuals (`stats::lm` ≡ `statsmodels.OLS`)
  * Connected components (`igraph` ≡ `networkx`)
* **Equivalent up to rotation/reflection**
  * 2D embeddings (PCA / MDS / UMAP)
* **Equivalent in expectation, not point-for-point**
  * Random-forest feature rankings
  * Isolation-forest outlier flags
  * Train/test splits when seeded
* **Never compare raw**
  * UMAP coordinates between R and Python
  * RF probability scores between R and Python

For cross-language comparisons, compare downstream *qualitative* results
(cluster *membership*, outlier sets, FP/FN classifications) — not
embedding coordinates or raw probabilities.

## 1. UMAP — `uwot` (R) vs `umap-learn` (Python)

**What it does.** Non-linear dimensionality reduction; in spqrp we feed
it a precomputed distance matrix and use the 2D coordinates only for
visualisation.

**Python.** `umap-learn` (McInnes et al.). Adam-like optimizer,
spectral initialization with one Laplacian normalization.

**R.** `uwot::umap(X = as.dist(d), metric = "precomputed", seed = ...)`.
Stochastic gradient descent on the cross-entropy objective, spectral
initialization with a *different* Laplacian normalization.

**Why they diverge.** Different optimizer schedules and slightly
different defaults for `set_op_mix_ratio` and minimum-distance smoothing.
Same seed produces deterministic but library-specific results.

**Empirical magnitude.** On the bundled 30-sample mock cohort the
Procrustes-aligned RMSE between R and Python UMAP coords is on the
order of 0.1–0.3 distance units. Cluster topology (which samples are
neighbours) is preserved.

**What to do if you care about parity.**
* Use `method = "PCA"` (closed-form, identical math both sides).
* Or stick with UMAP and compare *graph membership*, not coords:

```r
# Suppose r_assignments / py_assignments are two named integer vectors of
# cluster IDs over the same samples.
ari_score <- function(r, py) {
  samples <- intersect(names(r), names(py))
  tab <- table(r[samples], py[samples])
  n <- sum(tab)
  a <- sum(choose(rowSums(tab), 2))
  b <- sum(choose(colSums(tab), 2))
  ab <- a * b / choose(n, 2)
  (sum(choose(tab, 2)) - ab) / (0.5 * (a + b) - ab)
}
```

ARI of 1 means identical membership (regardless of cluster IDs).

## 2. MDS — `stats::cmdscale` (R) vs `sklearn.manifold.MDS` (Python)

**Python.** sklearn's `MDS(metric=True)` runs iterative SMACOF.

**R.** `stats::cmdscale` is classical (Torgerson) MDS — the closed-form
double-centred eigendecomposition. This is mathematically equivalent
to the `J · (-0.5 · d²) · J` then PCA construction in the Python source
(`core.py`), just expressed in one cmdscale call.

**Why they diverge.** Iterative SMACOF (Python) optimizes a stress
function via majorization; cmdscale (R) solves it in closed form. They
agree exactly when the distance matrix is Euclidean and diverge
slightly on non-Euclidean (e.g. correlation) inputs.

**What to do for sklearn-faithful behaviour.**

```r
# Optional dependency
# install.packages("smacof")
res <- cluster_samples_iteratively(result_filtered, df,
                                    method = "MDS",
                                    mds_backend = "smacof")
```

## 3. Training pipeline — `randomForest`/`ranger` (R) vs `imblearn.BalancedRandomForestClassifier` (Python)

### Why bit-identical training results are impossible

`train_with_normalise()` runs the same algorithm in both ports: filter
proteins by occurrence, log-transform, split by patient, normalise per
sample, optionally plate-correct, isolate outliers, build pairwise
features (`X[i] - X[j]`, label = same patient?), fit a random-forest
classifier, extract Gini-impurity importance, strip the `diff_` prefix
to get a per-protein ranking. The *pipeline shape* is identical. But
at the implementation level several steps use libraries that cannot
share a common RNG, and the random forest itself is not the same
algorithm:

| Step | Python | R |
|---|---|---|
| Train/test split | `sklearn.model_selection.train_test_split(random_state=42)` | base R `sample.int()` after `set.seed(42)` |
| Outlier filtering | `sklearn.ensemble.IsolationForest` | `solitude::isolationForest` (see Section 4) |
| Pairwise NaN distance | `sklearn.metrics.nan_euclidean_distances` | hand-rolled with same formula (see [`compute_nan_euclidean_distance`](../R/protein_selection.R)) ✓ matches |
| Random forest | **`imblearn.BalancedRandomForestClassifier`** — per-tree balanced bootstrap, sklearn RNG | **`randomForest` (default), `ranger`, or `themis_smote`** — each different |
| Importance metric | `clf.feature_importances_` (impurity, normalised to sum to 1.0) | `randomForest::importance(..., type = 2)` (MeanDecreaseGini) — closest to sklearn; `ranger`'s `variable.importance` rolls class-weighting into the value. The R port **post-normalises the result to sum to 1.0** so the numeric scale matches sklearn. |

The first two rows mean: even with `random_state = 42` and `set.seed(42)`,
Python and R draw a different sequence of indices and split patients
into different train/test groups. There is no seed value that recovers
parity, because the Mersenne-Twister implementations are not
portable across the two languages even when starting from the same
numeric seed.

The fourth row is the bigger one. `BalancedRandomForestClassifier`
**undersamples the majority class within each tree's bootstrap draw**.
No R Random Forest package implements that exact contract. The closest
match is `randomForest::randomForest(sampsize = c(min_n, min_n),
strata = labels)`, which produces a per-tree balanced bootstrap with
the same minority-size cap. `ranger` does not bootstrap-balance at all
— it re-weights the impurity calculation by class frequency, which is
qualitatively different and produces visibly different importance
distributions.

The fifth row is subtler but matters when comparing raw importance
*values*. sklearn normalises `feature_importances_` so the vector sums
to `1.0`; both R backends return the un-normalised raw impurity
decrease (`randomForest`'s `MeanDecreaseGini` and `ranger`'s
`variable.importance`). The R port therefore **normalises the
importance vector to sum to `1.0` after extracting it from the
classifier**, so the values you see in `retrieve_ranking()` use the
same scale as sklearn's output (typical bar height: `0.001`–`0.30`,
not `5.2`–`50.0+`). Rank order is preserved by the normalisation.

### What we did to make the defaults match in spirit

1. **Default backend switched to `randomForest`.** Previously `"ranger"`
   for speed; now `"randomForest"` because per-tree balanced bootstrap
   reproduces sklearn's `BalancedRandomForestClassifier` mechanism. Pass
   `classifier_backend = "ranger"` explicitly if you need the speed of
   the C++ implementation and accept a wider divergence from Python.
2. **`randomForest` promoted to `Imports`.** It was a `Suggests` (since
   `"ranger"` was the default). It is now a hard dependency — installing
   the package installs `randomForest` automatically.
3. **NaN-aware Euclidean distance already matches sklearn's formula.**
   [`compute_nan_euclidean_distance`](../R/protein_selection.R)
   implements `sqrt( n / m * sum_finite (x_i - y_i)^2 )` exactly as
   `sklearn.metrics.nan_euclidean_distances` does. No change needed.
4. **Importance values normalised to sum to 1.0.** Both
   `randomForest::importance(type = 2)` and `ranger$variable.importance`
   return un-normalised impurity decreases; the R port post-processes
   them so the returned vector sums to 1.0 across features, matching
   sklearn's `clf.feature_importances_` convention. Rank order is
   preserved, so any downstream code using `arrange(desc(Importance))`
   or `head(rk, n)` is unaffected. Code that hard-coded an absolute
   threshold like `Importance > 1` against the old un-normalised scale
   needs updating; use ranks or percentile cutoffs instead.
5. **Train/test split is a documented escape hatch.** Pass
   `train_individuals` and `test_individuals` to `train_with_normalise`
   to bypass the language-local `sample()`/`train_test_split` RNG.
   Recommended workflow for cross-language comparison: dump the
   per-patient train/test assignment from one side to JSON, load it on
   the other, hand both ports the same split (see "Practical guidance"
   below).

### What still doesn't match, even after all of the above

* **Tree topologies.** Even with `randomForest` and a fixed split, the
  R and Python RFs draw different per-tree feature/threshold pairs
  because the RNG streams remain incompatible. Rankings agree on rank,
  not on numeric importance.
* **Absolute importance magnitudes.** Always compare **ranks** (top-N
  membership, Spearman rank correlation), never raw importance values.
* **Borderline-importance proteins.** Two proteins with nearly identical
  importance in Python may swap positions in R (and vice versa). This
  is the largest residual divergence.
* **Upstream filtering composition.** If `outlier_removal = TRUE`,
  whichever samples each port flags affects train/test composition and
  cascades into the RF. See Section 4 for the outlier-filtering
  divergence.

### Practical guidance

```r
# 1. Default — recommended. Closest behaviour to Python:
results <- train_with_normalise(df)
retrieve_ranking(results)

# 2. Faster backend (older default) — useful on large cohorts where
#    randomForest is the bottleneck. Diverges more from Python:
results_fast <- train_with_normalise(df, classifier_backend = "ranger")

# 3. Heavy class imbalance — synthesise minority via SMOTE first:
results_smote <- train_with_normalise(df, classifier_backend = "themis_smote")

# 4. Compare two rankings: top-20 overlap and Spearman rank correlation
top_n_overlap <- function(rank_a, rank_b, n = 20) {
  length(intersect(head(rank_a$Protein, n), head(rank_b$Protein, n)))
}
spearman <- function(rank_a, rank_b) {
  joined <- merge(rank_a, rank_b, by = "Protein",
                  suffixes = c("_a", "_b"))
  stats::cor(joined$Importance_a, joined$Importance_b, method = "spearman")
}

# 5. Pin the train/test split across languages.
#    In Python:
#        json.dump({"train": train_ids, "test": test_ids}, open("split.json","w"))
#    In R:
#        split <- jsonlite::fromJSON("split.json")
#        results <- train_with_normalise(df,
#                       train_individuals = split$train,
#                       test_individuals  = split$test)
#    This removes split RNG as a source of divergence so any remaining
#    ranking difference is purely the RF backend.
```

The top-`n` proteins agree across backends and across languages much
more often than the absolute importance numbers do. **Always compare
via Jaccard or Spearman rank correlation on the top 10-20 proteins, not
raw importance values.**

## 4. Isolation Forest — `solitude` (R) vs `sklearn.ensemble.IsolationForest`

### Why the two cannot produce identical outlier lists

R uses [`solitude::isolationForest`](https://github.com/talegari/solitude),
a pure-R port of the original Liu et al. (2008) Isolation Forest that
wraps `ranger`. Python uses `sklearn.ensemble.IsolationForest`, scikit-
learn's bespoke C implementation of the same algorithm. They share the
*idea* — build random trees that isolate samples and compute the
anomaly score `2^(-E[h(x)]/c(n))` — but **at every parameter level
below the public API they make different concrete choices**, none of
which the user can reconcile:

| Aspect | sklearn | solitude (via ranger) |
|---|---|---|
| RNG library | NumPy / sklearn's C RNG stream | ranger's C++ `std::mt19937` |
| `seed` semantics | `random_state = 0` is a valid deterministic seed | `seed = 0` means *"ignore seed, randomise"* — we pass `1L` |
| Features per split (`max_features` / `mtry`) | `1` — sklearn's `IsolationForest` hard-codes `max_features=1` on the internal `ExtraTreeRegressor`, so each split sees exactly one random feature | `ncol - 1` (hardcoded inside solitude's `fit`; not user-settable) |
| Split rule | The one random feature × one random threshold drawn from the **node-local** [min, max] is used **as is** — no quality test (truly random) | For each of the `ncol - 1` candidate features, draw one random threshold from the **node-local** [min, max] and pick the **best** of those `ncol - 1` random (feature, threshold) pairs by variance decrease on a random-permutation target |
| `max_depth` | `ceil(log2(max_samples))` | `ceil(log2(sample_size))` (same formula ✓) |
| `bootstrap` | `False` | `False` ✓ |
| `min_samples_leaf` | `1` | `1` ✓ |

The first two rows mean: even with identical data, identical sample
counts and identical formulas, the two libraries draw a **completely
different sequence of random feature/threshold pairs**. There is no
seed value that you could pass to both to recover the same forest.

The next two rows are deeper. Both libraries draw split thresholds
from the **node-local** range of the candidate feature — verified by
reading ranger's `TreeRegression.cpp::findBestSplitValueExtraTrees`,
which calls `data->getMinMaxValues(min, max, sampleIDs, varID,
start_pos[nodeID], end_pos[nodeID])` and then `uniform_real_distribution(min, max)`,
exactly matching sklearn's `RandomSplitter`. So "where the threshold
comes from" is **not** the divergence. The real divergence is:

* sklearn evaluates **one** random (feature, threshold) pair per split
  and accepts it unconditionally. The tree is genuinely random — no
  quality criterion.
* ranger (via solitude) evaluates **`ncol - 1`** random (feature, threshold)
  pairs per split and keeps the **best** one by variance decrease on
  a dummy regression target `yy = sample.int(n)` (a random permutation).
  Even though the dummy target is meaningless for isolation, picking
  "best of `ncol - 1` random splits" still pushes solitude's trees toward
  more informative splits than pure sklearn-style randomness.

In short: solitude's trees are *less* random than sklearn's because
solitude inherited ranger's extratrees machinery, which is designed
for supervised problems. Path lengths come out slightly shorter on
average → anomaly scores come out slightly higher → the same numeric
cutoff has a different meaning.

The cumulative effect: **solitude's anomaly-score distribution is
systematically shifted upward versus sklearn's**. On clean cohorts,
solitude's inlier scores cluster between roughly `0.55` and `0.60`,
where sklearn's would cluster around `0.40`–`0.50`. A given numeric
threshold therefore does not mean the same thing in the two ports.

### What we did to make the defaults match in spirit

Even though the score values cannot match bit-for-bit, the **business
outcome** can: "flag near-zero samples on clean data, flag the
genuinely anomalous samples on dirty data." We tuned the R defaults so
that calling `by_isolation_forest(df)` or `remove_outlier_samples(df)`
with **no arguments** reproduces sklearn's "contamination=auto"
behaviour qualitatively:

1. **Same engine family.** We switched from `isotree`
   (a different IF variant — Cortes' extended IF) to `solitude`, which
   implements the same Liu et al. algorithm sklearn uses. This shrinks
   the systematic gap dramatically; everything below addresses what
   remains.
2. **Matched parameters where possible:**
   `sample_size = nrow(mat) %/% 2L` (= Python's
   `max_samples = len // 2`), `num_trees = 100`,
   single-threaded (`nproc = 1L` / Python `n_jobs = 1` for
   determinism), median imputation per column (matches pandas
   `fillna(df.median())`).
3. **Forced deterministic row and column order** in `long_to_wide()`
   (codepoint sort), so the subsample drawn by ranger is reproducible
   across runs and matches what pandas would feed sklearn.
4. **Re-calibrated the `outlier_threshold` default empirically.**
   sklearn's `contamination="auto"` corresponds to a fixed cutoff of
   `0.5` on the anomaly score. Because solitude's inliers sit at
   ~0.55-0.60, that cutoff would flag everything. We chose
   `outlier_threshold = 0.6` empirically as the value that:

   * leaves clean small mock cohorts (~8 samples, random noise) with
     **0 flagged samples**, mirroring sklearn's `contamination="auto"`
     on that data, and
   * still isolates the genuinely anomalous sample in a
     "1 outlier in 30" stress test (the injected outlier's score
     consistently lands above 0.6; inliers stay below).

   The docstring for `by_isolation_forest()` records this rationale.

### What still doesn't match, even after all of the above

Outlier *lists* will still differ on:

* **Borderline samples.** A sample whose true anomaly score sits near
  the cutoff in one library may be on the other side of it in the
  other. This is the largest residual divergence; nothing portable
  fixes it.
* **Very small cohorts (n < 10).** Score variance is high relative to
  the spread; expect 1-2 samples of disagreement at random.

### Practical guidance

* For **typical use** (just clean a cohort before downstream
  clustering), call `remove_outlier_samples(df)` with no arguments and
  trust the default. Setting `quiet = FALSE` (the default) prints the
  flagged `Sample_ID`s so you can sanity-check what got removed, and
  `result$anomaly_plot` is a `plotly` bar chart of the per-sample
  scores you can print or save to inspect the result visually.
* For **cross-language reporting**, never quote outlier *scores* as
  matching Python; quote outlier *sample IDs* and accept the boundary
  drift described above.
* For **tighter or looser detection**, adjust `outlier_threshold`
  (lower → more flags, higher → fewer) **or** pass `contamination`
  explicitly to switch to a percentile rule
  (`contamination = 0.05` flags the top 5% by score regardless of
  threshold; this also mirrors sklearn's `contamination = 0.05` API).

```r
# Default — recommended for most cohorts:
res <- remove_outlier_samples(df)
res$df             # filtered data (samples removed)
res$outlier_list   # which Sample_IDs were flagged
res$anomaly_plot   # auto-prints in RStudio Viewer / notebook output

# Tighter, for very conservative removal:
res <- remove_outlier_samples(df, outlier_threshold = 0.65)

# Percentile-based, mirrors sklearn's contamination = 0.05:
res <- remove_outlier_samples(df, contamination = 0.05)
```

## 5. Random seeds and train/test splits

Even with `set.seed(42)` in R and `random_state = 42` in Python, the
Mersenne-Twister state mapping differs, and the splitting algorithms
themselves differ:

* Python `sklearn.model_selection.train_test_split` shuffles indices
  with a permutation and partitions.
* R `sample()` uses a different permutation under the same seed.

The result: patient assignments to train vs test differ across
languages. For reproducibility:

* **Fix the split explicitly** by passing `train_individuals` and
  `test_individuals` to `train_with_normalise`.
* Or compare metrics averaged over multiple seeds.

## 6. Things that ARE identical

These can be `all.equal`'d across R and Python:

* `percentile_cutoff(x, p) == numpy.percentile(x, p)` for any non-NaN `x`.
* `get_distances(df, "euclidean" | "manhattan" | "minkowski")`'s output
  equals `sklearn.metrics.pairwise_distances` on the same numeric
  matrix.
* `get_evaluation_metrics` produces identical TP/FP/FN/TN/precision/
  recall/F1 for the same pairwise classification.
* `plate_correct_residuals_by_protein`'s residuals equal
  `statsmodels.OLS(...).resid` to floating-point precision.
* Connected-component membership from `igraph::components()` equals
  `networkx.connected_components()`.

## Decision matrix

* **"I want the closest match to my previous Python results"**
  → `classifier_backend = "randomForest"`, `method = "MDS"`,
  `mds_backend = "smacof"`. Fix the train/test split.
* **"I'm starting fresh; speed matters"**
  → defaults (`"ranger"`, `"UMAP"`, `"cmdscale"`).
* **"My cohort is heavily imbalanced (rare patient class)"**
  → `classifier_backend = "themis_smote"`.

## References

* Liu et al. (2008) — Isolation Forest.
* McInnes et al. (2018) — UMAP.
* Wright & Ziegler (2017) — `ranger`.
* `uwot` README — explicit notes on differences from `umap-learn`.
