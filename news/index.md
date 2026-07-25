# Changelog

## prolfqua 1.7.0

- Added the `binomial_nested` facade for nested LFQ data. It reuses the
  Firth backend’s completed peptide-detection representation, fits
  protein-by-sample detected/undetected peptide counts with a
  quasibinomial GLM, and reports moderated protein-level log odds-ratio
  contrasts. A symmetric pseudo-count stabilizes separated fits, and
  posterior dispersion is bounded below by one by default
  (`binomial_bound = FALSE` disables the bound).
- **Variance moderation now uses
  [`limma::squeezeVar()`](https://rdrr.io/pkg/limma/man/squeezeVar.html)
  directly.**
  [`moderated_p_limma()`](https://wolski.github.io/prolfqua/reference/moderated_p_limma.md)
  — and therefore every moderated Wald-test facade (`lm`, `rlm`, `rfit`,
  `lm_impute`, `lm_missing`, `rfit_impute`, `lmer_nested`) and the DEqMS
  small-group fallback — now calls
  [`limma::squeezeVar()`](https://rdrr.io/pkg/limma/man/squeezeVar.html)
  instead of the vendored `squeezeVarRob()` fork (limma is already a
  dependency). Modern limma estimates the empirical-Bayes prior robustly
  on the fractional / low residual degrees of freedom these estimators
  produce, without the `min_df` workaround prolfqua inherited from
  MSqRob. This slightly changes moderated p-values, FDR and confidence
  bounds; on the IonStar spike-in benchmark AUC is unchanged, FDR
  calibration is maintained, and ≥99.7% of significance calls are
  identical. See the `SqueezeVar_comparison` vignette in
  `prolfquabenchmark`. (`ropeca_nested` is unaffected — it aggregates
  raw peptide-level p-values and never used variance moderation. A
  contrast with fewer than three features now returns unmoderated
  results, matching limma, where the empirical-Bayes prior is
  undefined.)
- **Fixed the scale used to moderate `rlm` contrasts.**
  `StrategyRLM$sigma()` now returns the robust scale
  [`MASS::rlm`](https://rdrr.io/pkg/MASS/man/rlm.html) builds its
  standard errors from (`fit$s`) instead of
  [`stats::sigma()`](https://rdrr.io/r/stats/sigma.html) (the
  ordinary-residual scale). Because the two differ,
  `ContrastsModerated`’s `sigma / sqrt(var.post)` rescaling was leaving
  a spurious per-protein factor in moderated `rlm` statistics, p-values,
  FDR and confidence bounds. Only the `rlm` facade’s moderated output is
  affected — `lm` and `rfit` were already coherent, and `rlm`’s
  unmoderated statistics (`diff`, `std.error`) were always correct.
- **Removed the exported `squeezeVarRob()`** and its internal helpers
  (`fitFDist_LG`, `fitFDistRobustly_LG`); the whole vendored
  `squeezeVarRob.R` file is gone. Use
  [`limma::squeezeVar()`](https://rdrr.io/pkg/limma/man/squeezeVar.html)
  instead. DEqMS moderation now calls
  [`limma::trigammaInverse()`](https://rdrr.io/pkg/limma/man/trigammainverse.html)
  (identical to the copy that was removed). This is a breaking change
  only for code that called `prolfqua::squeezeVarRob()` directly — the
  ecosystem does not. Also dropped the now-unused `statmod` from
  `Imports`.
- User-facing errors at prolfqua’s public entry points now carry typed
  condition classes (all inheriting from `prolfqua_error`), so callers
  and tests can catch failures by class instead of matching message
  text. Covered:
  [`setup_analysis()`](https://wolski.github.io/prolfqua/reference/setup_analysis.md),
  [`build_contrast_analysis()`](https://wolski.github.io/prolfqua/reference/build_contrast_analysis.md),
  the `LFQData$set_data()` / `set_config_value()` mutators,
  `LFQData$get_Aggregator()`, the aggregators’ “must aggregate a
  multi-level hierarchy” / “aggregate first” checks,
  `LFQDataTransformer$center_to_reference()`, and
  `ContrastsSimpleImpute$get_contrasts()`. `LFQData$set_data()` and
  `set_config_value()` additionally now reject data that no longer
  contains the columns required by the current configuration, instead of
  accepting invalid state silently.
- Abundance-density plots now carry per-sample Plotly highlight keys,
  allowing interactive reports to fade non-hovered sample curves while
  preserving the existing ggplot output.
- **Breaking — contrast schema.** The `modelName` column of
  `get_contrasts()` output is now the selected facade key (`lm`, `rlm`,
  `rfit`, `lm_impute`, `lm_missing`, `limma`, `limma_impute`,
  `limma_voom`, `limma_voom_impute`, `limpa`, `deqms`, `deqms_voom`,
  `firth`, `lmer_nested`, `ropeca_nested`, `firth_nested`,
  `limpa_nested`) instead of the testing-schema label
  (`WaldTest_moderated`, `*_DEqMS`, `*_imputed`, …). Rescue/imputation
  state moved to a new `estimate_type` column with values `observed`,
  `lod_imputed`, or `missing_fallback`. The redundant `facade` column
  was removed (its role is now played by `modelName`). The
  moderated-Wald-test wording belongs in methods text, not per-row
  labels. Downstream code that filtered on
  `modelName == "WaldTest_moderated"` (or the `_imputed`/`_DEqMS`
  variants) or read the `facade` column must migrate to the facade key
  and `estimate_type`.
- [`build_contrast_analysis()`](https://wolski.github.io/prolfqua/reference/build_contrast_analysis.md)
  and the exported `FACADE_REGISTRY` now derive method dispatch and the
  available-method list from a single seeded registry
  ([`list_facades()`](https://wolski.github.io/prolfqua/reference/list_facades.md)),
  removing three hand-maintained copies that could drift. New exported
  base class `ContrastsFacadeBase` holds the shared facade plumbing; the
  19 built-in facades are now thin subclasses.
- `ContrastsPlotter` colours volcano/MA/score plots by `estimate_type`
  when present (and the colour column was left at its default), keeping
  LOD-imputed / fallback rows visually distinct now that `modelName` is
  uniform per run.
- Correctness fixes: `AggregateTopN` now validates `func` via
  `match.arg` (an invalid value errors instead of silently meaning
  `mean`); `Model$get_anova()` drops degenerate rows by `is.na(factor)`
  rather than a never-matching `"NULL"` string;
  [`sim_lfq_data()`](https://wolski.github.io/prolfqua/reference/sim_lfq_data.md)
  now honours its `mean_prot` argument;
  `LFQDataTransformer`/`LFQDataSummariser` deep-clone their input like
  the other decorators (no more mutating the caller’s object);
  `ContrastsFirth$get_linfct()` returns the linfct-annotated model
  copies instead of mutating the shared model;
  [`contrasts_fisher_exact()`](https://wolski.github.io/prolfqua/reference/contrasts_fisher_exact.md)
  pre-allocates `nrow(x)` results; `ContrastsMissing` now actually
  validates its output schema.
- Fixed the vectorized contrast path
  (`options(prolfqua.vectorize = TRUE)`) to match the default loop path
  on rank-deficient (aliased) model fits: a contrast whose weights fall
  on non-estimable coefficients now returns `NA` instead of a spurious
  fold-change of 0 with a `NaN` p-value. The invalid-row guard now
  counts nonzero weights rather than summing signed weights, so
  contrasts with canceling `+1`/`-1` weights on missing coefficients are
  correctly flagged. The default (non-vectorized) path was never
  affected.
- The vectorized contrast path now assigns stable `contrast_1`,
  `contrast_2`, … labels to unnamed contrasts, matching the documented
  naming behavior instead of leaking an expression into the result row
  names.
- [`get_contrast()`](https://wolski.github.io/prolfqua/reference/get_contrast.md)
  now derives `group_1`/`group_2` from the contrast’s left/right side
  expressions, fixing mislabeled per-group columns for averaging
  contrasts such as `(group_A + group_B)/2 - group_Ctrl` (previously it
  used the first two group tokens). Simple `A - B` contrasts are
  unaffected. A contrast that is not a difference `LHS - RHS` now errors
  with a clear message instead of silently extracting the first two
  group tokens. Nested contrasts that reference an earlier contrast by
  name (e.g. `Interaction = "AvsB_gv_X - AvsB_gv_Z"`) remain supported
  in the `ContrastsMissing` / `lm_missing` path.
- [`setup_analysis()`](https://wolski.github.io/prolfqua/reference/setup_analysis.md)
  now stops with an informative error (listing the offending keys) when
  a hierarchy-key/sample combination has more than one observation,
  instead of silently returning a different-schema count table that
  crashed downstream. Pass `debug = TRUE` to recover the old behaviour
  and return the count table for inspection.
- Removed the unused `impute_with_zcomp()`, `estimate_lod_global()`, and
  `function_lod_quantile()` exports (and the `zCompositions`
  dependency). For missing-value imputation use
  `AggregateLimpa$new(lfqdata, impute_only = TRUE)$aggregate()`.
- Hardened
  [`plot_pca()`](https://wolski.github.io/prolfqua/reference/plot_pca.md):
  errors early on duplicated sample names, an all-missing matrix, or too
  few samples instead of returning `NULL` (which broke `pca_plotly()`);
  joins scores to annotation with an explicit `by`; makes
  `prcomp(center = TRUE, scale. = FALSE)` explicit.
- Hardened abundance heatmaps for sparse significant-feature subsets:
  when row or column distances are non-finite because of missing values,
  [`plot_heatmap()`](https://wolski.github.io/prolfqua/reference/plot_heatmap.md)
  now falls back to the input order instead of returning a
  `ComplexHeatmap` object that fails during drawing.
- `LFQDataPlotter$heatmap()` now shows only the `top_n` most variable
  features (default 1000), ranked by the prolfqua per-feature statistic
  (CV for untransformed data, sd for transformed, via `LFQDataStats`).
  Row clustering uses
  [`stats::hclust`](https://rdrr.io/r/stats/hclust.html), which errors
  above 65536 features, so peptide-list / entrapment searches with tens
  of thousands of degenerate protein groups no longer crash the QC
  heatmap. Pass `top_n = NULL` (or `Inf`) to keep every feature.
- `StrategyLogistf` now uses Wald confidence intervals instead of
  profiling every coefficient, preventing `firth_nested` analyses from
  stalling for days on proteins with hundreds or thousands of peptides.

## prolfqua 1.6.1

- Added a `check-bioc` Makefile target and a Bioconductor Docker check
  image.
- Updated BiocCheck-related documentation and vignette metadata.
- `ContrastsPlotter$volcano()` and
  [`volcano_plotly()`](https://wolski.github.io/prolfqua/reference/volcano_Plotly.md)
  no longer cap positive FDR/p-value scores at `1e-4` by default.
- Renamed `LFQData$to_wide()` to `data_wide()` and removed deprecated
  compatibility wrappers.
- Split `nr_obs_experiment()` into
  [`nr_children_experiment()`](https://wolski.github.io/prolfqua/reference/nr_children_experiment.md)
  and
  [`nr_features_experiment()`](https://wolski.github.io/prolfqua/reference/nr_features_experiment.md).
- Made `LFQData` data/config internals private and moved callers to
  accessor methods.
- Refactored transformer and aggregation helpers around current
  `LFQData` APIs.
