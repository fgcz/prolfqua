# Changelog

## prolfqua 1.6.3

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
