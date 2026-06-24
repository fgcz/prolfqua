# prolfqua 1.6.3

- `get_contrast()` now derives `group_1`/`group_2` from the contrast's left/right side expressions, fixing mislabeled per-group columns for averaging contrasts such as `(group_A + group_B)/2 - group_Ctrl` (previously it used the first two group tokens). Simple `A - B` contrasts are unaffected. A contrast that is not a difference `LHS - RHS` now errors with a clear message instead of silently extracting the first two group tokens. Nested contrasts that reference an earlier contrast by name (e.g. `Interaction = "AvsB_gv_X - AvsB_gv_Z"`) remain supported in the `ContrastsMissing` / `lm_missing` path.
- `setup_analysis()` now stops with an informative error (listing the offending keys) when a hierarchy-key/sample combination has more than one observation, instead of silently returning a different-schema count table that crashed downstream. Pass `debug = TRUE` to recover the old behaviour and return the count table for inspection.
- Removed the unused `impute_with_zcomp()`, `estimate_lod_global()`, and `function_lod_quantile()` exports (and the `zCompositions` dependency). For missing-value imputation use `AggregateLimpa$new(lfqdata, impute_only = TRUE)$aggregate()`.
- Hardened `plot_pca()`: errors early on duplicated sample names, an all-missing matrix, or too few samples instead of returning `NULL` (which broke `pca_plotly()`); joins scores to annotation with an explicit `by`; makes `prcomp(center = TRUE, scale. = FALSE)` explicit.

# prolfqua 1.6.1

- Added a `check-bioc` Makefile target and a Bioconductor Docker check image.
- Updated BiocCheck-related documentation and vignette metadata.
- `ContrastsPlotter$volcano()` and `volcano_plotly()` no longer cap positive FDR/p-value scores at `1e-4` by default.
- Renamed `LFQData$to_wide()` to `data_wide()` and removed deprecated compatibility wrappers.
- Split `nr_obs_experiment()` into `nr_children_experiment()` and `nr_features_experiment()`.
- Made `LFQData` data/config internals private and moved callers to accessor methods.
- Refactored transformer and aggregation helpers around current `LFQData` APIs.
