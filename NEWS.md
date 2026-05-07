# prolfqua 1.6.1

- Added a `check-bioc` Makefile target and a Bioconductor Docker check image.
- Updated BiocCheck-related documentation and vignette metadata.
- `ContrastsPlotter$volcano()` and `volcano_plotly()` no longer cap positive FDR/p-value scores at `1e-4` by default.
- Renamed `LFQData$to_wide()` to `data_wide()` and removed deprecated compatibility wrappers.
- Split `nr_obs_experiment()` into `nr_children_experiment()` and `nr_features_experiment()`.
- Made `LFQData` data/config internals private and moved callers to accessor methods.
- Refactored transformer and aggregation helpers around current `LFQData` APIs.
