# Changelog

## prolfqua 1.6.1

- Added a `check-bioc` Makefile target and a Bioconductor Docker check
  image.
- Updated BiocCheck-related documentation and vignette metadata.
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
