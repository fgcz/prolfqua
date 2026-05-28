# Shorten heatmap sample labels

## Problem

The tabbed SE report correlation and abundance heatmaps become unreadable for
metabolomics datasets with many long sample names. The prefixes are often
shared and the suffix carries the useful replicate/group information.

## Plan

- [x] Add an internal suffix-label helper that keeps the last N characters.
- [x] Update `plot_heatmap_cor()` to pass shortened column labels to
  `pheatmap` while preserving the matrix and annotation row names.
- [x] Update abundance `plot_heatmap()` and `plot_raster()` to use suffix
  labels for sample columns while keeping middle-truncated row labels for long
  feature IDs.
- [x] Keep the labels configurable with a conservative default of 20 characters.
- [x] Expose the option through `LFQDataPlotter$heatmap()` and
  `LFQDataPlotter$raster()`.
- [x] Add regression tests that long sample names are displayed as suffixes in
  correlation, abundance heatmap, and raster plots.
- [x] Regenerate docs and run targeted tests.

## Verification

- `make document`
- `Rscript -e "devtools::test(filter = '^plotting_functions$')"`
- Re-ran after extending the same suffix-label behavior to abundance heatmap
  and raster plots.
