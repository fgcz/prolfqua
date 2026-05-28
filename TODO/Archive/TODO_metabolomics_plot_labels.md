# Improve metabolomics plot labels

## Problem

Metabolomics datasets can have many samples and long feature IDs. Current
plotting defaults can produce unreadable reports:

- density plots show one legend entry per sample unless each caller suppresses
  the legend explicitly;
- raster/heatmap row labels use full row names when enabled, so long feature IDs
  consume too much plot space.

## Plan

- [x] Add internal helpers for automatic sample legend decisions and middle
  truncation of display labels.
- [x] Update density plots to hide the sample legend automatically above a
  configurable sample count threshold while preserving explicit `legend = TRUE`
  and `legend = FALSE` behavior.
- [x] Update raster and heatmap plotting to pass shortened display labels to
  `pheatmap` when row names are shown, without changing matrix row names or the
  underlying data.
- [x] Expose the row label length through `LFQDataPlotter$raster()` and
  `LFQDataPlotter$heatmap()` with conservative defaults.
- [x] Add focused tests for automatic density legend suppression and heatmap
  label truncation.
- [x] Run targeted tests and whitespace checks.

## Verification

- `Rscript -e "devtools::test(filter = '^plotting_functions$')"`
- `make document`
- `Rscript -e "devtools::test()"`
