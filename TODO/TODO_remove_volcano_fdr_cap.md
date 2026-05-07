# Remove Default Volcano FDR Cap

## Problem

Volcano plots in reports are capped at `-log10(FDR) = 4` because
`ContrastsPlotter$volcano()` and `$volcano_plotly()` default to
`min_score = 0.0001`. Users expect smaller non-zero FDR values to remain visible
instead of being flattened at the top of the plot.

## Plan

- [x] Change the default volcano significance floor from `1e-4` to no floor for
  positive finite values.
- [x] Keep exact zero and negative significance values finite by replacing only
  those with the smallest positive observed value in the same score column.
- [x] Preserve the explicit `min_score` argument for callers that still want a
  cap/floor.
- [x] Apply the same behavior to the exported `volcano_plotly()` helper.
- [x] Add regression tests and regenerate documentation.
