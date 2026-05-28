# Named Volcano Scores for Contrast Plotters

## Context

`prolfquapp` report templates access `contrast_plotter$volcano()$FDR`.
For SAINT, the statistical column is `BFDR`, and
`prolfquasaint::ContrastsSAINTexpress$get_Plotter()` correctly configures the
volcano score as `BFDR`. This makes the generic report path miss `$FDR`, or
forces downstream aliases that obscure the SAINT-specific axis label.

## Plan

- Extend `prolfqua::ContrastsPlotter` volcano specifications with an optional
  `name` field.
- Keep the plotted score column unchanged, so `score = "BFDR"` still produces a
  y-axis labelled `-log10(BFDR)`.
- Use `name` only as the returned list key, so
  `list(score = "BFDR", name = "FDR")` returns `volcano()$FDR`.
- Update/extend tests for the generic plotter behavior.
- Update `prolfquasaint::ContrastsSAINTexpress$get_Plotter()` to request
  `name = "FDR"` for its BFDR volcano plot.
- Update SAINT tests/docs so callers can rely on `$FDR` while retaining BFDR
  semantics in the plot.

## Validation

- Run targeted `prolfqua` tests for `ContrastsPlotter`.
- Run targeted `prolfquasaint` tests for `ContrastsSAINTexpress`.
- Re-document and reinstall affected packages if the behavior is used from the
  installed CLI/report path.
