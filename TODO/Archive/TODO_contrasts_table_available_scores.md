# TODO: Make ContrastsTable plots use available score columns

## Problem

`ContrastsTable$get_Plotter()` always configures volcano and histogram plots for both `p.value` and `FDR`.
Table-backed contrast results can legitimately contain `FDR` without `p.value`, for example SAINTexpress results
exported through the SummarizedExperiment/Quarto report path. Calling `contrast_plotter$volcano()$FDR` still tries to
build every configured volcano score first, so the missing `p.value` column aborts the whole plot.

## Plan

- Update `ContrastsTable$get_Plotter()` so volcano and histogram specifications are built only from score columns present
  in `self$contrast_result`.
- Prefer `p.value` and `FDR` when available, preserving existing output for normal DEA tables.
- Keep `statistic` score plots only when the `statistic` column is present.
- Add a focused test that a table with `FDR` but no `p.value` can render the FDR volcano and does not configure a p-value
  volcano.
- Run targeted `ContrastsTable`/plotter tests, format touched files, then install `prolfqua` with vignettes and rerender
  the SAINT example.
