# Fix prolfqua limpa Vignette Build Failure

## Problem

Windows bundle builds fail while installing `prolfqua` with `build_vignettes = TRUE`. Two upstream `prolfqua`/`limpa`
compatibility issues have appeared during vignette builds:

1. `AggregateLimpa` called `limpa::dpc(expr_matrix, b1.upper = 1)` unconditionally, but the installed limpa version on
   Windows did not accept `b1.upper`.
2. After fixing that call, `vignettes/limpa_example.Rmd` still fails in chunk `lowlevel_dpc_hyper` because it assumes
   `dpc_est$mu.prior`, `dpc_est$df.prior`, and `dpc_est$s2.prior` are numeric values and calls `round()` directly. On
   the Windows build host at least one of these fields is non-numeric or unavailable.

## Root Cause

`prolfqua` assumes one specific `limpa::dpc()` return shape. The vignette should present the hyperparameter diagnostics
only when those fields are numeric, and otherwise report that the field is not available in the installed limpa return
object. This is an upstream `prolfqua` documentation/build issue, not a `pd_metaboanalyst` bundle issue.

## Fix

- Keep the committed `AggregateLimpa` compatibility fix: call `limpa::dpc()` with `b1.upper` only when that formal
  argument exists.
- Fix `vignettes/limpa_example.Rmd` so `lowlevel_dpc_hyper` formats DPC hyperparameters through a small helper that only
  rounds numeric scalar values.
- Keep `pd_metaboanalyst` bundle behavior unchanged: it should still build `prolfqua` and `prolfquapp` vignettes.

## Validation

- Render `vignettes/limpa_example.Rmd` locally.
- Run `git diff --check`.
- Re-run the Windows bundle build after committing the upstream fix.
