# Fix prolfqua limpa Vignette Build Failure

## Problem

Windows bundle builds fail while installing `prolfqua` with `build_vignettes = TRUE` because `prolfqua` calls
`limpa::dpc(expr_matrix, b1.upper = 1)` from `AggregateLimpa`. The installed limpa version on the Windows machine does
not accept `b1.upper`, so vignettes that exercise the limpa facade fail with:

```text
unused argument (b1.upper = 1)
```

## Root Cause

`prolfqua` assumes a newer `limpa::dpc()` API and unconditionally passes an optional argument that is not available in
all supported/installed limpa versions.

## Fix

- Fix `prolfqua` upstream by calling `limpa::dpc()` with `b1.upper` only when the installed limpa version supports that
  formal argument.
- Add a focused regression test for both old-style and new-style `dpc()` function signatures.
- Keep `pd_metaboanalyst` bundle behavior unchanged: it should still build `prolfqua` and `prolfquapp` vignettes.
