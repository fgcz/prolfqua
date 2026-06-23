# Fix: rfit crash on Rfit null-model fallback (unnamed coefficients)

## Status: ARCHIVED 2026-06-23 — source fix landed and covered

The `prolfqua` source work is complete: the fix is on `main`, the regression
test is present, and `Rfit` is declared in `DESCRIPTION` `Suggests`. The
`prolfquapp` image build has also been updated to install `Rfit`.

No local staged WU347631/32824 input was available to prove the exact original
workunit end-to-end from this checkout. Treat any future real-WU rerun as an
operations/deployment validation task, not an active `prolfqua` source TODO.

Fix commit: `b58cc714`; merged to `main` in `c085acf9`.
Found while following up on **bfabric WU347631** (FGCZ order/container 32824),
an A414 DEA job run with `model = NONE` → `model_extra = rfit`.

Follow-up on 2026-06-22:

- The fix and regression test are present on `prolfqua/main`.
- Local source validation passes when the working tree is loaded:

  ```bash
  Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-ContrastsRfit.R")'
  ```

  Result: 47 passed, 0 failed, 0 skipped. One expected warning remains from the rank-deficiency test
  (`lsfit(): 'X' matrix was collinear`).
- Running `testthat::test_file()` without `devtools::load_all()` still exercises the installed prolfqua and failed here
  because the installed package is stale. Install/rebuild prolfqua before any end-to-end DEA validation or container test.
- No staged WU347631/32824 input directory was found under `/Users/wolski/projects/prolfqua_fml`; only the
  `slurmworker/config/A414_DEA` configuration is present. The real-job validation therefore remains open.

## Symptom

The DEA aborted in the contrast step with:

```
model formula: normalized_abundance ~ G_
determine linear functions:
ERROR Caused by error in `tibble::as_tibble()`:
! Columns 1 and 2 must be named.  ✖ Empty names found at locations 1 and 2.
```

`prolfqua_dea.sh` exited 0 despite the R error (wrapper masks the failure),
so bfabric saw an empty/failed run with no clear cause.

## Root cause

Inside `Rfit::rfit.default()`:

```r
bhat  <- lsfit(x, yhat, intercept = FALSE)$coef   # NAMED by colnames(x)
bhat0 <- c(alphahat0, rep(0, length(bhat) - 1))    # UNNAMED, slopes = 0
if (disp(bhat, x, y, scores) > disp(bhat0, x, y, scores)) {
    bhat <- bhat0          # null-model fallback -> coefficients() is UNNAMED
}
res <- list(coefficients = bhat, ...)
```

When a protein's groups **do not separate** (the rank fit is no better than
the intercept-only model), `rfit()` returns `c(median(y), 0, ...)` — an
**unnamed** coefficient vector with the slope(s) hard-zeroed, **at full QR
rank**.

The existing guard in `StrategyRfit$model_fun` only catches QR rank
deficiency (`fit$qrx1$rank < ncol(fit$x)`). The null fallback keeps full rank,
so it slipped through. Downstream:

`linfct_from_model()` → `.model_coeff_matrix()` does
`colnames(coeff_matrix) <- names(coeffs)` → `NULL` → `.coeff_weights_factor_levels()`
→ `tibble::as_tibble()` throws "Columns must be named".

**One degenerate protein aborts the entire DEA for all proteins.**

## Why it was never caught

`tests/testthat/test-ContrastsRfit.R` is gated on
`skip_if_not_installed("Rfit")`, and Rfit is not in the CI / dev R library,
so the rfit backend has effectively never been exercised in CI.

## Fix applied

`R/rfit.R`, in `StrategyRfit$model_fun`, immediately after the `rfit()` call
and before the rank-deficiency guard:

```r
if (is.null(names(fit$coefficients))) {
  names(fit$coefficients) <- colnames(fit$x)
}
```

`fit$x` (the model matrix) always carries the correct column names
(`(Intercept)`, `G_T`, …) and is length-matched to the coefficients, so this
restores the contract `linfct_from_model()` / `vcov.rfit_prolfqua()` rely on.
The affected protein is then reported with `diff ≈ 0` / non-significant — the
honest result for a group showing no rank-based effect — instead of crashing
or being silently dropped.

Regression test added: `rfit restores names when it falls back to the null
model` — fits a non-separating two-group design, asserts the names are
restored, the slope is 0, and `linfct_from_model()` succeeds.

## Reproduction (for re-validation)

A minimal trigger (needs Rfit installed):

```r
d <- data.frame(y = c(13,15,17, 14,15,16), G_ = factor(rep(c("T","CT"), each = 3)))
fit <- Rfit::rfit(y ~ G_, data = d)
stopifnot(is.null(names(fit$coefficients)))   # unnamed, coef = c(15, 0)
```

End-to-end against the real job: the staged input from WU347631 was run
through `prolfqua/prolfquapp:2.0.22` with
`prolfqua_dea.sh -y config.yaml --software prolfquapp.DIANN -m rfit ...` and
reproduced the crash exactly. The data has 4 tissue groups (SM, T, AT, CT);
most proteins are observed in only 2–3 of them, and the offending protein was
observed in T & CT with non-separating abundances.

## Remaining work (handoff)

1. **Install/rebuild before validating.** The local installed prolfqua can be
   stale even when the source tree contains the fix. Use the package Makefile
   or a rebuilt prolfquapp image before re-running WU347631.
2. **Finish the patched end-to-end run.** The mechanism is confirmed in
   isolation and source-level tests pass, but the full patched DEA against
   WU347631 still needs to be run with the staged input. Re-run and confirm it
   produces a populated result directory (xlsx / contrasts) with the degenerate
   protein present and diff ≈ 0.
3. **Ship it through the image.** The job runs the bundled
   `prolfquapp:2.0.22` Docker image. The prolfqua commit does not change the
   deployed container — rebuild prolfqua into a new prolfquapp image (or test
   via a dev image mounting the patched package) before re-running WU347631.
4. **Keep rfit tests active in CI.** `Rfit` is declared in `DESCRIPTION`
   `Suggests`, and it is installed locally. Confirm the GitHub Actions
   dependency setup installs `Suggests` so `test-ContrastsRfit.R` actually
   runs there; otherwise regressions here stay invisible. Consider also
   exercising the *partial-model* path
   (`.linfct_partial_model`) with real per-protein group dropout.
5. **Optional, separate:** `prolfqua_dea.sh` returns exit 0 on R failure —
   consider propagating the non-zero status so bfabric/app-runner sees the
   real failure instead of an empty result dir. (Lives in prolfquapp /
   slurmworker, not prolfqua.)
6. **Consider** whether other backends share this assumption that
   `coef()` is always named (the lm path uses NA, not unnamed zeros, so it is
   fine; limma/limpa paths not audited here).
