# Implementation Plan: `Rfit::rfit()` Rank-Regression Backend

## Status: DONE

Landed on `Modelling2R6` in commits `277325df` (formatting prep) and
`19f224c3` (rfit backend), with the `@examples` follow-up + the rank-
deficient test added on top.

Outcomes vs. plan:

- **Step 0 spike** — resolved decisively. `Rfit::rfit()` exposes
  `coef`/`vcov`/`tauhat` with names identical to `lm`, but the fit
  carries no `$model`, no `terms()`, and `vcov.rfit` returns an unnamed
  matrix.
- **Step 2a primary path used** — augmenting the fit in `model_fun` with
  `model.frame()`, `terms()`, and the `rfit_prolfqua` subclass plus three
  S3 methods (`vcov`/`df.residual`/`sigma`) restores the contract that
  `linfct_from_model()` / `compute_contrast()` rely on. **Step 2b lm
  fallback was not needed** — no core changes required.
- **Rank-deficient designs**: `rfit` silently zeros unestimable
  coefficients (unlike `lm`'s NA), and `vcov()` then Cholesky-fails. The
  adapter detects this in `model_fun` via `qrx1$rank < ncol(x)` and
  routes the affected proteins through `.error_handler`, so they surface
  via `get_missing()` rather than crashing the contrast loop.
- **Open questions resolved**:
  - `tauhat`-as-`sigma` adopted (mirrors `rlm`); documented as a
    caveat in the facade.
  - Double-fit cost is irrelevant — fallback path unused.
- **Tests**: 42 passing in `test-ContrastsRfit.R` (S3 glue, linfct
  parity with `lm`, sign/correlation agreement, 2-factor + interactions,
  rank deficiency, `build_contrast_analysis` dispatch, registry).
  Existing facade suite: 263 passing, no regression.
- **Docs**: all five Rd files (`StrategyRfit`, `strategy_rfit`,
  `ContrastsRfitFacade`, and the three S3 methods) carry unguarded
  runnable `@examples`; `Rfit` added to `Suggests`.

Archived 2026-05-28.

## Goal

Add an `lm`-like rank-based regression backend using `Rfit::rfit()`, reusing the
classic `build_model()` / `compute_contrast()` path and the Wald contrast
machinery. User-facing entry: `method = "rfit"` via `build_contrast_analysis()`.

## Strategy: clone the `StrategyRLM` shape

`StrategyRLM` (`R/tidyMS_R6_Modelling.R:257`) is the template — `MASS::rlm` is
already a non-`lm` backend that reuses `compute_contrast`, `stats::df.residual`,
`stats::sigma`, and works with `linfct_from_model()`. `StrategyRfit` mirrors it.

The only real difference from rlm: an `rfit` object does **not** carry a
`$model` frame or `terms`, which `linfct_from_model()` → `.model_coeff_matrix()`
(`R/tidyMS_contrasts.R:3`) requires (`m$model`, `coef(m)`, `terms(m)`).

## Step 0 — Verification spike (do first, ~30 min, throwaway)

Install Rfit, then on a simulated 2-factor protein confirm:

```r
f  <- Rfit::rfit(Sepal.Length ~ Species, data = iris)
names(f); class(f)
coef(f); names(coef(f))          # must match lm coef names
vcov(f)                          # coefficient-wise covariance, named?
tryCatch(terms(f),  error = \(e) "no terms")
f$model                          # likely NULL
m_lm <- lm(Sepal.Length ~ Species, data = iris)
identical(names(coef(f)), names(coef(m_lm)))   # alignment check
```

Decision gate:
- **If** attaching `$model` + `$terms` to the fit makes `linfct_from_model(f)`
  return the same scaffold as the lm fit → **primary path** (Step 2a).
- **Else** → **lm-fallback path** (Step 2b): fit a parallel per-protein `lm`
  purely to obtain the linfct.

## Step 1 — `StrategyRfit` + `strategy_rfit()`

New file `R/rfit.R` (keeps Rfit isolated; `Rfit` in `Suggests`, guarded with
`requireNamespace`). Mirror `StrategyRLM` fields: `formula`, `model_name`,
`report_columns`, `is_mixed = FALSE`, `anova_df = get_anova_df(test = "F")`.

Methods:
- `isSingular(model)` — NA coef or `df_residual < 2` (copy rlm logic).
- `contrast_fun(...)` — `compute_contrast(...)`.
- `df_residual(model)` — `length(model$residuals) - length(coef(model))`.
- `sigma(model)` — rank-regression scale `model$tauhat` (confirm field name in
  spike). Used for downstream eBayes moderation; if unsuitable, document that
  rfit results should not be limma-moderated.

`strategy_rfit()` wrapper following `strategy_lm()` (`:236`).

## Step 2a — `model_fun` (primary: augment the fit)

```r
model_fun = function(x, pb, get_formula = FALSE) {
  if (get_formula) return(self$formula)
  if (!missing(pb)) pb$tick()
  tryCatch({
    fit <- Rfit::rfit(self$formula, data = x)
    fit$model <- stats::model.frame(self$formula, data = x)  # for .model_coeff_matrix
    fit$terms <- stats::terms(self$formula, data = x)        # for terms(fit)
    fit
  }, error = .error_handler)
}
```

Failures must return a character string via `.error_handler` (so
`model_analyse()`'s `!is.character(x)` success test holds). Verify `terms(fit)`
and `coef(fit)`/`vcov(fit)` then satisfy `linfct_from_model()` unchanged.

## Step 2b — lm fallback (only if 2a fails the gate)

Keep `model_fun` returning the plain rfit fit, but build the linfct from a
parallel `lm(self$formula, data = x)` and feed it to `compute_contrast(rfit_fit,
linfct = <from lm>)`. This needs a small hook where the per-model linfct is
computed; locate it in `build_model()` (`R/tidyMS_build_model.R:133`) /
`model_analyse()` (`:625`) and confirm the linfct source can be overridden per
strategy without disturbing the lm/rlm paths. Add an assertion that
`names(coef(rfit_fit))` matches the lm linfct columns; define behavior when an
rfit design is rank-deficient differently from lm.

## Step 3 — `ContrastConfiguration`

rfit emits the standard LM schema (`diff`, `std.error`, `statistic`, `df`,
`p.value`, `FDR`), so the default LM-flavoured `ContrastConfiguration` applies —
no custom config needed. Confirm the facade exposes it via
`self$config <- self$contrast$get_config()`.

## Step 4 — Facade `ContrastsRfitFacade`

In `R/ContrastsFacades.R`, clone `ContrastsLMFacade` (`:483`):
- `.assert_aggregated_facade_input()`,
- prepend response to formula,
- `strat <- strategy_rfit(full_formula, ...)`,
- `self$model <- build_model(lfqdata, strat)`,
- build `Contrasts$new(self$model, contrasts)`,
- set `self$config <- self$contrast$get_config()`,
- add a `facade` column.

Decide weights: `rfit()` has no `weights` argument, so `nr_children` weighting
is **not** supported — document this divergence from `lm`/`limma`. `needs = "same"`.

## Step 5 — Register

Add `rfit = .builtin_facade_entry("ContrastsRfitFacade", "same")` to
`FACADE_REGISTRY` (`R/ContrastsFacades.R:1165`).

## Step 6 — Tests (`tests/testthat/test-ContrastsFacades.R` + new file)

- construction: `strategy_rfit()` returns `StrategyRfit`.
- `linfct_from_model()` on an rfit fit equals the lm-derived scaffold
  (the core regression test for the primary path).
- `get_contrasts()` returns the standard schema.
- two-factor design + interactions; multiple contrasts.
- rank deficiency (missing factor level) and missingness — defined behavior.
- `to_wide()`, `get_Plotter()` work.
- invariant: rfit fold-changes track lm signs; not required to match magnitude.

## Step 7 — Docs

Roxygen on `StrategyRfit`, `strategy_rfit`, `ContrastsRfitFacade` with runnable
`@examples` (no `\dontrun`). Add `Rfit` to `DESCRIPTION` `Suggests`. Run
`make document`. Note in CLAUDE.md facade tables that `rfit` is `needs = "same"`,
unweighted.

## Open questions / risks

- `tauhat` as `sigma` — is rank scale meaningful for eBayes moderation? If not,
  expose rfit only as an unmoderated Wald backend.
- rfit handling of rank-deficient designs vs lm — the main edge-case risk.
- Double-fit cost in the fallback path (lm + rfit per protein).
