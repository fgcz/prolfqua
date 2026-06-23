# TODO: `rfit_impute`

## Goal

Add an `rfit_impute` facade to `prolfqua`: the rank-regression analogue of
`lm_impute`.

The intent is to rescue proteins where plain `rfit` cannot estimate one or
more requested contrasts because the per-protein fit failed, became singular,
or was fitted on an incomplete design. The rescue should follow the same
principle as `lm_impute`: refit the affected protein after LOD imputation, but
do not trust the artificially small variance from the imputed fit. Instead,
borrow uncertainty from successful non-imputed donor fits.

This belongs in upstream `prolfqua`, not in `prolfquapp`.

## Current reference points

- Plain `rfit` is implemented through `StrategyRfit`, `strategy_rfit()`, and
  `ContrastsRfitFacade`.
- `StrategyRfit` already augments `Rfit::rfit()` objects with a model frame,
  terms, coefficient names, and the `rfit_prolfqua` subclass so the standard
  `Contrasts` path can use `linfct_from_model()` and `compute_contrast()`.
- `lm_impute` is implemented through `build_model_impute()`,
  `impute_refit_singular()`, and the `lm_imputed` S3 wrapper, which overrides
  `vcov()`, `sigma()`, and `df.residual()`.

The `rfit_impute` implementation should preserve this same separation:
`Contrasts`, `ContrastsModerated`, `Model`, and `ContrastsInterface` should
remain unaware of imputation.

## Important `Rfit` caveat: ties

`Rfit` computes residual ranks inside its dispersion calculation with
`rank(e, ties.method = "first")`. LOD imputation can create many identical
values, especially for proteins with whole missing groups. Therefore
`rfit_impute` can create tied residuals whose rank order depends on row order.

Implementation requirements:

- Keep sample row ordering deterministic when completing data with the sample
  template.
- Document that LOD-created ties follow `Rfit`'s first-order tie handling.
- Add a regression test with LOD-created ties that confirms the facade returns
  stable, non-crashing contrast results.

This is not a reason to avoid `rfit_impute`, but it is a reason to keep the
first implementation conservative and explicit.

## Proposed implementation

### 1. Add an imputed `rfit` S3 wrapper

Add a constructor analogous to `new_lm_imputed()`:

```r
new_rfit_imputed <- function(model, borrowed_vcov, borrowed_sigma, borrowed_df, n_observed) {
  attr(model, "borrowed_vcov") <- borrowed_vcov
  attr(model, "borrowed_sigma") <- borrowed_sigma
  attr(model, "borrowed_df") <- borrowed_df
  attr(model, "n_observed") <- n_observed
  class(model) <- c("rfit_imputed", class(model))
  model
}
```

Add S3 methods:

- `vcov.rfit_imputed()` returns the borrowed covariance matrix.
- `sigma.rfit_imputed()` returns the borrowed rank-scale estimate.
- `df.residual.rfit_imputed()` returns the chosen imputed degrees of freedom.

The coefficients must stay those from the LOD-imputed `Rfit::rfit()` refit.
The borrowed values affect uncertainty only.

### 2. Reuse the LOD completion/refit workflow

Reuse the useful parts of `build_model_impute()` / `impute_refit_singular()`:

- first fit all proteins normally;
- identify failed, singular, or incomplete-coefficient fits;
- build a sample template from annotation columns;
- complete each affected protein to all expected sample rows;
- impute missing response values with LOD;
- clamp observed values with `pmax(response, lod)`;
- refit using `strategy_rfit()`;
- replace only successfully rescued model rows.

Keep `build_model()` unchanged.

### 3. Borrow covariance conservatively

Do not copy the `lm_impute` `"sigma"` method blindly.

For `lm`, `borrow_method = "sigma"` can use:

```r
borrowed_sigma^2 * summary(new_model)$cov.unscaled
```

`rfit` does not expose the same `lm`-style `cov.unscaled` object. The first
`rfit_impute` implementation should therefore support only a full covariance
borrowing mode:

- collect successful, non-singular `rfit_prolfqua` donor fits;
- extract named `vcov()` matrices;
- require matching coefficient names and dimensions;
- compute the element-wise median covariance matrix;
- borrow the median `sigma()` value, which is `tauhat` for `rfit_prolfqua`;
- use either observed or borrowed residual df, mirroring `lm_impute`.

If donor covariance matrices cannot be aligned by coefficient names, fail the
imputed rescue for those rows rather than silently falling back to an
`lm`-specific approximation.

### 4. Add the facade

Add `ContrastsRfitImputeFacade` in the same style as
`ContrastsLMImputeFacade`:

```r
strat <- strategy_rfit(full_formula, ...)
self$model <- build_model_rfit_impute(
  lfqdata,
  strat,
  lod = lod,
  df_method = df_method
)
self$contrast <- ContrastsModerated$new(Contrasts$new(self$model, contrasts))
self$config <- self$contrast$get_config()
```

Register it as:

```r
rfit_impute = .builtin_facade_entry("ContrastsRfitImputeFacade", "same")
```

The output schema should match plain `rfit` and `lm_impute`, including the
standard columns:

- `modelName`
- `contrast`
- `avgAbd`
- `diff`
- `FDR`
- `statistic`
- `std.error`
- `df`
- `p.value`
- `conf.low`
- `conf.high`
- `sigma`

Rescued rows should be tagged through the existing imputed model-name path,
ending up as `WaldTest_moderated_imputed` after moderation.

## Constraints

- Do not add observation weights. `Rfit::rfit()` has no `weights` argument, so
  `nr_children` weighting is not supported for `rfit_impute`, matching plain
  `rfit`.
- Do not add broad public API beyond the new method/facade unless the shared
  code extraction clearly requires it.
- Keep backend-specific details in the builder/wrapper layer. Downstream
  report and contrast code should see the normal contrast schema.
- Treat covariance borrowing failures as failed rescues, not as a reason to
  produce unlabelled or underestimated uncertainty.

## Suggested tests

Add tests near the existing rfit and imputation coverage:

1. `new_rfit_imputed()` preserves `rfit_prolfqua` behavior while overriding
   `vcov()`, `sigma()`, and `df.residual()`.
2. The imputed builder rescues at least some failed/singular/incomplete `rfit`
   rows in a fixture with missing groups.
3. `ContrastsRfitImputeFacade$get_contrasts()` returns the standard contrast
   schema and includes the `facade = "rfit_impute"` column.
4. Rescued rows are tagged as `WaldTest_moderated_imputed`.
5. `get_missing()` for `rfit_impute` is less than or equal to plain `rfit` on
   the same fixture.
6. A fixture with LOD-created ties returns stable, non-crashing results.
7. The existing `rfit` null-model fallback coefficient-name fix still works
   for imputed and non-imputed fits.
8. `to_wide()` and `get_Plotter()` work on the facade.

## Open decision

Whether `df_method = "observed"` should remain the default, as in `lm_impute`,
or whether `rfit_impute` should default to borrowed df because rank-regression
fixtures can have very small observed df after missingness. Start with the
`lm_impute` default for consistency unless validation shows unstable p-values.
