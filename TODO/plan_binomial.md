# Plan: native `binomial_nested` model

## Objective

Add a `binomial_nested` facade for differential peptide/precursor detection without adding `msqrob2` or another
dependency. Reuse the existing Firth missingness preparation and prolfqua's generic model, contrast, and moderation
machinery wherever possible.

The result is a protein-level log odds ratio for detection. Combining it with an intensity model in a hurdle facade is
out of scope for this change.

## Reuse-first design

The Firth path already performs the important data preparation:

1. `LFQData$complete_cases()` creates the missing child-feature/sample rows.
2. `encode_bin_resp()` converts the response to detected (`1`) or missing (`0`).
3. The completed rows retain the parent, child, sample, isotope, and experimental-factor columns.

Generalize `.prepare_logistf_lfqdata()` to a backend-neutral internal helper such as
`.prepare_detection_lfqdata()`. Both Firth and `binomial_nested` will use this helper. Do not create a second pipeline
that reconstructs child identities or sample combinations independently.

The only new data operation is a summary of those prepared rows by parent and sample:

```r
detected <- sum(bin_resp)
available <- dplyr::n()
undetected <- available - detected
```

Because `complete_cases()` establishes one row per configured hierarchy/sample combination, `available` comes directly
from the existing completed representation. Grouping will retain the sample, isotope, and factor columns required by
the model formula.

## Statistical implementation

Add a small `StrategyBinomial`/`strategy_binomial()` implementation compatible with `build_model()` and
`model_analyse()`. For each parent it will fit:

```r
stats::glm(
  cbind(detected + prior_count, undetected + prior_count) ~ design,
  family = stats::quasibinomial(),
  data = x
)
```

Using `quasibinomial()` gives the same coefficient fit as the binomial model while exposing the residual dispersion,
Student-t coefficient table, covariance matrix, and residual degrees of freedom in the forms expected by prolfqua's
existing `Model` and `Contrasts` classes.

The strategy will provide only the normal strategy contract:

- `model_fun`
- `isSingular`
- `contrast_fun = compute_contrast`
- `model_name = "binomial_nested"`
- `anova_df = get_anova_df(test = "F")`
- `df_residual = stats::df.residual`
- `sigma = stats::sigma`

Public parameters:

- `prior_count = 0.1`
- `binomial_bound = TRUE`

No new `ModelBinomialNested` or `ContrastsBinomialNested` class is planned. The existing `Model`, `Contrasts`,
`ContrastsModerated`, `ContrastsPlotter`, and `pivot_model_contrasts_to_wide()` implementations should satisfy the
public interfaces.

## Minimal generalization for dispersion bounding

`ContrastsModerated` already applies `limma::squeezeVar()` across proteins. Generalize that existing path with an
optional posterior-variance floor, defaulting to no floor so all current models remain unchanged.

For `binomial_nested`:

- `binomial_bound = TRUE` passes a floor of `1`.
- `binomial_bound = FALSE` passes no floor.
- A posterior variance below the floor is set to the floor and its posterior degrees of freedom become infinite.

This keeps empirical-Bayes moderation in one shared implementation rather than reproducing it in a binomial-specific
adapter.

## Facade

Add only `ContrastsBinomialNestedFacade`:

1. Validate nested `LFQData`.
2. Call the shared detection preparation helper.
3. Summarize the completed binary rows into detected/undetected counts.
4. Fit with `strategy_binomial()` and `build_model()`.
5. Construct `Contrasts$new()` and wrap it in `ContrastsModerated$new()` with the requested variance floor.
6. Expose the wrapped contrast configuration.
7. Register `binomial_nested` with `needs = "nested"`.

The standard contrast output will use:

- `modelName = "binomial_nested"`
- `estimate_type = "observed"`
- `diff` as the log odds ratio
- `avgAbd` as the average linear predictor on the log-odds scale

Document the facade in `build_contrast_analysis()` and add a NEWS entry. No DESCRIPTION dependency change is needed.

## Tests

Add focused tests for:

- the shared Firth preparation remaining unchanged after generalization;
- count summaries matching the completed Firth binary rows exactly;
- implicit missing rows becoming zero detections;
- sample/factor alignment;
- coefficients and raw dispersion matching an independent `stats::glm(..., family = quasibinomial())` fit;
- moderation matching `limma::squeezeVar()`;
- `prior_count` and `binomial_bound` behaviour;
- log odds contrasts and standard errors matching direct matrix calculations;
- multiple contrasts, multifactor designs, interactions, absent levels, rank deficiency, and insufficient residual df;
- the reused `ModelInterface` and `ContrastsInterface` methods, including plots and wide output;
- facade dispatch and nested-input validation;
- absence of `msqrob2` from source and package metadata.

## Validation

1. Run `make document` and inspect the generated exports and Rd changes.
2. Run focused tests and lintr for the changed source and tests.
3. Run `make test`, `make check-fast`, and `make install` in prolfqua.
4. Run `make test` and `make build` in prolfquapp.
5. Confirm unrelated worktree changes remain intact.

## Completion criteria

- [ ] Firth and binomial backends share the missingness preparation.
- [ ] No new package dependency is introduced.
- [ ] Only a binomial strategy, thin facade, count summary, and small moderation generalization are added.
- [ ] Existing generic model and contrast adapters provide the public contracts.
- [ ] Log odds contrasts and moderated inference agree with independent reference calculations.
- [ ] Core and downstream validation pass.
