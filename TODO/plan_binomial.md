# Plan: native `binomial_nested` model

## Objective

Add a `binomial_nested` facade for differential peptide/precursor detection without adding `msqrob2` or another
dependency. Reuse the existing Firth missingness preparation and prolfqua's generic model, contrast, and moderation
machinery wherever possible.

The result is a protein-level log odds ratio for detection. Combining it with an intensity model in a hurdle facade is
out of scope for this change.

## Method choices and limitations

The method intentionally follows the count component discussed for the msqrob2 hurdle approach, implemented with base
R and prolfqua's existing limma dependency:

- Nested child-feature rows are collapsed to detected/undetected counts per parent and sample.
- Child features are treated as exchangeable binomial trials conditional on the sample design; child identity and
  within-parent correlation are not modelled.
- In the name `binomial_nested`, `nested` describes the required input shape and child-to-parent aggregation, not a
  hierarchical or mixed-effects model. `firth_nested` remains the alternative that models binary child rows and retains
  child identity as a covariate.
- `prior_count = 0.1` is retained as a symmetric stabilizing pseudocount. It is not equivalent to Firth's Jeffreys
  penalty and shrinks low-child-count proteins more strongly. A denominator-scaled prior is deliberately not used
  because it would define a different method. Separation tests will pin the resulting finite estimate and document the
  role of this parameter.
- Quasibinomial dispersions are noisy at small residual degrees of freedom and especially weak for single-child
  proteins. Moderating them with `limma::squeezeVar()` is a pragmatic empirical-Bayes approximation rather than an
  exact quasibinomial likelihood result. The default lower bound of one prevents sub-binomial overconfidence.

## Reuse the Firth data preparation

The Firth path already performs the required preparation:

1. `LFQData$complete_cases()` creates missing child-feature/sample rows.
2. `encode_bin_resp()` converts the response to detected (`1`) or missing (`0`).
3. The completed rows retain parent, child, sample, isotope, and experimental-factor columns.

Generalize `.prepare_logistf_lfqdata()` to a backend-neutral internal helper such as
`.prepare_detection_lfqdata()`. Both Firth and `binomial_nested` will use it. Update the existing direct internal test
reference when the helper is renamed.

Do not create a second pipeline that reconstructs child identities or sample combinations independently.

## Count summary

The only new data preparation is a summary of the completed binary rows by parent, sample, isotope label, and design
factors:

```r
detected <- sum(bin_resp)
available <- dplyr::n()
undetected <- available - detected
```

Because `complete_cases()` establishes one row per configured hierarchy/sample combination, the denominator comes
directly from the existing representation.

The facade will create this parent-by-sample data frame before fitting and pass it directly to `build_model()` with
`subject_id = lfqdata$subject_id()`. A new `LFQData` wrapper is unnecessary because `build_model()` already accepts a
data frame. Consequently, residual degrees of freedom are based on samples, not child rows.

## Binomial strategy

Add `StrategyBinomial` and `strategy_binomial()` using the standard strategy contract. For each parent, fit:

```r
stats::glm(
  cbind(detected + prior_count, undetected + prior_count) ~ design,
  family = stats::quasibinomial(),
  data = x
)
```

`quasibinomial()` provides the same coefficient fit as the binomial model while exposing a Student-t coefficient table,
covariance matrix, residual degrees of freedom, and dispersion in forms understood by the existing `Model` and
`Contrasts` classes.

The strategy will provide:

- `model_fun`
- `isSingular`
- `contrast_fun = compute_contrast`
- `model_name = "binomial_nested"`
- `report_columns` consistent with the other strategies
- `is_mixed = FALSE`
- `anova_df = get_anova_df(test = "F")`
- `df_residual = stats::df.residual`
- `sigma = stats::sigma`

The scale must match the Pearson dispersion embedded in `vcov(glm)`. R's `stats::sigma.glm()` method returns
`sqrt(summary(model)$dispersion)`, so ordinary `stats::sigma()` dispatch is correct and reusable here. A regression test
will pin this S3-dispatch contract and its covariance-scale coherence.

Public facade parameters:

- `prior_count = 0.1`
- `binomial_bound = TRUE`

No new `ModelBinomialNested` or `ContrastsBinomialNested` adapter is planned. The existing `Model`, `Contrasts`,
`ContrastsModerated`, `ContrastsPlotter`, and `pivot_model_contrasts_to_wide()` implementations should provide the
public contracts.

## Moderation-floor generalization

Thread an optional posterior-variance floor through the existing shared path:

1. `moderated_p_limma()` accepts the optional floor.
2. `moderated_p_limma_long()` forwards it.
3. `ContrastsModerated` stores it and passes it from `get_contrasts()`.

The default is no floor, leaving all existing facades unchanged. For `binomial_nested`:

- `binomial_bound = TRUE` supplies a floor of one.
- `binomial_bound = FALSE` supplies no floor.

Apply the floor after the existing infinite-prior fallback. For rows below the floor, set posterior variance to the
floor and total posterior degrees of freedom to `Inf` before calculating statistics, p-values, and confidence bounds.
Explicitly test and document that the standard `df` output can therefore be infinite.

## Facade and contrast configuration

Add only `ContrastsBinomialNestedFacade`, following `ContrastsLmerNestedFacade` as the structural template:

1. Validate nested `LFQData`.
2. Call the shared detection preparation helper.
3. Summarize completed binary rows into the parent-by-sample count data frame.
4. Fit that data frame with `strategy_binomial()` and `build_model(..., subject_id = parent)`.
5. Construct `Contrasts$new()` and wrap it in `ContrastsModerated$new()` with the requested variance floor.
6. Set `facade_name = "binomial_nested"`.
7. Register the facade with `needs = "nested"`.

Set an explicit `ContrastConfiguration` using the standard result columns and `supports_dea_qc = FALSE`. The current
configuration API maps column roles but has no unit-label fields, so log-odds semantics belong in the facade and result
documentation rather than a new configuration API.

The standard output will use:

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
- sample/factor alignment and parent `subject_id` selection;
- coefficients and Pearson dispersion matching an independent quasibinomial fit;
- the invariant `strategy$sigma(fit) * sqrt(diag(summary(fit)$cov.unscaled)) == sqrt(diag(vcov(fit)))`;
- moderated statistics matching a direct calculation using `sqrt(summary(fit)$dispersion)`;
- moderation and the optional floor matching `limma::squeezeVar()` plus the documented bound;
- `prior_count` and `binomial_bound` behaviour;
- complete separation producing a finite, pinned pseudocount-dependent estimate;
- a single-child protein and proteins with insufficient residual degrees of freedom;
- log odds contrasts and standard errors matching direct matrix calculations;
- multiple contrasts, multifactor designs, interactions, absent levels, and rank deficiency;
- reused `ModelInterface` and `ContrastsInterface` methods, including plots and wide output;
- explicit configuration roles and `supports_dea_qc = FALSE`;
- facade dispatch and nested-input validation;
- absence of `msqrob2` from source and package metadata.

## Validation

1. Run `make document` and inspect generated exports and Rd changes.
2. Run focused tests and lintr for changed source and tests.
3. Run `make test`, `make check-fast`, and `make install` in prolfqua.
4. Run `make test` and `make build` in prolfquapp.
5. Confirm unrelated worktree changes remain intact.

## Completion criteria

- [ ] Firth and binomial backends share missingness preparation.
- [ ] No new package dependency is introduced.
- [ ] The strategy uses the Pearson dispersion scale embedded in `vcov()`.
- [ ] Only a binomial strategy, thin facade, count summary, and small moderation generalization are added.
- [ ] Existing generic model and contrast adapters provide the public contracts.
- [ ] Log odds contrasts and moderated inference agree with independent reference calculations.
- [ ] Method assumptions and limitations are documented.
- [ ] Core and downstream validation pass.
