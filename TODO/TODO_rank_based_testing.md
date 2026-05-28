# TODO: Rank-Based Testing Backends

## Question

Can prolfqua support a factorial rank-based test with an API similar to the
current `lm` backend, including formula handling, contrasts, and facade
integration?

## Findings

There is no true factorial version of the Wilcoxon signed-rank test that shares
the `lm` modelling API. The signed-rank test is a paired two-condition test, not
a coefficient-based model with arbitrary contrasts, residual degrees of freedom,
and interaction terms.

Rank-based methods therefore fall into two separate backend types:

1. `Rfit::rfit()` as an LM-like rank-regression backend that can plausibly reuse
   the classic `build_model()` / `compute_contrast()` path.
2. A dedicated Wilcoxon contrast backend for paired two-condition comparisons,
   implemented as a `Contrasts*` adapter driven by a `ContrastConfiguration`.

## Candidate Backends

### `Rfit::rfit()` — primary target

Closest candidate for an `lm`-like rank-based backend.

- Formula interface similar to `lm`; supports k-way factorial designs.
- Exposes `coef.rfit`, `vcov.rfit`, and `summary.rfit` (Estimate / Std.Error /
  t-value / p-value), plus a `tauhat` scale estimate.
- `drop.test()` gives reduction-in-dispersion tests with an F-distribution and
  explicit `df1` / `df2`, usable for ANOVA-style omnibus output.

Resolved checks (previously "open" — Rfit satisfies these):

- `coef()` returns a named coefficient vector matching the design matrix.
- `vcov()` is available and gives coefficient-wise standard errors.
- A usable scale (`tauhat`) and residual df (`n - p`) exist for Wald-style
  contrast output.
- Contrast statistics can be returned in the standard schema (`diff`,
  `std.error`, `statistic`, `df`, `p.value`, `FDR`).

**The decisive open check — contrast scaffold construction.** The classic path
does *not* build contrasts from `coef()` / `vcov()` alone. `linfct_from_model()`
→ `.model_coeff_matrix()` (`R/tidyMS_contrasts.R`) reconstructs factor-level and
interaction structure by reaching into the fitted **model frame** and **terms**:

```r
data <- m$model                                   # needs the model frame
coeffs <- coef(m)
interaction_columns <- intersect(attributes(terms(m))$term.labels, colnames(data))
```

An `rfit` object is built around `$coefficients`, `$residuals`,
`$fitted.values`, and the `tauhat` scale; it does **not** appear to carry a
`$model` frame, and `terms()` support is uncertain. (Verify on an installed
Rfit — it is not currently in the library: `names(rfit(...))`, `terms(fit)`,
`fit$model`.) If either is absent, `linfct_from_model()` will fail or misbuild
the scaffold even when `coef()` / `vcov()` are correct.

Mitigation (cheap, keeps the classic path): have `strategy_rfit()`'s `model_fun`
attach the model frame and a terms attribute to the returned fit, e.g. store
`model.frame(formula, data)` as `$model` and `terms(formula)` on the object, so
`m$model`, `coef(m)`, `vcov(m)`, and `terms(m)` all resolve. If that proves
fragile, write a dedicated `linfct` path that derives the scaffold from the
design matrix / coefficient names instead (this is the skill's Step 5 rule).

### Wilcoxon signed-rank — dedicated paired backend

Useful for paired two-condition comparisons, not a general factorial model.

- Do not route through `StrategyLM`.
- Implement as a dedicated `ContrastsWilcoxon` adapter.
- Report the Hodges-Lehmann estimate (or paired median difference) as the
  effect.
- Handle pairing, missing pairs, ties, and zero differences explicitly.
- `coin::wilcoxsign_test()` is a suitable engine.

**Use `ContrastConfiguration` rather than a bespoke schema.** This backend's
native columns deliberately diverge from the LM schema (no `std.error` / `df` /
`sigma`; an extra `n_pairs`). This is exactly the SAINTexpress-style case the
`ContrastConfiguration` mechanism handles (`R/ContrastsInterface.R`,
`R/ContrastConfiguration.R`). The adapter should:

1. keep its native columns in `get_contrasts()`, and
2. set a Wilcoxon-flavoured `ContrastConfiguration` in `initialize()`, mapping
   `contrast_col`, `effect_col` (HL estimate), `score_col` (the test statistic),
   `pvalue_col`, `fdr_col`, and `avg_abundance_col` to the native column names.

Then `filter_significant`, `get_ora`, `contrast_summary_table`, and `get_rank`
are inherited — no overrides needed. Because Wilcoxon has a p-value, `get_rank`
falls back correctly to `sign(effect) * -log10(p.value)`.

Native output columns:

- `modelName`, `contrast`, `avgAbd`, `diff` (HL estimate), `FDR`, `statistic`,
  `p.value`, `n_pairs`

### `ARTool::art()` — possible later backend

Factorial rank-based ANOVA via aligned rank transform.

- Supports factorial designs and interactions.
- Better suited to omnibus ANOVA than LM-style contrasts on the original
  response — ART contrasts are not ordinary linear functions of the response, so
  contrast handling needs careful adapter design.
- Not the first integration target.

### `coin` — nonparametric test framework

`coin::wilcoxsign_test()` and `coin::independence_test()`.

- Strong for exact/asymptotic rank tests and blocked designs.
- Not coefficient-based; better as the engine behind a dedicated test backend
  (see Wilcoxon above) than as input to the LM-style contrast machinery.

## Recommended Direction

Start with `Rfit::rfit()` if the goal is an LM-like rank-regression backend:

- Add `strategy_rfit()` and reuse `build_model()` / `compute_contrast()`.
- The gating risk is **not** `vcov()` reliability but model-frame / `terms()`
  availability for `linfct_from_model()`. Attach them in `model_fun`, or write a
  dedicated `linfct` path.
- Add tests for two-factor designs, interactions, rank deficiency, missingness,
  and multiple contrasts.

Implement Wilcoxon separately for paired nonparametric testing:

- Add a dedicated builder and `ContrastsWilcoxon` adapter.
- Do not expose it as a full factorial model.
- Drive divergent columns through `ContrastConfiguration` and inherit the
  default interface methods.
- Make pairing requirements explicit in facade validation.

## Interface Notes

Any new backend must still satisfy the public contracts:

- `ModelInterface`, if it is model-like.
- `ContrastsInterface`, always — inherit the defaults and populate
  `ContrastConfiguration` in `initialize()` rather than hard-coding columns.
- A facade registered in `R/ContrastsFacades.R` (built-in via
  `.builtin_facade_entry`) if it should be reachable via
  `build_contrast_analysis()`. Facades classify as `needs = "same"` or
  `needs = "nested"`.

Keep backend-specific differences inside the adapter, and return the standard
contrast schema where possible; where columns genuinely diverge, resolve them
generically through `ContrastConfiguration` accessors (`cfg$effect_col`, etc.).
Avoid adding public API beyond the requested backend.
