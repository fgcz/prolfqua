# Compute contrasts of a mixed-effects model with lmerTest

Compute contrasts of a mixed-effects model with lmerTest

## Usage

``` r
compute_lmer_contrast(model, linfct, ddf = c("Satterthwaite", "Kenward-Roger"))
```

## Arguments

- model:

  mixed effects model

- linfct:

  linear function

- ddf:

  method to determine denominator degrees of freedom

## Value

a tibble with one row per row of `linfct` and columns `lhs`, `estimate`,
`std.error`, `df`, `statistic`, `conf.low`, `conf.high`, `p.value` and
`sigma`; `NA` (with a warning) if `linfct` has no rows or the model is
rank-deficient

## Note

Deprecated in prolfqua 1.9.5; use
[`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md).
See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
