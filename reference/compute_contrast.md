# Compute contrasts for possibly incomplete models

only keeps non NA coefficients.

## Usage

``` r
compute_contrast(m, linfct, confint = 0.95)
```

## Arguments

- m:

  linear model generated using lm

- linfct:

  linear function

- confint:

  confidence interval default 0.95

## Value

a data.frame with one row per row of `linfct` and columns `lhs`,
`sigma`, `df`, `estimate`, `std.error`, `statistic`, `p.value`,
`conf.low` and `conf.high`; the estimates are `NA` for rows that use a
coefficient the model could not estimate

## Note

Deprecated in prolfqua 1.9.5; use
[`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md).
See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
