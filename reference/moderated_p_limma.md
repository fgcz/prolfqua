# Moderate p-values - limma approach

Moderate p-values - limma approach

## Usage

``` r
moderated_p_limma(
  contrast_df,
  df = "df",
  estimate = "diff",
  robust = FALSE,
  confint = 0.95,
  variance_floor = NULL
)
```

## Arguments

- contrast_df:

  contrast table with the columns `sigma`, the residual degrees of
  freedom and the estimate

- df:

  name of the residual degrees-of-freedom column

- estimate:

  name of the effect-estimate column

- robust:

  passed to
  [`limma::squeezeVar()`](https://rdrr.io/pkg/limma/man/squeezeVar.html):
  estimate the prior robustly

- confint:

  confidence level of the moderated intervals

- variance_floor:

  optional positive lower bound for posterior variances

## Value

`contrast_df`, ungrouped, with the added columns `moderated.var.post`,
`moderated.var.prior`, `moderated.df.prior`, `moderated.std.error`,
`moderated.statistic`, `moderated.df.total`, `moderated.p.value`,
`moderated.conf.low` and `moderated.conf.high`

## Note

Deprecated in prolfqua 1.9.5; use
[`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md).
See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
