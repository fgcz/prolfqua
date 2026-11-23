# DEqMS-style count-dependent variance moderation

Applies count-dependent empirical Bayes variance shrinkage to a contrast
result table. Proteins quantified from many peptides get less shrinkage;
proteins from few peptides get more.

## Usage

``` r
moderated_p_deqms(
  mm,
  count_col,
  df = "df",
  estimate = "diff",
  loess_span = 0.75,
  confint = 0.95
)
```

## Arguments

- mm:

  data.frame from one contrast group with columns: sigma, df, statistic,
  std.error, and the estimate column

- count_col:

  name of column with peptide/PSM count per protein

- df:

  name of the degrees of freedom column

- estimate:

  name of the fold change column

- loess_span:

  span parameter for LOESS fit (default 0.75)

- confint:

  confidence level for intervals (default 0.95)

## Value

data.frame with added moderated.\* columns

## Note

Deprecated in prolfqua 1.9.5; use
`build_contrast_analysis(method = "deqms")`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
