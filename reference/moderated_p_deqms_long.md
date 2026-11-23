# DEqMS-style moderation for long contrast table

Splits by contrast group and applies
[`moderated_p_deqms`](https://fgcz.github.io/prolfqua/reference/moderated_p_deqms.md)
to each group independently.

## Usage

``` r
moderated_p_deqms_long(
  mm,
  count_col,
  group_by_col = "contrast",
  estimate = "diff",
  loess_span = 0.75
)
```

## Arguments

- mm:

  contrast result data.frame (long format)

- count_col:

  name of column with peptide/PSM count per protein

- group_by_col:

  column to group by (default "contrast")

- estimate:

  name of the fold change column

- loess_span:

  span parameter for LOESS fit (default 0.75)

## Value

data.frame with moderated columns added

## Note

Deprecated in prolfqua 1.9.5; use
`build_contrast_analysis(method = "deqms")`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
