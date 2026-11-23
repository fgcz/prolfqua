# Compute protein-level fold changes and p-values using the beta distribution

takes p-value of the scaled p-value

## Usage

``` r
summary_ropeca_median_p_scaled(
  contrasts_data,
  contrast = "contrast",
  subject_id = "protein_Id",
  estimate = "diff",
  statistic = "statistic",
  p.value = "moderated.p.value",
  max.n = 10
)
```

## Arguments

- contrasts_data:

  data frame

- contrast:

  name of column with contrast identifier

- subject_id:

  name of column with typically protein Id

- estimate:

  name of column with effect size estimate

- statistic:

  statistic name of column with statistic (typically t-statistics)

- p.value:

  name of column with moderated.p.value

- max.n:

  used to limit the number of peptides in probability computation.

## Value

data.frame with columns

## Note

Deprecated in prolfqua 1.9.5; use
`build_contrast_analysis(method = "ropeca_nested")`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
