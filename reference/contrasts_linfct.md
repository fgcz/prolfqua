# Compute group averages or contrasts for each model

used in p2621, p2109

## Usage

``` r
contrasts_linfct(
  models,
  linfct,
  subject_id = "protein_Id",
  contrastfun = .compute_lmer_contrast
)
```

## Arguments

- models:

  model table with one fit per subject in the column `linear_model`,
  e.g. the `model_df` of
  [`model_analyse()`](https://fgcz.github.io/prolfqua/reference/model_analyse.md)

- linfct:

  linear functions: a matrix with one row per contrast and one column
  per model coefficient, or a list with one such matrix per model

- subject_id:

  columns identifying a subject

- contrastfun:

  function computing the contrasts of one fit, e.g.
  [`compute_lmer_contrast()`](https://fgcz.github.io/prolfqua/reference/compute_lmer_contrast.md)

## Value

a tibble with one row per subject and row of `linfct`, holding the
`subject_id` columns, the columns returned by `contrastfun`, and
`isSingular`, `sigma.model` and `df.residual.model`; subjects whose
contrasts failed are dropped with a warning

## Note

Deprecated in prolfqua 1.9.5; use
[`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md).
See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
