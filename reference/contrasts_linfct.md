# Compute group averages or contrasts for each model

used in p2621, p2109

## Usage

``` r
contrasts_linfct(
  models,
  linfct,
  subject_id = "protein_Id",
  contrastfun = prolfqua::compute_lmer_contrast
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

## Examples

``` r
modelSummary_A <- sim_build_models_lm()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, group, mean, N, Replicate.
#> setup_analysis: added 59 rows with NA response for feature and sample combinations absent from the input.
#> fit (Model): 59 missing responses dropped (e.g. 0EfVhX~0087, 7cbcrd~5725, 9VUkAq~4703, BEJI92~5282, CGzoYe~2147, ...): NA response (missing observation) left out of the fit
#> See $provenance().
m <- get_complete_model_fit(modelSummary_A$model_df)

linfct <- linfct_from_model(m$linear_model[[1]], as_list = FALSE)
factor_contrasts <- linfct_matrix_contrasts(linfct, c(A_vs_B = "TreatmentA - TreatmentB"))
#> Warning: The `x` argument of `as_tibble.matrix()` must have unique column names if
#> `.name_repair` is omitted as of tibble 2.0.0.
#> ℹ Using compatibility `.name_repair`.
#> ℹ The deprecated feature was likely used in the prolfqua package.
#>   Please report the issue at <https://github.com/fgcz/prolfqua/issues>.

factor_levelContrasts <- contrasts_linfct( m,
        factor_contrasts,
        subject_id = "protein_Id",
        contrastfun = prolfqua::compute_contrast)
```
