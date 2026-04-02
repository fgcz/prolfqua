# compute group averages

used in p2621, p2109

## Usage

``` r
contrasts_linfct(
  models,
  linfct,
  subject_Id = "protein_Id",
  contrastfun = prolfqua::compute_lmer_contrast
)
```

## Examples

``` r
modelSummary_A <- sim_build_models_lm()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
m <- get_complete_model_fit(modelSummary_A$modelDF)

factor_contrasts <- linfct_factors_contrasts( m$linear_model[[1]])

factor_levelContrasts <- contrasts_linfct( m,
        factor_contrasts,
        subject_Id = "protein_Id",
        contrastfun = prolfqua::compute_contrast)
#> contrasts_linfct
```
