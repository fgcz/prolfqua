# compute group averages

used in p2621, p2109

## Usage

``` r
contrasts_linfct(
  models,
  linfct,
  subject_Id = "protein_Id",
  contrastfun = prolfqua::my_contest
)
```

## Examples

``` r
modelSummary_A <- sim_build_models_lm()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
#> Joining with `by = join_by(protein_Id)`
m <- get_complete_model_fit(modelSummary_A$modelDF)

factor_contrasts <- linfct_factors_contrasts( m$linear_model[[1]])

factor_levelContrasts <- contrasts_linfct( m,
        factor_contrasts,
        subject_Id = "protein_Id",
        contrastfun = prolfqua::my_contrast_V2)
#> contrasts_linfct
```
