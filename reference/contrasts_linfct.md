# compute group averages

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

## Examples

``` r
modelSummary_A <- sim_build_models_lm()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
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
#> contrasts_linfct
```
