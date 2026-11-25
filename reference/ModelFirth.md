# Firth logistic regression models of detection

Holds one Firth logistic regression fit (logistf) per subject, modelling
whether a protein is detected in a sample, with the model table in
`models$model_df`. Contrasts are computed from it with
[ContrastsFirth](https://fgcz.github.io/prolfqua/reference/ContrastsFirth.md).

## Value

An R6 class generator.

## See also

Other modelling:
[`AnovaExtractor`](https://fgcz.github.io/prolfqua/reference/AnovaExtractor.md),
[`Contrasts`](https://fgcz.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsComparison`](https://fgcz.github.io/prolfqua/reference/ContrastsComparison.md),
[`ContrastsDEqMSFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsDEqMSVoomFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsDEqMSVoomFacade.md),
[`ContrastsFacadeBase`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.md),
[`ContrastsFirth`](https://fgcz.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsFirthFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsHurdleFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsHurdleFacade.md),
[`ContrastsHurdleNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsHurdleNestedFacade.md),
[`ContrastsLMFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsLMImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLMImputeFacade.md),
[`ContrastsLimma`](https://fgcz.github.io/prolfqua/reference/ContrastsLimma.md),
[`ContrastsLimmaFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimmaFacade.md),
[`ContrastsLimmaImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimmaImputeFacade.md),
[`ContrastsLimmaVoomFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimmaVoomFacade.md),
[`ContrastsLimmaVoomImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimmaVoomImputeFacade.md),
[`ContrastsLimpaFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimpaFacade.md),
[`ContrastsLimpaNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimpaNestedFacade.md),
[`ContrastsLmerNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLmerNestedFacade.md),
[`ContrastsModerated`](https://fgcz.github.io/prolfqua/reference/ContrastsModerated.md),
[`ContrastsModeratedDEqMS`](https://fgcz.github.io/prolfqua/reference/ContrastsModeratedDEqMS.md),
[`ContrastsPlotter`](https://fgcz.github.io/prolfqua/reference/ContrastsPlotter.md),
[`ContrastsRLMFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsRLMFacade.md),
[`ContrastsROPECA`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECA.md),
[`ContrastsROPECANestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECANestedFacade.md),
[`ContrastsROPECAV2`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAV2.md),
[`ContrastsROPECAV2NestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAV2NestedFacade.md),
[`ContrastsRfitFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsRfitFacade.md),
[`ContrastsRfitImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsRfitImputeFacade.md),
[`ContrastsTable`](https://fgcz.github.io/prolfqua/reference/ContrastsTable.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`Model`](https://fgcz.github.io/prolfqua/reference/Model.md),
[`ModelLimma`](https://fgcz.github.io/prolfqua/reference/ModelLimma.md),
[`StrategyLM`](https://fgcz.github.io/prolfqua/reference/StrategyLM.md),
[`StrategyLimma`](https://fgcz.github.io/prolfqua/reference/StrategyLimma.md),
[`StrategyLimpa`](https://fgcz.github.io/prolfqua/reference/StrategyLimpa.md),
[`StrategyLmer`](https://fgcz.github.io/prolfqua/reference/StrategyLmer.md),
[`StrategyLogistf`](https://fgcz.github.io/prolfqua/reference/StrategyLogistf.md),
[`StrategyRLM`](https://fgcz.github.io/prolfqua/reference/StrategyRLM.md),
[`StrategyRfit`](https://fgcz.github.io/prolfqua/reference/StrategyRfit.md),
[`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md),
[`build_model()`](https://fgcz.github.io/prolfqua/reference/build_model.md),
[`build_model_glm_protein()`](https://fgcz.github.io/prolfqua/reference/build_model_glm_protein.md),
[`build_model_impute()`](https://fgcz.github.io/prolfqua/reference/build_model_impute.md),
[`build_model_limma()`](https://fgcz.github.io/prolfqua/reference/build_model_limma.md),
[`build_model_limma_impute()`](https://fgcz.github.io/prolfqua/reference/build_model_limma_impute.md),
[`build_model_limma_voom()`](https://fgcz.github.io/prolfqua/reference/build_model_limma_voom.md),
[`build_model_limma_voom_impute()`](https://fgcz.github.io/prolfqua/reference/build_model_limma_voom_impute.md),
[`build_model_limpa()`](https://fgcz.github.io/prolfqua/reference/build_model_limpa.md),
[`compute_borrowed_variance()`](https://fgcz.github.io/prolfqua/reference/compute_borrowed_variance.md),
[`compute_borrowed_variance_limma()`](https://fgcz.github.io/prolfqua/reference/compute_borrowed_variance_limma.md),
[`df.residual.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/df.residual.rfit_prolfqua.md),
[`group_label()`](https://fgcz.github.io/prolfqua/reference/group_label.md),
[`impute_from_model()`](https://fgcz.github.io/prolfqua/reference/impute_from_model.md),
[`impute_refit_singular()`](https://fgcz.github.io/prolfqua/reference/impute_refit_singular.md),
[`list_facades()`](https://fgcz.github.io/prolfqua/reference/list_facades.md),
[`lookup_facade()`](https://fgcz.github.io/prolfqua/reference/lookup_facade.md),
[`lr_test()`](https://fgcz.github.io/prolfqua/reference/lr_test.md),
[`merge_contrasts_results()`](https://fgcz.github.io/prolfqua/reference/merge_contrasts_results.md),
[`new_imputed_model()`](https://fgcz.github.io/prolfqua/reference/new_imputed_model.md),
[`nobs.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/nobs.rfit_prolfqua.md),
[`pivot_model_contrasts_to_wide()`](https://fgcz.github.io/prolfqua/reference/pivot_model_contrasts_to_wide.md),
[`prolfqua_analysis_methods`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md),
[`register_facade()`](https://fgcz.github.io/prolfqua/reference/register_facade.md),
[`sigma.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/sigma.rfit_prolfqua.md),
[`strategy_limma()`](https://fgcz.github.io/prolfqua/reference/strategy_limma.md),
[`strategy_limpa()`](https://fgcz.github.io/prolfqua/reference/strategy_limpa.md),
[`strategy_logistf()`](https://fgcz.github.io/prolfqua/reference/strategy.md),
[`vcov.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/vcov.rfit_prolfqua.md)

## Super class

[`ModelInterface`](https://fgcz.github.io/prolfqua/reference/ModelInterface.md)
-\> `ModelFirth`

## Public fields

- `models`:

  list with `model_df`, one fitted model per subject, and the `strategy`

- `model_name`:

  name of model

- `subject_id`:

  e.g. protein_Id

- `p.adjust`:

  function to adjust p-values

## Methods

### Public methods

- [`ModelFirth$new()`](#method-ModelFirth-initialize)

- [`ModelFirth$get_coefficients()`](#method-ModelFirth-get_coefficients)

- [`ModelFirth$get_anova()`](#method-ModelFirth-get_anova)

- [`ModelFirth$anova_histogram()`](#method-ModelFirth-anova_histogram)

- [`ModelFirth$clone()`](#method-ModelFirth-clone)

Inherited methods

- [`ModelInterface$coef_histogram()`](https://fgcz.github.io/prolfqua/reference/ModelInterface.html#method-coef_histogram)
- [`ModelInterface$coef_pairs()`](https://fgcz.github.io/prolfqua/reference/ModelInterface.html#method-coef_pairs)
- [`ModelInterface$coef_volcano()`](https://fgcz.github.io/prolfqua/reference/ModelInterface.html#method-coef_volcano)
- [`ModelInterface$get_provenance_log()`](https://fgcz.github.io/prolfqua/reference/ModelInterface.html#method-get_provenance_log)
- [`ModelInterface$provenance()`](https://fgcz.github.io/prolfqua/reference/ModelInterface.html#method-provenance)

------------------------------------------------------------------------

### `ModelFirth$new()`

initialize

#### Usage

    ModelFirth$new(
      models,
      model_name = "modelFirth",
      subject_id = "protein_Id",
      p.adjust = prolfqua::adjust_p_values
    )

#### Arguments

- `models`:

  list with a `model_df` data.frame (one fitted model per subject) and
  the `strategy`

- `model_name`:

  name of model

- `subject_id`:

  subject column name

- `p.adjust`:

  method to adjust p-values

------------------------------------------------------------------------

### `ModelFirth$get_coefficients()`

return model coefficient table

#### Usage

    ModelFirth$get_coefficients()

------------------------------------------------------------------------

### `ModelFirth$get_anova()`

return anova table

#### Usage

    ModelFirth$get_anova()

------------------------------------------------------------------------

### `ModelFirth$anova_histogram()`

histogram of ANOVA results

#### Usage

    ModelFirth$anova_histogram(what = "p.value")

#### Arguments

- `what`:

  one case-sensitive choice: "p.value" (default) or "FDR".

------------------------------------------------------------------------

### `ModelFirth$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ModelFirth$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 10, with_missing = TRUE,
  weight_missing = 0.5, seed = 3)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 5 unconfigured columns: average_prot_abundance, sd, mean, N, Replicate.
#> setup_analysis: added 23 rows with NA response for feature and sample combinations absent from the input.
mod <- build_model_glm_protein(LFQData$new(istar$data, istar$config), "~ group_")
tmp <- mod$get_coefficients()
stopifnot(nrow(tmp) == 30)
mod$coef_histogram()
#> $plot

#> 
#> $name
#> [1] "Coef_Histogram_modelFirth.pdf"
#> 
mod$coef_pairs()
#> $plot
#> # A tibble: 10 × 4
#>    subject_id  `(Intercept)`   group_B group_Ctrl
#>    <chr>               <dbl>     <dbl>      <dbl>
#>  1 7IZdVV~0841         2.20   9.70e-16  -2.20e+ 0
#>  2 AZPG26~2091         2.20   1.26e-15  -4.27e-16
#>  3 AoNKbb~6908         2.20  -1.69e-16  -1.35e+ 0
#>  4 EnbuYk~0350         0.847 -8.47e- 1   4.14e-16
#>  5 JnEDsz~4938        -0.847 -1.35e+ 0  -4.18e-17
#>  6 NNASGK~1417         2.20   1.26e-15  -4.27e-16
#>  7 XYrp6h~5793         2.20   1.26e-15  -4.27e-16
#>  8 quTD7H~4566         2.20  -1.35e+ 0  -2.20e+ 0
#>  9 tCZCHm~6695         2.20  -1.35e+ 0  -2.16e-15
#> 10 tHE075~8231         0.847  1.35e+ 0   1.59e-16
#> 
#> $name
#> [1] "Coef_Pairsplot_modelFirth.pdf"
#> 
mod$get_anova()
#> Warning: `get_anova()` is not implemented for Firth models; returning `NULL`.
#> NULL
mod$coef_volcano()
#> $plot

#> 
#> $name
#> [1] "Coef_volcano_plot_modelFirth.pdf"
#> 
mod$anova_histogram()
#> Warning: `anova_histogram()` is not implemented for Firth models; returning `NULL`.
#> NULL
```
