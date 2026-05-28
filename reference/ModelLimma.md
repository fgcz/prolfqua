# R6 class representing a limma modelling result

R6 class representing a limma modelling result

R6 class representing a limma modelling result

## Value

An R6 class generator.

## Details

Same API as
[`Model`](https://wolski.github.io/prolfqua/reference/Model.md):
`get_anova()`, `get_coefficients()`, `coef_histogram()`,
`coef_volcano()`, `anova_histogram()`.

## See also

Other modelling:
[`AnovaExtractor`](https://wolski.github.io/prolfqua/reference/AnovaExtractor.md),
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsDEqMSVoomFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSVoomFacade.md),
[`ContrastsFirth`](https://wolski.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsFirthFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsFirthNestedFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthNestedFacade.md),
[`ContrastsLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsLMImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMImputeFacade.md),
[`ContrastsLMMissingFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMMissingFacade.md),
[`ContrastsLimma`](https://wolski.github.io/prolfqua/reference/ContrastsLimma.md),
[`ContrastsLimmaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaFacade.md),
[`ContrastsLimmaImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaImputeFacade.md),
[`ContrastsLimmaVoomFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaVoomFacade.md),
[`ContrastsLimmaVoomImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaVoomImputeFacade.md),
[`ContrastsLimpaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimpaFacade.md),
[`ContrastsLimpaNestedFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimpaNestedFacade.md),
[`ContrastsLmerNestedFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLmerNestedFacade.md),
[`ContrastsMissing`](https://wolski.github.io/prolfqua/reference/ContrastsMissing.md),
[`ContrastsModerated`](https://wolski.github.io/prolfqua/reference/ContrastsModerated.md),
[`ContrastsModeratedDEqMS`](https://wolski.github.io/prolfqua/reference/ContrastsModeratedDEqMS.md),
[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.md),
[`ContrastsRLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsRLMFacade.md),
[`ContrastsROPECA`](https://wolski.github.io/prolfqua/reference/ContrastsROPECA.md),
[`ContrastsROPECANestedFacade`](https://wolski.github.io/prolfqua/reference/ContrastsROPECANestedFacade.md),
[`ContrastsRfitFacade`](https://wolski.github.io/prolfqua/reference/ContrastsRfitFacade.md),
[`ContrastsTable`](https://wolski.github.io/prolfqua/reference/ContrastsTable.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`LR_test()`](https://wolski.github.io/prolfqua/reference/LR_test.md),
[`Model`](https://wolski.github.io/prolfqua/reference/Model.md),
[`ModelFirth`](https://wolski.github.io/prolfqua/reference/ModelFirth.md),
[`StrategyLM`](https://wolski.github.io/prolfqua/reference/StrategyLM.md),
[`StrategyLimma`](https://wolski.github.io/prolfqua/reference/StrategyLimma.md),
[`StrategyLimpa`](https://wolski.github.io/prolfqua/reference/StrategyLimpa.md),
[`StrategyLmer`](https://wolski.github.io/prolfqua/reference/StrategyLmer.md),
[`StrategyLogistf`](https://wolski.github.io/prolfqua/reference/StrategyLogistf.md),
[`StrategyRLM`](https://wolski.github.io/prolfqua/reference/StrategyRLM.md),
[`StrategyRfit`](https://wolski.github.io/prolfqua/reference/StrategyRfit.md),
[`build_contrast_analysis()`](https://wolski.github.io/prolfqua/reference/build_contrast_analysis.md),
[`build_model()`](https://wolski.github.io/prolfqua/reference/build_model.md),
[`build_model_glm_peptide()`](https://wolski.github.io/prolfqua/reference/build_model_glm_peptide.md),
[`build_model_glm_protein()`](https://wolski.github.io/prolfqua/reference/build_model_glm_protein.md),
[`build_model_impute()`](https://wolski.github.io/prolfqua/reference/build_model_impute.md),
[`build_model_limma()`](https://wolski.github.io/prolfqua/reference/build_model_limma.md),
[`build_model_limma_impute()`](https://wolski.github.io/prolfqua/reference/build_model_limma_impute.md),
[`build_model_limma_voom()`](https://wolski.github.io/prolfqua/reference/build_model_limma_voom.md),
[`build_model_limma_voom_impute()`](https://wolski.github.io/prolfqua/reference/build_model_limma_voom_impute.md),
[`build_model_limpa()`](https://wolski.github.io/prolfqua/reference/build_model_limpa.md),
[`build_model_logistf()`](https://wolski.github.io/prolfqua/reference/build_model_logistf.md),
[`compute_borrowed_variance()`](https://wolski.github.io/prolfqua/reference/compute_borrowed_variance.md),
[`compute_borrowed_variance_limma()`](https://wolski.github.io/prolfqua/reference/compute_borrowed_variance_limma.md),
[`compute_contrast()`](https://wolski.github.io/prolfqua/reference/compute_contrast.md),
[`compute_lmer_contrast()`](https://wolski.github.io/prolfqua/reference/compute_lmer_contrast.md),
[`contrasts_fisher_exact()`](https://wolski.github.io/prolfqua/reference/contrasts_fisher_exact.md),
[`df.residual.rfit_prolfqua()`](https://wolski.github.io/prolfqua/reference/df.residual.rfit_prolfqua.md),
[`get_anova_df()`](https://wolski.github.io/prolfqua/reference/get_anova_df.md),
[`get_complete_model_fit()`](https://wolski.github.io/prolfqua/reference/get_complete_model_fit.md),
[`get_p_values_pbeta()`](https://wolski.github.io/prolfqua/reference/get_p_values_pbeta.md),
[`group_label()`](https://wolski.github.io/prolfqua/reference/group_label.md),
[`impute_refit_singular()`](https://wolski.github.io/prolfqua/reference/impute_refit_singular.md),
[`is_singular_lm()`](https://wolski.github.io/prolfqua/reference/is_singular_lm.md),
[`linfct_all_possible_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_all_possible_contrasts.md),
[`linfct_factors_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_factors_contrasts.md),
[`linfct_from_model()`](https://wolski.github.io/prolfqua/reference/linfct_from_model.md),
[`linfct_matrix_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_matrix_contrasts.md),
[`list_facades()`](https://wolski.github.io/prolfqua/reference/list_facades.md),
[`lookup_facade()`](https://wolski.github.io/prolfqua/reference/lookup_facade.md),
[`merge_contrasts_results()`](https://wolski.github.io/prolfqua/reference/merge_contrasts_results.md),
[`model_analyse()`](https://wolski.github.io/prolfqua/reference/model_analyse.md),
[`model_summary()`](https://wolski.github.io/prolfqua/reference/model_summary.md),
[`moderated_p_deqms()`](https://wolski.github.io/prolfqua/reference/moderated_p_deqms.md),
[`moderated_p_deqms_long()`](https://wolski.github.io/prolfqua/reference/moderated_p_deqms_long.md),
[`moderated_p_limma()`](https://wolski.github.io/prolfqua/reference/moderated_p_limma.md),
[`moderated_p_limma_long()`](https://wolski.github.io/prolfqua/reference/moderated_p_limma_long.md),
[`new_lm_imputed()`](https://wolski.github.io/prolfqua/reference/new_lm_imputed.md),
[`pivot_model_contrasts_to_wide()`](https://wolski.github.io/prolfqua/reference/pivot_model_contrasts_to_wide.md),
[`plot_lmer_peptide_predictions()`](https://wolski.github.io/prolfqua/reference/plot_lmer_peptide_predictions.md),
[`register_facade()`](https://wolski.github.io/prolfqua/reference/register_facade.md),
[`sigma.rfit_prolfqua()`](https://wolski.github.io/prolfqua/reference/sigma.rfit_prolfqua.md),
[`sim_build_models_lm()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lm.md),
[`sim_build_models_lmer()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lmer.md),
[`sim_build_models_logistf()`](https://wolski.github.io/prolfqua/reference/sim_build_models_logistf.md),
[`sim_make_model_lm()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lm.md),
[`sim_make_model_lmer()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lmer.md),
[`strategy_limma()`](https://wolski.github.io/prolfqua/reference/strategy_limma.md),
[`strategy_limpa()`](https://wolski.github.io/prolfqua/reference/strategy_limpa.md),
[`strategy_logistf()`](https://wolski.github.io/prolfqua/reference/strategy.md),
[`summary_ROPECA_median_p.scaled()`](https://wolski.github.io/prolfqua/reference/summary_ROPECA_median_p.scaled.md),
[`unregister_facade()`](https://wolski.github.io/prolfqua/reference/unregister_facade.md),
[`vcov.rfit_prolfqua()`](https://wolski.github.io/prolfqua/reference/vcov.rfit_prolfqua.md)

## Super class

[`prolfqua::ModelInterface`](https://wolski.github.io/prolfqua/reference/ModelInterface.md)
-\> `ModelLimma`

## Public fields

- `fit`:

  limma MArrayLM object from lmFit

- `design`:

  design matrix

- `formula`:

  model formula

- `subject_id`:

  protein ID column name(s)

- `model_name`:

  model name

- `rowdata`:

  protein ID mapping from to_wide()\$rowdata

- `trend`:

  passed to eBayes

- `robust`:

  passed to eBayes

- `dummy_model`:

  one fitted lm for linfct extraction

- `p.adjust`:

  function to adjust p-values

- `imputed_proteins`:

  character vector of protein ids (matching
  `rownames(fit$coefficients)`) that were rescued via LOD imputation.
  `ContrastsLimma$get_contrasts()` reads this to tag refit rows with a
  `"_imputed"` modelName suffix so the rescue is visible downstream.
  Empty for builders that do not impute.

## Methods

### Public methods

- [`ModelLimma$new()`](#method-ModelLimma-new)

- [`ModelLimma$get_coefficients()`](#method-ModelLimma-get_coefficients)

- [`ModelLimma$get_anova()`](#method-ModelLimma-get_anova)

- [`ModelLimma$coef_histogram()`](#method-ModelLimma-coef_histogram)

- [`ModelLimma$coef_volcano()`](#method-ModelLimma-coef_volcano)

- [`ModelLimma$coef_pairs()`](#method-ModelLimma-coef_pairs)

- [`ModelLimma$anova_histogram()`](#method-ModelLimma-anova_histogram)

- [`ModelLimma$clone()`](#method-ModelLimma-clone)

------------------------------------------------------------------------

### Method [`new()`](https://rdrr.io/r/methods/new.html)

initialize ModelLimma

#### Usage

    ModelLimma$new(
      fit,
      design,
      formula,
      subject_id,
      model_name,
      rowdata,
      trend = FALSE,
      robust = FALSE,
      dummy_model = NULL,
      p.adjust = prolfqua::adjust_p_values,
      imputed_proteins = character(0)
    )

#### Arguments

- `fit`:

  limma MArrayLM from lmFit

- `design`:

  design matrix

- `formula`:

  model formula

- `subject_id`:

  protein ID column name(s)

- `model_name`:

  model name

- `rowdata`:

  protein ID mapping

- `trend`:

  passed to eBayes

- `robust`:

  passed to eBayes

- `dummy_model`:

  one fitted lm for linfct extraction

- `p.adjust`:

  function to adjust p-values

- `imputed_proteins`:

  character vector of LOD-rescued protein ids

------------------------------------------------------------------------

### Method `get_coefficients()`

return model coefficient table in long format

#### Usage

    ModelLimma$get_coefficients()

#### Returns

data.frame

------------------------------------------------------------------------

### Method `get_anova()`

return anova table (F-test per protein across all non-intercept
coefficients)

#### Usage

    ModelLimma$get_anova()

#### Returns

data.frame

------------------------------------------------------------------------

### Method `coef_histogram()`

histogram of model coefficient p-values

#### Usage

    ModelLimma$coef_histogram()

------------------------------------------------------------------------

### Method `coef_volcano()`

volcano plot of non-intercept coefficients

#### Usage

    ModelLimma$coef_volcano()

------------------------------------------------------------------------

### Method `coef_pairs()`

pairs plot of coefficients

#### Usage

    ModelLimma$coef_pairs()

------------------------------------------------------------------------

### Method `anova_histogram()`

histogram of ANOVA F-test p-values

#### Usage

    ModelLimma$anova_histogram(what = c("p.value", "FDR"))

#### Arguments

- `what`:

  show either "p.value" or "FDR"

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ModelLimma$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- sim_lfq_data_protein_config(Nprot = 50)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lProt <- LFQData$new(istar$data, istar$config)
lProt$rename_response("transformedIntensity")
strat <- strategy_limma("transformedIntensity ~ group_")
mod <- build_model_limma(lProt, strat)
#> Warning: Partial NA coefficients for 1 probe(s)

coeffs <- mod$get_coefficients()
head(coeffs)
#> # A tibble: 6 × 6
#>   protein_Id  factor      Estimate Std..Error t.value  Pr...t..
#>   <chr>       <chr>          <dbl>      <dbl>   <dbl>     <dbl>
#> 1 0EfVhX~7161 (Intercept)     20.3      0.482    42.1 5.89e-143
#> 2 0m5WN4~3543 (Intercept)     20.5      0.482    42.6 1.17e-144
#> 3 76k03k~9735 (Intercept)     20.2      0.481    41.9 3.21e-142
#> 4 7QuTub~5556 (Intercept)     22.8      0.482    47.4 7.08e-159
#> 5 7cbcrd~0495 (Intercept)     17.2      0.681    25.2 3.32e- 82
#> 6 7soopj~3451 (Intercept)     26.3      0.481    54.7 3.79e-179
anova_tbl <- mod$get_anova()
head(anova_tbl)
#> # A tibble: 6 × 5
#>   protein_Id  F.value  p.value factor                  FDR
#>   <chr>         <dbl>    <dbl> <chr>                 <dbl>
#> 1 0EfVhX~7161    9.95 4.81e- 5 group_B+group_Ctrl 2.14e- 4
#> 2 0m5WN4~3543   26.8  2.57e-12 group_B+group_Ctrl 1.80e-11
#> 3 76k03k~9735    1.57 2.09e- 1 group_B+group_Ctrl 5.38e- 1
#> 4 7QuTub~5556    1.40 2.46e- 1 group_B+group_Ctrl 5.48e- 1
#> 5 7cbcrd~0495    6.12 2.22e- 3 group_B+group_Ctrl 8.35e- 3
#> 6 7soopj~3451    1.51 2.21e- 1 group_B+group_Ctrl 5.41e- 1
mod$coef_histogram()
#> $plot
#> Warning: Removed 1 row containing non-finite outside the scale range (`stat_bin()`).

#> 
#> $name
#> [1] "Coef_Histogram_limma.pdf"
#> 
mod$coef_volcano()
#> $plot
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_point()`).

#> 
#> $name
#> [1] "Coef_VolcanoPlot_limma.pdf"
#> 
mod$anova_histogram()
#> $plot
#> Warning: Removed 1 row containing non-finite outside the scale range (`stat_bin()`).

#> 
#> $name
#> [1] "Anova_p.values_limma.pdf"
#> 
```
