# Build limma-voom model with LOD imputation for failed proteins

Combines vooma precision weights with LOD imputation for proteins that
have entire missing groups (NA coefficients). Mirrors
[`build_model_limma_impute`](https://fgcz.github.io/prolfqua/reference/build_model_limma_impute.md)
but uses vooma weights.

## Usage

``` r
build_model_limma_voom_impute(
  lfqdata,
  strategy,
  model_name = paste0(strategy$model_name, "Imputed"),
  lod = NULL,
  df_method = c("observed", "borrowed"),
  span = 0.5,
  plot = FALSE
)
```

## Arguments

- lfqdata:

  an [`LFQData`](https://fgcz.github.io/prolfqua/reference/LFQData.md)
  object (aggregated to protein level)

- strategy:

  output of
  [`strategy_limma`](https://fgcz.github.io/prolfqua/reference/strategy_limma.md)

- model_name:

  name of model (default: strategy name + "Imputed")

- lod:

  numeric limit of detection; if NULL, auto-computed from data

- df_method:

  how to set degrees of freedom for imputed proteins: `"observed"`
  (default) uses `max(n_observed - p, 1)`; `"borrowed"` uses the median
  df from successful fits

- span:

  lowess smoother span for vooma trend (default 0.5)

- plot:

  logical; if TRUE, plot the mean-variance trend

## Value

a
[`ModelLimma`](https://fgcz.github.io/prolfqua/reference/ModelLimma.md)
object with a hybrid fit

## See also

Other modelling:
[`AnovaExtractor`](https://fgcz.github.io/prolfqua/reference/AnovaExtractor.md),
[`Contrasts`](https://fgcz.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsDEqMSVoomFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsDEqMSVoomFacade.md),
[`ContrastsFacadeBase`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.md),
[`ContrastsFirth`](https://fgcz.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsFirthFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsFirthNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsFirthNestedFacade.md),
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
[`ContrastsRfitFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsRfitFacade.md),
[`ContrastsRfitImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsRfitImputeFacade.md),
[`ContrastsTable`](https://fgcz.github.io/prolfqua/reference/ContrastsTable.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`LR_test()`](https://fgcz.github.io/prolfqua/reference/LR_test.md),
[`Model`](https://fgcz.github.io/prolfqua/reference/Model.md),
[`ModelFirth`](https://fgcz.github.io/prolfqua/reference/ModelFirth.md),
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
[`build_model_glm_peptide()`](https://fgcz.github.io/prolfqua/reference/build_model_glm_peptide.md),
[`build_model_glm_protein()`](https://fgcz.github.io/prolfqua/reference/build_model_glm_protein.md),
[`build_model_impute()`](https://fgcz.github.io/prolfqua/reference/build_model_impute.md),
[`build_model_limma()`](https://fgcz.github.io/prolfqua/reference/build_model_limma.md),
[`build_model_limma_impute()`](https://fgcz.github.io/prolfqua/reference/build_model_limma_impute.md),
[`build_model_limma_voom()`](https://fgcz.github.io/prolfqua/reference/build_model_limma_voom.md),
[`build_model_limpa()`](https://fgcz.github.io/prolfqua/reference/build_model_limpa.md),
[`build_model_logistf()`](https://fgcz.github.io/prolfqua/reference/build_model_logistf.md),
[`compute_borrowed_variance()`](https://fgcz.github.io/prolfqua/reference/compute_borrowed_variance.md),
[`compute_borrowed_variance_limma()`](https://fgcz.github.io/prolfqua/reference/compute_borrowed_variance_limma.md),
[`compute_contrast()`](https://fgcz.github.io/prolfqua/reference/compute_contrast.md),
[`compute_lmer_contrast()`](https://fgcz.github.io/prolfqua/reference/compute_lmer_contrast.md),
[`df.residual.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/df.residual.rfit_prolfqua.md),
[`get_anova_df()`](https://fgcz.github.io/prolfqua/reference/get_anova_df.md),
[`get_complete_model_fit()`](https://fgcz.github.io/prolfqua/reference/get_complete_model_fit.md),
[`get_p_values_pbeta()`](https://fgcz.github.io/prolfqua/reference/get_p_values_pbeta.md),
[`group_label()`](https://fgcz.github.io/prolfqua/reference/group_label.md),
[`impute_from_model()`](https://fgcz.github.io/prolfqua/reference/impute_from_model.md),
[`impute_refit_singular()`](https://fgcz.github.io/prolfqua/reference/impute_refit_singular.md),
[`linfct_from_model()`](https://fgcz.github.io/prolfqua/reference/linfct_from_model.md),
[`linfct_matrix_contrasts()`](https://fgcz.github.io/prolfqua/reference/linfct_matrix_contrasts.md),
[`list_facades()`](https://fgcz.github.io/prolfqua/reference/list_facades.md),
[`lookup_facade()`](https://fgcz.github.io/prolfqua/reference/lookup_facade.md),
[`merge_contrasts_results()`](https://fgcz.github.io/prolfqua/reference/merge_contrasts_results.md),
[`model_analyse()`](https://fgcz.github.io/prolfqua/reference/model_analyse.md),
[`model_summary()`](https://fgcz.github.io/prolfqua/reference/model_summary.md),
[`moderated_p_deqms()`](https://fgcz.github.io/prolfqua/reference/moderated_p_deqms.md),
[`moderated_p_deqms_long()`](https://fgcz.github.io/prolfqua/reference/moderated_p_deqms_long.md),
[`moderated_p_limma()`](https://fgcz.github.io/prolfqua/reference/moderated_p_limma.md),
[`new_imputed_model()`](https://fgcz.github.io/prolfqua/reference/new_imputed_model.md),
[`pivot_model_contrasts_to_wide()`](https://fgcz.github.io/prolfqua/reference/pivot_model_contrasts_to_wide.md),
[`register_facade()`](https://fgcz.github.io/prolfqua/reference/register_facade.md),
[`sigma.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/sigma.rfit_prolfqua.md),
[`sim_build_models_lm()`](https://fgcz.github.io/prolfqua/reference/sim_build_models_lm.md),
[`sim_build_models_lmer()`](https://fgcz.github.io/prolfqua/reference/sim_build_models_lmer.md),
[`sim_build_models_logistf()`](https://fgcz.github.io/prolfqua/reference/sim_build_models_logistf.md),
[`sim_make_model_lm()`](https://fgcz.github.io/prolfqua/reference/sim_make_model_lm.md),
[`sim_make_model_lmer()`](https://fgcz.github.io/prolfqua/reference/sim_make_model_lmer.md),
[`strategy_limma()`](https://fgcz.github.io/prolfqua/reference/strategy_limma.md),
[`strategy_limpa()`](https://fgcz.github.io/prolfqua/reference/strategy_limpa.md),
[`strategy_logistf()`](https://fgcz.github.io/prolfqua/reference/strategy.md),
[`summary_ROPECA_median_p_scaled()`](https://fgcz.github.io/prolfqua/reference/summary_ROPECA_median_p_scaled.md),
[`vcov.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/vcov.rfit_prolfqua.md)

## Examples

``` r
istar <- sim_lfq_data_protein_config(Nprot = 50, weight_missing = 0.5)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfqdata <- LFQData$new(istar$data, istar$config)
lfqdata$rename_response("transformedIntensity")

strat <- strategy_limma("transformedIntensity ~ group_")
mod <- build_model_limma_voom_impute(lfqdata, strat)
#> Warning: Partial NA coefficients for 2 probe(s)
mod$get_coefficients()
#> # A tibble: 147 × 6
#>    protein_Id  factor      Estimate Std..Error t.value  Pr...t..
#>    <chr>       <chr>          <dbl>      <dbl>   <dbl>     <dbl>
#>  1 0EfVhX~7161 (Intercept)     20.6      0.528    39.0 7.04e-121
#>  2 0m5WN4~3543 (Intercept)     20.8      0.529    39.4 4.55e-122
#>  3 76k03k~9735 (Intercept)     20.2      0.466    43.2 1.40e-132
#>  4 7QuTub~5556 (Intercept)     22.8      0.528    43.2 2.09e-132
#>  5 7cbcrd~0495 (Intercept)     16.7      0.184    90.6 1.32e-223
#>  6 7soopj~3451 (Intercept)     26.3      0.475    55.5 3.74e-162
#>  7 9VUkAq~8655 (Intercept)     22.2      0.508    43.8 4.98e-134
#>  8 At886V~0359 (Intercept)     17.1      0.826    20.7 2.96e- 60
#>  9 BEJI92~5483 (Intercept)     15.9      0.827    19.2 1.95e- 54
#> 10 CGzoYe~1248 (Intercept)     18.2      0.585    31.2 2.78e- 97
#> # ℹ 137 more rows
```
