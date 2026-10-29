# Build limpa vooma model from aggregated LFQData

Takes protein-level LFQData produced by
[`AggregateLimpa`](https://fgcz.github.io/prolfqua/reference/AggregateLimpa.md)
(which contains standard error and observation count columns) and fits a
limpa vooma model via
[`limpa::voomaLmFitWithImputation`](https://rdrr.io/pkg/limpa/man/voomaLmFitWithImputation.html).
The standard errors from DPC quantification are used as a bivariate
predictor for the vooma precision weights, and the observation counts
identify imputed entries for degree-of-freedom correction.

## Usage

``` r
build_model_limpa(lfqdata, strategy, model_name = strategy$model_name)
```

## Arguments

- lfqdata:

  [`LFQData`](https://fgcz.github.io/prolfqua/reference/LFQData.md)
  object. Must have `config$nr_children` set (used to flag imputed
  observations via `nr_children == 0`). `config$opt_se` is optional:
  when set (e.g. from
  [`AggregateLimpa`](https://fgcz.github.io/prolfqua/reference/AggregateLimpa.md)),
  the SE column is used as a vooma precision-weight predictor; otherwise
  vooma is fit without an external predictor.

- strategy:

  [`StrategyLimpa`](https://fgcz.github.io/prolfqua/reference/StrategyLimpa.md)
  object

- model_name:

  model name (default from strategy)

## Value

[`ModelLimma`](https://fgcz.github.io/prolfqua/reference/ModelLimma.md)
object

## Details

Returns a
[`ModelLimma`](https://fgcz.github.io/prolfqua/reference/ModelLimma.md)
object (since the output is a standard `MArrayLM`), which can be passed
directly to
[`ContrastsLimma`](https://fgcz.github.io/prolfqua/reference/ContrastsLimma.md).

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
[`build_model_limma_voom_impute()`](https://fgcz.github.io/prolfqua/reference/build_model_limma_voom_impute.md),
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
if (requireNamespace("limpa", quietly = TRUE)) {
  istar <- sim_lfq_data_peptide_config(Nprot = 10)
  lfq <- LFQData$new(istar$data, istar$config)
  lfq <- lfq$get_Transformer()$log2()$lfq
  lfq <- AggregateLimpa$new(lfq, "protein")$aggregate()
  strat <- strategy_limpa("limpa ~ group_")
  mod <- build_model_limpa(lfq, strat)
  mod$model_name
}
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
#> Column added : log2_abundance
#> completing cases
#> [1] "limpa"
```
