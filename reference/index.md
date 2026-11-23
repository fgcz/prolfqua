# Package index

## Differential expression analysis

One call fits a model per protein and tests the contrasts; each method
is a facade class.

- [`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md)
  : Build a contrast analysis using one of several statistical methods
- [`print(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`summary(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`print(`*`<summary.ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`coef(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`vcov(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`confint(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`formula(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`nobs(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`fitted(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`predict(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`residuals(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`plot(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  : Inspect a fitted prolfqua analysis
- [`list_facades()`](https://fgcz.github.io/prolfqua/reference/list_facades.md)
  : List currently registered contrast facades
- [`lookup_facade()`](https://fgcz.github.io/prolfqua/reference/lookup_facade.md)
  : Look up a contrast facade by short name
- [`register_facade()`](https://fgcz.github.io/prolfqua/reference/register_facade.md)
  : Register a contrast facade class
- [`ContrastsFacadeBase`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.md)
  : Base class for contrast analysis facades
- [`FACADE_REGISTRY`](https://fgcz.github.io/prolfqua/reference/FACADE_REGISTRY.md)
  : Registry of available contrast facade classes
- [`ContrastsBinomialFirthNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsBinomialFirthNestedFacade.md)
  : Firth-penalized detection-count facade for nested input
- [`ContrastsBinomialNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsBinomialNestedFacade.md)
  : Quasibinomial detection-count facade for nested input
- [`ContrastsDEqMSFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsDEqMSFacade.md)
  : DEqMS contrast analysis facade
- [`ContrastsDEqMSVoomFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsDEqMSVoomFacade.md)
  : DEqMS contrast analysis with vooma weights facade
- [`ContrastsFirthFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsFirthFacade.md)
  : Firth logistic missingness contrast analysis facade
- [`ContrastsHurdleFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsHurdleFacade.md)
  : Hurdle contrasts: intensity and detection combined with Fisher's
  method
- [`ContrastsHurdleNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsHurdleNestedFacade.md)
  : Hurdle contrast facade for nested input
- [`ContrastsLMFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLMFacade.md)
  : LM contrast analysis facade
- [`ContrastsLMImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLMImputeFacade.md)
  : LM contrast analysis with LOD imputation facade
- [`ContrastsLimmaFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimmaFacade.md)
  : Limma contrast analysis facade
- [`ContrastsLimmaImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimmaImputeFacade.md)
  : Limma contrast analysis with LOD imputation facade
- [`ContrastsLimmaVoomFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimmaVoomFacade.md)
  : Limma-voom contrast analysis facade
- [`ContrastsLimmaVoomImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimmaVoomImputeFacade.md)
  : Limma-voom contrast analysis with LOD imputation facade
- [`ContrastsLimpaFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimpaFacade.md)
  : Limpa contrast analysis facade
- [`ContrastsLimpaNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimpaNestedFacade.md)
  : Limpa contrast analysis facade for nested input
- [`ContrastsLmerNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLmerNestedFacade.md)
  : Lmer (mixed-model) contrast analysis facade for nested input
- [`ContrastsRLMFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsRLMFacade.md)
  : RLM contrast analysis facade
- [`ContrastsROPECANestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECANestedFacade.md)
  : ROPECA contrast analysis facade for nested input
- [`ContrastsROPECAV2NestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAV2NestedFacade.md)
  : ROPECA V2 contrast analysis facade for nested input
- [`ContrastsRfitFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsRfitFacade.md)
  : Rfit rank-based regression contrast analysis facade
- [`ContrastsRfitImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsRfitImputeFacade.md)
  : Rfit rank-regression contrast analysis with LOD imputation facade
- [`merge_contrasts_results()`](https://fgcz.github.io/prolfqua/reference/merge_contrasts_results.md)
  : Merge contrast results coming from two different model.
- [`pivot_model_contrasts_to_wide()`](https://fgcz.github.io/prolfqua/reference/pivot_model_contrasts_to_wide.md)
  : Pivot long-format model contrasts to wide format
- [`impute_from_model()`](https://fgcz.github.io/prolfqua/reference/impute_from_model.md)
  : Fill missing responses with predictions from per-subject linear
  models

## Data container

LFQData and the objects returned by its get\_ methods.

- [`AggregateLimpa`](https://fgcz.github.io/prolfqua/reference/AggregateLimpa.md)
  : Aggregate or impute intensities with limpa's DPC quantification
- [`AggregateMedpolish`](https://fgcz.github.io/prolfqua/reference/AggregateMedpolish.md)
  : Aggregate peptide intensities by median polish
- [`AggregateRlm`](https://fgcz.github.io/prolfqua/reference/AggregateRlm.md)
  : Aggregate peptide intensities by robust linear regression
- [`AggregateTopN`](https://fgcz.github.io/prolfqua/reference/AggregateTopN.md)
  : Aggregate peptide intensities from the top N peptides
- [`LFQData`](https://fgcz.github.io/prolfqua/reference/LFQData.md) :
  LFQData R6 class
- [`LFQDataPlotter`](https://fgcz.github.io/prolfqua/reference/LFQDataPlotter.md)
  : LFQDataPlotter —- Create various visualization of the LFQdata
- [`LFQDataStats`](https://fgcz.github.io/prolfqua/reference/LFQDataStats.md)
  : Decorates LFQData with methods to compute statistics of interactions
- [`LFQDataSummariser`](https://fgcz.github.io/prolfqua/reference/LFQDataSummariser.md)
  : Summarize LFQData
- [`lfqdata_from_wide()`](https://fgcz.github.io/prolfqua/reference/lfqdata_from_wide.md)
  : Build LFQData from wide matrices
- [`lfqdata_to_summarized_experiment()`](https://fgcz.github.io/prolfqua/reference/lfqdata_to_summarized_experiment.md)
  : Convert an LFQData object to a SummarizedExperiment
- [`plot(`*`<LFQData>`*`)`](https://fgcz.github.io/prolfqua/reference/plot.LFQData.md)
  : Plot an LFQData object
- [`LFQDataTransformer`](https://fgcz.github.io/prolfqua/reference/LFQDataTransformer.md)
  : Decorate LFQData with Methods for transforming Intensities
- [`MissingHelpers`](https://fgcz.github.io/prolfqua/reference/MissingHelpers.md)
  : Group statistics and limit of detection

## Configuration

- [`AnalysisConfiguration`](https://fgcz.github.io/prolfqua/reference/AnalysisConfiguration.md)
  : Analysis Configuration

- [`ContrastConfiguration`](https://fgcz.github.io/prolfqua/reference/ContrastConfiguration.md)
  : Contrast Configuration

- [`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md)
  : Internal Functions by category

- [`list_to_contrast_configuration()`](https://fgcz.github.io/prolfqua/reference/list_to_contrast_configuration.md)
  : Rebuild a ContrastConfiguration from a list

- [`r6_extract_values()`](https://fgcz.github.io/prolfqua/reference/r6_extract_values.md)
  : Extract all value slots in an R6 object

- [`setup_analysis()`](https://fgcz.github.io/prolfqua/reference/setup_analysis.md)
  :

  Setup a tidy table compatible with a `AnalysisConfiguration`

- [`list_to_analysis_configuration()`](https://fgcz.github.io/prolfqua/reference/list_to_analysis_configuration.md)
  : Reconstruct a configuration from minimal YAML

- [`center_to_reference_cfg()`](https://fgcz.github.io/prolfqua/reference/center_to_reference_cfg.md)
  : Center data to a reference subset

## Contrast definitions

- [`annotation_add_contrasts()`](https://fgcz.github.io/prolfqua/reference/annotation_add_contrasts.md)
  : DRY function: process and export annotated contrasts
- [`generate_contrasts()`](https://fgcz.github.io/prolfqua/reference/generate_contrasts.md)
  : Combined generate_contrasts
- [`generate_contrasts_for_factor()`](https://fgcz.github.io/prolfqua/reference/generate_contrasts_for_factor.md)
  : Single-factor contrasts (pairwise comparisons)
- [`interaction_contrasts()`](https://fgcz.github.io/prolfqua/reference/interaction_contrasts.md)
  : Interaction contrasts (difference of differences)
- [`level_specific_contrasts()`](https://fgcz.github.io/prolfqua/reference/level_specific_contrasts.md)
  : Level-specific contrasts (per secondary level)
- [`main_effect_contrasts()`](https://fgcz.github.io/prolfqua/reference/main_effect_contrasts.md)
  : Main effects contrasts

## Models, strategies and contrasts

The steps behind the facades, for inspecting or changing one of them.

- [`AnovaExtractor`](https://fgcz.github.io/prolfqua/reference/AnovaExtractor.md)
  : R6 class for extracting ANOVA results as a data frame
- [`Contrasts`](https://fgcz.github.io/prolfqua/reference/Contrasts.md)
  : Estimate contrasts using Wald Test
- [`ContrastsDEqMSFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsDEqMSFacade.md)
  : DEqMS contrast analysis facade
- [`ContrastsDEqMSVoomFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsDEqMSVoomFacade.md)
  : DEqMS contrast analysis with vooma weights facade
- [`ContrastsFacadeBase`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.md)
  : Base class for contrast analysis facades
- [`ContrastsFirth`](https://fgcz.github.io/prolfqua/reference/ContrastsFirth.md)
  : Estimate contrasts using Wald Test
- [`ContrastsFirthFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsFirthFacade.md)
  : Firth logistic missingness contrast analysis facade
- [`ContrastsHurdleFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsHurdleFacade.md)
  : Hurdle contrasts: intensity and detection combined with Fisher's
  method
- [`ContrastsHurdleNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsHurdleNestedFacade.md)
  : Hurdle contrast facade for nested input
- [`ContrastsLMFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLMFacade.md)
  : LM contrast analysis facade
- [`ContrastsLMImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLMImputeFacade.md)
  : LM contrast analysis with LOD imputation facade
- [`ContrastsLimma`](https://fgcz.github.io/prolfqua/reference/ContrastsLimma.md)
  : Limma-based contrasts (direct limma pipeline)
- [`ContrastsLimmaFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimmaFacade.md)
  : Limma contrast analysis facade
- [`ContrastsLimmaImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimmaImputeFacade.md)
  : Limma contrast analysis with LOD imputation facade
- [`ContrastsLimmaVoomFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimmaVoomFacade.md)
  : Limma-voom contrast analysis facade
- [`ContrastsLimmaVoomImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimmaVoomImputeFacade.md)
  : Limma-voom contrast analysis with LOD imputation facade
- [`ContrastsLimpaFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimpaFacade.md)
  : Limpa contrast analysis facade
- [`ContrastsLimpaNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimpaNestedFacade.md)
  : Limpa contrast analysis facade for nested input
- [`ContrastsLmerNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLmerNestedFacade.md)
  : Lmer (mixed-model) contrast analysis facade for nested input
- [`ContrastsModerated`](https://fgcz.github.io/prolfqua/reference/ContrastsModerated.md)
  : Limma moderated contrasts
- [`ContrastsModeratedDEqMS`](https://fgcz.github.io/prolfqua/reference/ContrastsModeratedDEqMS.md)
  : DEqMS count-dependent moderated contrasts
- [`ContrastsPlotter`](https://fgcz.github.io/prolfqua/reference/ContrastsPlotter.md)
  : Plot contrasts
- [`ContrastsRLMFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsRLMFacade.md)
  : RLM contrast analysis facade
- [`ContrastsROPECA`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECA.md)
  : ROPECA reproducibility-optimization method
- [`ContrastsROPECANestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECANestedFacade.md)
  : ROPECA contrast analysis facade for nested input
- [`ContrastsROPECAV2`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAV2.md)
  : ROPECA V2: median peptide z-score with a correlation-aware null
- [`ContrastsROPECAV2NestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAV2NestedFacade.md)
  : ROPECA V2 contrast analysis facade for nested input
- [`ContrastsRfitFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsRfitFacade.md)
  : Rfit rank-based regression contrast analysis facade
- [`ContrastsRfitImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsRfitImputeFacade.md)
  : Rfit rank-regression contrast analysis with LOD imputation facade
- [`ContrastsTable`](https://fgcz.github.io/prolfqua/reference/ContrastsTable.md)
  : Hold results when contrasts are added
- [`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md)
  : Internal Functions by category
- [`Model`](https://fgcz.github.io/prolfqua/reference/Model.md) : R6
  class representing modelling result
- [`ModelFirth`](https://fgcz.github.io/prolfqua/reference/ModelFirth.md)
  : R6 class representing modelling result
- [`ModelLimma`](https://fgcz.github.io/prolfqua/reference/ModelLimma.md)
  : R6 class representing a limma modelling result
- [`StrategyLM`](https://fgcz.github.io/prolfqua/reference/StrategyLM.md)
  : Linear model strategy (R6 class)
- [`StrategyLimma`](https://fgcz.github.io/prolfqua/reference/StrategyLimma.md)
  : R6 class for limma modelling strategy
- [`StrategyLimpa`](https://fgcz.github.io/prolfqua/reference/StrategyLimpa.md)
  : R6 class for limpa modelling strategy
- [`StrategyLmer`](https://fgcz.github.io/prolfqua/reference/StrategyLmer.md)
  : Linear mixed-effects model strategy (R6 class)
- [`StrategyLogistf`](https://fgcz.github.io/prolfqua/reference/StrategyLogistf.md)
  : Firth's logistic regression strategy (R6 class)
- [`StrategyRLM`](https://fgcz.github.io/prolfqua/reference/StrategyRLM.md)
  : Robust linear model strategy (R6 class)
- [`StrategyRfit`](https://fgcz.github.io/prolfqua/reference/StrategyRfit.md)
  : Rank-based regression strategy (R6 class)
- [`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md)
  : Build a contrast analysis using one of several statistical methods
- [`build_model()`](https://fgcz.github.io/prolfqua/reference/build_model.md)
  : Build protein models from data
- [`build_model_glm_protein()`](https://fgcz.github.io/prolfqua/reference/build_model_glm_protein.md)
  : Build Firth logistic model for aggregated LFQData
- [`build_model_impute()`](https://fgcz.github.io/prolfqua/reference/build_model_impute.md)
  : Build protein models with LOD imputation for failed fits
- [`build_model_limma()`](https://fgcz.github.io/prolfqua/reference/build_model_limma.md)
  : Build limma model from LFQData
- [`build_model_limma_impute()`](https://fgcz.github.io/prolfqua/reference/build_model_limma_impute.md)
  : Build limma model with LOD imputation for failed proteins
- [`build_model_limma_voom()`](https://fgcz.github.io/prolfqua/reference/build_model_limma_voom.md)
  : Build limma model with vooma precision weights (proteomics)
- [`build_model_limma_voom_impute()`](https://fgcz.github.io/prolfqua/reference/build_model_limma_voom_impute.md)
  : Build limma-voom model with LOD imputation for failed proteins
- [`build_model_limpa()`](https://fgcz.github.io/prolfqua/reference/build_model_limpa.md)
  : Build limpa vooma model from aggregated LFQData
- [`df.residual(`*`<rfit_prolfqua>`*`)`](https://fgcz.github.io/prolfqua/reference/df.residual.rfit_prolfqua.md)
  : Residual degrees of freedom for an augmented rfit fit
- [`group_label()`](https://fgcz.github.io/prolfqua/reference/group_label.md)
  : Build group labels from two factor levels
- [`impute_from_model()`](https://fgcz.github.io/prolfqua/reference/impute_from_model.md)
  : Fill missing responses with predictions from per-subject linear
  models
- [`list_facades()`](https://fgcz.github.io/prolfqua/reference/list_facades.md)
  : List currently registered contrast facades
- [`lookup_facade()`](https://fgcz.github.io/prolfqua/reference/lookup_facade.md)
  : Look up a contrast facade by short name
- [`lr_test()`](https://fgcz.github.io/prolfqua/reference/lr_test.md) :
  Likelihood ratio test
- [`merge_contrasts_results()`](https://fgcz.github.io/prolfqua/reference/merge_contrasts_results.md)
  : Merge contrast results coming from two different model.
- [`nobs(`*`<rfit_prolfqua>`*`)`](https://fgcz.github.io/prolfqua/reference/nobs.rfit_prolfqua.md)
  : Number of observations in an augmented rfit fit
- [`pivot_model_contrasts_to_wide()`](https://fgcz.github.io/prolfqua/reference/pivot_model_contrasts_to_wide.md)
  : Pivot long-format model contrasts to wide format
- [`print(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`summary(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`print(`*`<summary.ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`coef(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`vcov(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`confint(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`formula(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`nobs(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`fitted(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`predict(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`residuals(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  [`plot(`*`<ContrastsFacadeBase>`*`)`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md)
  : Inspect a fitted prolfqua analysis
- [`register_facade()`](https://fgcz.github.io/prolfqua/reference/register_facade.md)
  : Register a contrast facade class
- [`sigma(`*`<rfit_prolfqua>`*`)`](https://fgcz.github.io/prolfqua/reference/sigma.rfit_prolfqua.md)
  : Scale estimate for an augmented rfit fit
- [`strategy_logistf()`](https://fgcz.github.io/prolfqua/reference/strategy.md)
  [`strategy_rfit()`](https://fgcz.github.io/prolfqua/reference/strategy.md)
  [`strategy_lmer()`](https://fgcz.github.io/prolfqua/reference/strategy.md)
  [`strategy_lm()`](https://fgcz.github.io/prolfqua/reference/strategy.md)
  [`strategy_rlm()`](https://fgcz.github.io/prolfqua/reference/strategy.md)
  : Create Firth's logistic regression strategy
- [`strategy_limma()`](https://fgcz.github.io/prolfqua/reference/strategy_limma.md)
  : Create limma modelling strategy
- [`strategy_limpa()`](https://fgcz.github.io/prolfqua/reference/strategy_limpa.md)
  : Create limpa modelling strategy
- [`vcov(`*`<rfit_prolfqua>`*`)`](https://fgcz.github.io/prolfqua/reference/vcov.rfit_prolfqua.md)
  : Named variance-covariance matrix for an augmented rfit fit
- [`ModelInterface`](https://fgcz.github.io/prolfqua/reference/ModelInterface.md)
  : R6 interface class representing modelling result
- [`ContrastsInterface`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.md)
  : Base class for all Contrasts classes
- [`strategy_binomial()`](https://fgcz.github.io/prolfqua/reference/strategy_binomial.md)
  : Create a quasibinomial detection-count strategy
- [`strategy_binomial_firth()`](https://fgcz.github.io/prolfqua/reference/strategy_binomial_firth.md)
  : Create a Firth-penalized quasibinomial detection-count strategy
- [`StrategyBinomial`](https://fgcz.github.io/prolfqua/reference/StrategyBinomial.md)
  : Quasibinomial detection-count strategy
- [`StrategyBinomialFirth`](https://fgcz.github.io/prolfqua/reference/StrategyBinomialFirth.md)
  : Firth-penalized quasibinomial detection-count strategy
- [`adjust_p_values()`](https://fgcz.github.io/prolfqua/reference/adjust_p_values.md)
  : Adjust p-values for multiple testing

## Statistics and sample size

- [`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md)
  : Internal Functions by category
- [`fdr_sig_level()`](https://fgcz.github.io/prolfqua/reference/fdr_sig_level.md)
  : Per-test significance level for FDR-controlled sample-size planning
- [`null_split_check()`](https://fgcz.github.io/prolfqua/reference/null_split_check.md)
  : Null-split calibration check of normalization and model
- [`pooled_variance_explained()`](https://fgcz.github.io/prolfqua/reference/pooled_variance_explained.md)
  : Variance explained by each term, pooled over features
- [`power_fdr_prior_n()`](https://fgcz.github.io/prolfqua/reference/power_fdr_prior_n.md)
  : Plan sample size from a pilot variance prior and a planned effect
  distribution
- [`summarize_variance_explained()`](https://fgcz.github.io/prolfqua/reference/summarize_variance_explained.md)
  : Variance explained by the design, per feature
- [`HierarchyCountsSample`](https://fgcz.github.io/prolfqua/reference/HierarchyCountsSample.md)
  : Count distinct elements for each level of hierarchy per sample
- [`nr_children_experiment()`](https://fgcz.github.io/prolfqua/reference/nr_children_experiment.md)
  : Max nr_children across samples per hierarchy unit
- [`nr_features_experiment()`](https://fgcz.github.io/prolfqua/reference/nr_features_experiment.md)
  : Count distinct child features per hierarchy unit
- [`nr_obs_sample()`](https://fgcz.github.io/prolfqua/reference/nr_obs_sample.md)
  : Sum of observed children per hierarchy unit and sample
- [`summarize_stats_factors()`](https://fgcz.github.io/prolfqua/reference/summarize_stats_factors.md)
  : Compute variance, sd and related statistics for all factor levels
- [`which_missing()`](https://fgcz.github.io/prolfqua/reference/which_missing.md)
  : Draw missing values for an intensity vector based on its values

## Plotting

- [`ContrastsPlotter`](https://fgcz.github.io/prolfqua/reference/ContrastsPlotter.md)
  : Plot contrasts
- [`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md)
  : Internal Functions by category
- [`theme_prolfqua()`](https://fgcz.github.io/prolfqua/reference/theme_prolfqua.md)
  : Prism-style ggplot theme used by all prolfqua plots
- [`upset_interaction_missing_stats()`](https://fgcz.github.io/prolfqua/reference/upset_interaction_missing_stats.md)
  : UpSetR plot from interaction_missing_stats
- [`volcano_plotly()`](https://fgcz.github.io/prolfqua/reference/volcano_plotly.md)
  : Interactive volcano plots with plotly

## Simulated data

- [`sim_lfq_data()`](https://fgcz.github.io/prolfqua/reference/sim_lfq_data.md)
  : Simulate protein-level data with two groups
- [`sim_lfq_data_2factor_config()`](https://fgcz.github.io/prolfqua/reference/sim_lfq_data_2factor_config.md)
  : Simulate data, protein, with config with 2 factors Treatment and
  Background
- [`sim_lfq_data_peptide_config()`](https://fgcz.github.io/prolfqua/reference/sim_lfq_data_peptide_config.md)
  : Simulate data, protein and peptide, with config
- [`sim_lfq_data_protein_config()`](https://fgcz.github.io/prolfqua/reference/sim_lfq_data_protein_config.md)
  : Simulate data, protein, with config

## Utilities

- [`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md)
  : Internal Functions by category
- [`get_uniprot_id_from_fasta_header()`](https://fgcz.github.io/prolfqua/reference/get_uniprot_id_from_fasta_header.md)
  : Extracts uniprot ID
- [`is_contaminant()`](https://fgcz.github.io/prolfqua/reference/is_contaminant.md)
  : Detect contaminant identifiers
- [`is_decoy()`](https://fgcz.github.io/prolfqua/reference/is_decoy.md)
  : Detect decoy / reverse-database identifiers
- [`table_facade()`](https://fgcz.github.io/prolfqua/reference/table_facade.md)
  : Table facade to easily switch implementations
- [`find_package_file()`](https://fgcz.github.io/prolfqua/reference/find_package_file.md)
  : Find a file stored in a package

## Data

- [`PACKAGE_DATA`](https://fgcz.github.io/prolfqua/reference/PACKAGE_DATA.md)
  : Package Data
- [`prolfqua_data()`](https://fgcz.github.io/prolfqua/reference/prolfqua_data.md)
  : Load a dataset from prolfqua
