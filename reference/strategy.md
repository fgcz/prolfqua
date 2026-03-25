# Create Firth's logistic regression strategy

Convenience wrapper that creates a
[`StrategyLogistf`](https://wolski.github.io/prolfqua/reference/StrategyLogistf.md)
object.

Convenience wrapper that creates a
[`StrategyLmer`](https://wolski.github.io/prolfqua/reference/StrategyLmer.md)
object.

Convenience wrapper that creates a
[`StrategyLM`](https://wolski.github.io/prolfqua/reference/StrategyLM.md)
object.

Convenience wrapper that creates a
[`StrategyRLM`](https://wolski.github.io/prolfqua/reference/StrategyRLM.md)
object.

## Usage

``` r
strategy_logistf(
  modelstr,
  model_name = "logistf",
  report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value",
    "moderated.p.value.adjusted"),
  test = "Chisq"
)

strategy_lmer(
  modelstr,
  model_name = "Model",
  report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value",
    "moderated.p.value.adjusted")
)

strategy_lm(
  modelstr,
  model_name = "Model",
  report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value",
    "moderated.p.value.adjusted"),
  weights = NULL
)

strategy_rlm(
  modelstr,
  model_name = "Model",
  report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value",
    "moderated.p.value.adjusted")
)
```

## Arguments

- modelstr:

  model formula

- model_name:

  name of model

- report_columns:

  columns to report

- test:

  type of test statistic to use (e.g. "Chisq")

- weights:

  optional character string naming a column in the data containing
  per-observation weights, passed to
  [`lm`](https://rdrr.io/r/stats/lm.html). Default `NULL` (unweighted).

## Value

a
[`StrategyLogistf`](https://wolski.github.io/prolfqua/reference/StrategyLogistf.md)
object

a
[`StrategyLmer`](https://wolski.github.io/prolfqua/reference/StrategyLmer.md)
object

a
[`StrategyLM`](https://wolski.github.io/prolfqua/reference/StrategyLM.md)
object

a
[`StrategyRLM`](https://wolski.github.io/prolfqua/reference/StrategyRLM.md)
object

## See also

Other modelling:
[`AnovaExtractor`](https://wolski.github.io/prolfqua/reference/AnovaExtractor.md),
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsFirth`](https://wolski.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsFirthFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsLMImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMImputeFacade.md),
[`ContrastsLMMissingFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMMissingFacade.md),
[`ContrastsLimma`](https://wolski.github.io/prolfqua/reference/ContrastsLimma.md),
[`ContrastsLimmaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaFacade.md),
[`ContrastsLmerFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLmerFacade.md),
[`ContrastsMissing`](https://wolski.github.io/prolfqua/reference/ContrastsMissing.md),
[`ContrastsModerated`](https://wolski.github.io/prolfqua/reference/ContrastsModerated.md),
[`ContrastsModeratedDEqMS`](https://wolski.github.io/prolfqua/reference/ContrastsModeratedDEqMS.md),
[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.md),
[`ContrastsRLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsRLMFacade.md),
[`ContrastsROPECA`](https://wolski.github.io/prolfqua/reference/ContrastsROPECA.md),
[`ContrastsROPECAFacade`](https://wolski.github.io/prolfqua/reference/ContrastsROPECAFacade.md),
[`ContrastsTable`](https://wolski.github.io/prolfqua/reference/ContrastsTable.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`LR_test()`](https://wolski.github.io/prolfqua/reference/LR_test.md),
[`Model`](https://wolski.github.io/prolfqua/reference/Model.md),
[`ModelFirth`](https://wolski.github.io/prolfqua/reference/ModelFirth.md),
[`ModelLimma`](https://wolski.github.io/prolfqua/reference/ModelLimma.md),
[`StrategyLM`](https://wolski.github.io/prolfqua/reference/StrategyLM.md),
[`StrategyLimma`](https://wolski.github.io/prolfqua/reference/StrategyLimma.md),
[`StrategyLmer`](https://wolski.github.io/prolfqua/reference/StrategyLmer.md),
[`StrategyLogistf`](https://wolski.github.io/prolfqua/reference/StrategyLogistf.md),
[`StrategyRLM`](https://wolski.github.io/prolfqua/reference/StrategyRLM.md),
[`build_contrast_analysis()`](https://wolski.github.io/prolfqua/reference/build_contrast_analysis.md),
[`build_model()`](https://wolski.github.io/prolfqua/reference/build_model.md),
[`build_model_glm_peptide()`](https://wolski.github.io/prolfqua/reference/build_model_glm_peptide.md),
[`build_model_glm_protein()`](https://wolski.github.io/prolfqua/reference/build_model_glm_protein.md),
[`build_model_impute()`](https://wolski.github.io/prolfqua/reference/build_model_impute.md),
[`build_model_limma()`](https://wolski.github.io/prolfqua/reference/build_model_limma.md),
[`build_model_logistf()`](https://wolski.github.io/prolfqua/reference/build_model_logistf.md),
[`compute_borrowed_variance()`](https://wolski.github.io/prolfqua/reference/compute_borrowed_variance.md),
[`compute_contrast()`](https://wolski.github.io/prolfqua/reference/compute_contrast.md),
[`compute_lmer_contrast()`](https://wolski.github.io/prolfqua/reference/compute_lmer_contrast.md),
[`contrasts_fisher_exact()`](https://wolski.github.io/prolfqua/reference/contrasts_fisher_exact.md),
[`get_anova_df()`](https://wolski.github.io/prolfqua/reference/get_anova_df.md),
[`get_complete_model_fit()`](https://wolski.github.io/prolfqua/reference/get_complete_model_fit.md),
[`get_p_values_pbeta()`](https://wolski.github.io/prolfqua/reference/get_p_values_pbeta.md),
[`group_label()`](https://wolski.github.io/prolfqua/reference/group_label.md),
[`impute_refit_singular()`](https://wolski.github.io/prolfqua/reference/impute_refit_singular.md),
[`isSingular_lm()`](https://wolski.github.io/prolfqua/reference/isSingular_lm.md),
[`linfct_all_possible_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_all_possible_contrasts.md),
[`linfct_factors_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_factors_contrasts.md),
[`linfct_from_model()`](https://wolski.github.io/prolfqua/reference/linfct_from_model.md),
[`linfct_matrix_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_matrix_contrasts.md),
[`merge_contrasts_results()`](https://wolski.github.io/prolfqua/reference/merge_contrasts_results.md),
[`model_analyse()`](https://wolski.github.io/prolfqua/reference/model_analyse.md),
[`model_summary()`](https://wolski.github.io/prolfqua/reference/model_summary.md),
[`moderated_p_deqms()`](https://wolski.github.io/prolfqua/reference/moderated_p_deqms.md),
[`moderated_p_deqms_long()`](https://wolski.github.io/prolfqua/reference/moderated_p_deqms_long.md),
[`moderated_p_limma()`](https://wolski.github.io/prolfqua/reference/moderated_p_limma.md),
[`moderated_p_limma_long()`](https://wolski.github.io/prolfqua/reference/moderated_p_limma_long.md),
[`new_lm_imputed()`](https://wolski.github.io/prolfqua/reference/new_lm_imputed.md),
[`pivot_model_contrasts_2_Wide()`](https://wolski.github.io/prolfqua/reference/pivot_model_contrasts_2_Wide.md),
[`plot_lmer_peptide_predictions()`](https://wolski.github.io/prolfqua/reference/plot_lmer_peptide_predictions.md),
[`sim_build_models_lm()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lm.md),
[`sim_build_models_lmer()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lmer.md),
[`sim_build_models_logistf()`](https://wolski.github.io/prolfqua/reference/sim_build_models_logistf.md),
[`sim_make_model_lm()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lm.md),
[`sim_make_model_lmer()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lmer.md),
[`strategy_limma()`](https://wolski.github.io/prolfqua/reference/strategy_limma.md),
[`summary_ROPECA_median_p.scaled()`](https://wolski.github.io/prolfqua/reference/summary_ROPECA_median_p.scaled.md)

Other modelling:
[`AnovaExtractor`](https://wolski.github.io/prolfqua/reference/AnovaExtractor.md),
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsFirth`](https://wolski.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsFirthFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsLMImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMImputeFacade.md),
[`ContrastsLMMissingFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMMissingFacade.md),
[`ContrastsLimma`](https://wolski.github.io/prolfqua/reference/ContrastsLimma.md),
[`ContrastsLimmaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaFacade.md),
[`ContrastsLmerFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLmerFacade.md),
[`ContrastsMissing`](https://wolski.github.io/prolfqua/reference/ContrastsMissing.md),
[`ContrastsModerated`](https://wolski.github.io/prolfqua/reference/ContrastsModerated.md),
[`ContrastsModeratedDEqMS`](https://wolski.github.io/prolfqua/reference/ContrastsModeratedDEqMS.md),
[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.md),
[`ContrastsRLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsRLMFacade.md),
[`ContrastsROPECA`](https://wolski.github.io/prolfqua/reference/ContrastsROPECA.md),
[`ContrastsROPECAFacade`](https://wolski.github.io/prolfqua/reference/ContrastsROPECAFacade.md),
[`ContrastsTable`](https://wolski.github.io/prolfqua/reference/ContrastsTable.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`LR_test()`](https://wolski.github.io/prolfqua/reference/LR_test.md),
[`Model`](https://wolski.github.io/prolfqua/reference/Model.md),
[`ModelFirth`](https://wolski.github.io/prolfqua/reference/ModelFirth.md),
[`ModelLimma`](https://wolski.github.io/prolfqua/reference/ModelLimma.md),
[`StrategyLM`](https://wolski.github.io/prolfqua/reference/StrategyLM.md),
[`StrategyLimma`](https://wolski.github.io/prolfqua/reference/StrategyLimma.md),
[`StrategyLmer`](https://wolski.github.io/prolfqua/reference/StrategyLmer.md),
[`StrategyLogistf`](https://wolski.github.io/prolfqua/reference/StrategyLogistf.md),
[`StrategyRLM`](https://wolski.github.io/prolfqua/reference/StrategyRLM.md),
[`build_contrast_analysis()`](https://wolski.github.io/prolfqua/reference/build_contrast_analysis.md),
[`build_model()`](https://wolski.github.io/prolfqua/reference/build_model.md),
[`build_model_glm_peptide()`](https://wolski.github.io/prolfqua/reference/build_model_glm_peptide.md),
[`build_model_glm_protein()`](https://wolski.github.io/prolfqua/reference/build_model_glm_protein.md),
[`build_model_impute()`](https://wolski.github.io/prolfqua/reference/build_model_impute.md),
[`build_model_limma()`](https://wolski.github.io/prolfqua/reference/build_model_limma.md),
[`build_model_logistf()`](https://wolski.github.io/prolfqua/reference/build_model_logistf.md),
[`compute_borrowed_variance()`](https://wolski.github.io/prolfqua/reference/compute_borrowed_variance.md),
[`compute_contrast()`](https://wolski.github.io/prolfqua/reference/compute_contrast.md),
[`compute_lmer_contrast()`](https://wolski.github.io/prolfqua/reference/compute_lmer_contrast.md),
[`contrasts_fisher_exact()`](https://wolski.github.io/prolfqua/reference/contrasts_fisher_exact.md),
[`get_anova_df()`](https://wolski.github.io/prolfqua/reference/get_anova_df.md),
[`get_complete_model_fit()`](https://wolski.github.io/prolfqua/reference/get_complete_model_fit.md),
[`get_p_values_pbeta()`](https://wolski.github.io/prolfqua/reference/get_p_values_pbeta.md),
[`group_label()`](https://wolski.github.io/prolfqua/reference/group_label.md),
[`impute_refit_singular()`](https://wolski.github.io/prolfqua/reference/impute_refit_singular.md),
[`isSingular_lm()`](https://wolski.github.io/prolfqua/reference/isSingular_lm.md),
[`linfct_all_possible_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_all_possible_contrasts.md),
[`linfct_factors_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_factors_contrasts.md),
[`linfct_from_model()`](https://wolski.github.io/prolfqua/reference/linfct_from_model.md),
[`linfct_matrix_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_matrix_contrasts.md),
[`merge_contrasts_results()`](https://wolski.github.io/prolfqua/reference/merge_contrasts_results.md),
[`model_analyse()`](https://wolski.github.io/prolfqua/reference/model_analyse.md),
[`model_summary()`](https://wolski.github.io/prolfqua/reference/model_summary.md),
[`moderated_p_deqms()`](https://wolski.github.io/prolfqua/reference/moderated_p_deqms.md),
[`moderated_p_deqms_long()`](https://wolski.github.io/prolfqua/reference/moderated_p_deqms_long.md),
[`moderated_p_limma()`](https://wolski.github.io/prolfqua/reference/moderated_p_limma.md),
[`moderated_p_limma_long()`](https://wolski.github.io/prolfqua/reference/moderated_p_limma_long.md),
[`new_lm_imputed()`](https://wolski.github.io/prolfqua/reference/new_lm_imputed.md),
[`pivot_model_contrasts_2_Wide()`](https://wolski.github.io/prolfqua/reference/pivot_model_contrasts_2_Wide.md),
[`plot_lmer_peptide_predictions()`](https://wolski.github.io/prolfqua/reference/plot_lmer_peptide_predictions.md),
[`sim_build_models_lm()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lm.md),
[`sim_build_models_lmer()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lmer.md),
[`sim_build_models_logistf()`](https://wolski.github.io/prolfqua/reference/sim_build_models_logistf.md),
[`sim_make_model_lm()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lm.md),
[`sim_make_model_lmer()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lmer.md),
[`strategy_limma()`](https://wolski.github.io/prolfqua/reference/strategy_limma.md),
[`summary_ROPECA_median_p.scaled()`](https://wolski.github.io/prolfqua/reference/summary_ROPECA_median_p.scaled.md)

Other modelling:
[`AnovaExtractor`](https://wolski.github.io/prolfqua/reference/AnovaExtractor.md),
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsFirth`](https://wolski.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsFirthFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsLMImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMImputeFacade.md),
[`ContrastsLMMissingFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMMissingFacade.md),
[`ContrastsLimma`](https://wolski.github.io/prolfqua/reference/ContrastsLimma.md),
[`ContrastsLimmaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaFacade.md),
[`ContrastsLmerFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLmerFacade.md),
[`ContrastsMissing`](https://wolski.github.io/prolfqua/reference/ContrastsMissing.md),
[`ContrastsModerated`](https://wolski.github.io/prolfqua/reference/ContrastsModerated.md),
[`ContrastsModeratedDEqMS`](https://wolski.github.io/prolfqua/reference/ContrastsModeratedDEqMS.md),
[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.md),
[`ContrastsRLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsRLMFacade.md),
[`ContrastsROPECA`](https://wolski.github.io/prolfqua/reference/ContrastsROPECA.md),
[`ContrastsROPECAFacade`](https://wolski.github.io/prolfqua/reference/ContrastsROPECAFacade.md),
[`ContrastsTable`](https://wolski.github.io/prolfqua/reference/ContrastsTable.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`LR_test()`](https://wolski.github.io/prolfqua/reference/LR_test.md),
[`Model`](https://wolski.github.io/prolfqua/reference/Model.md),
[`ModelFirth`](https://wolski.github.io/prolfqua/reference/ModelFirth.md),
[`ModelLimma`](https://wolski.github.io/prolfqua/reference/ModelLimma.md),
[`StrategyLM`](https://wolski.github.io/prolfqua/reference/StrategyLM.md),
[`StrategyLimma`](https://wolski.github.io/prolfqua/reference/StrategyLimma.md),
[`StrategyLmer`](https://wolski.github.io/prolfqua/reference/StrategyLmer.md),
[`StrategyLogistf`](https://wolski.github.io/prolfqua/reference/StrategyLogistf.md),
[`StrategyRLM`](https://wolski.github.io/prolfqua/reference/StrategyRLM.md),
[`build_contrast_analysis()`](https://wolski.github.io/prolfqua/reference/build_contrast_analysis.md),
[`build_model()`](https://wolski.github.io/prolfqua/reference/build_model.md),
[`build_model_glm_peptide()`](https://wolski.github.io/prolfqua/reference/build_model_glm_peptide.md),
[`build_model_glm_protein()`](https://wolski.github.io/prolfqua/reference/build_model_glm_protein.md),
[`build_model_impute()`](https://wolski.github.io/prolfqua/reference/build_model_impute.md),
[`build_model_limma()`](https://wolski.github.io/prolfqua/reference/build_model_limma.md),
[`build_model_logistf()`](https://wolski.github.io/prolfqua/reference/build_model_logistf.md),
[`compute_borrowed_variance()`](https://wolski.github.io/prolfqua/reference/compute_borrowed_variance.md),
[`compute_contrast()`](https://wolski.github.io/prolfqua/reference/compute_contrast.md),
[`compute_lmer_contrast()`](https://wolski.github.io/prolfqua/reference/compute_lmer_contrast.md),
[`contrasts_fisher_exact()`](https://wolski.github.io/prolfqua/reference/contrasts_fisher_exact.md),
[`get_anova_df()`](https://wolski.github.io/prolfqua/reference/get_anova_df.md),
[`get_complete_model_fit()`](https://wolski.github.io/prolfqua/reference/get_complete_model_fit.md),
[`get_p_values_pbeta()`](https://wolski.github.io/prolfqua/reference/get_p_values_pbeta.md),
[`group_label()`](https://wolski.github.io/prolfqua/reference/group_label.md),
[`impute_refit_singular()`](https://wolski.github.io/prolfqua/reference/impute_refit_singular.md),
[`isSingular_lm()`](https://wolski.github.io/prolfqua/reference/isSingular_lm.md),
[`linfct_all_possible_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_all_possible_contrasts.md),
[`linfct_factors_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_factors_contrasts.md),
[`linfct_from_model()`](https://wolski.github.io/prolfqua/reference/linfct_from_model.md),
[`linfct_matrix_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_matrix_contrasts.md),
[`merge_contrasts_results()`](https://wolski.github.io/prolfqua/reference/merge_contrasts_results.md),
[`model_analyse()`](https://wolski.github.io/prolfqua/reference/model_analyse.md),
[`model_summary()`](https://wolski.github.io/prolfqua/reference/model_summary.md),
[`moderated_p_deqms()`](https://wolski.github.io/prolfqua/reference/moderated_p_deqms.md),
[`moderated_p_deqms_long()`](https://wolski.github.io/prolfqua/reference/moderated_p_deqms_long.md),
[`moderated_p_limma()`](https://wolski.github.io/prolfqua/reference/moderated_p_limma.md),
[`moderated_p_limma_long()`](https://wolski.github.io/prolfqua/reference/moderated_p_limma_long.md),
[`new_lm_imputed()`](https://wolski.github.io/prolfqua/reference/new_lm_imputed.md),
[`pivot_model_contrasts_2_Wide()`](https://wolski.github.io/prolfqua/reference/pivot_model_contrasts_2_Wide.md),
[`plot_lmer_peptide_predictions()`](https://wolski.github.io/prolfqua/reference/plot_lmer_peptide_predictions.md),
[`sim_build_models_lm()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lm.md),
[`sim_build_models_lmer()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lmer.md),
[`sim_build_models_logistf()`](https://wolski.github.io/prolfqua/reference/sim_build_models_logistf.md),
[`sim_make_model_lm()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lm.md),
[`sim_make_model_lmer()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lmer.md),
[`strategy_limma()`](https://wolski.github.io/prolfqua/reference/strategy_limma.md),
[`summary_ROPECA_median_p.scaled()`](https://wolski.github.io/prolfqua/reference/summary_ROPECA_median_p.scaled.md)

Other modelling:
[`AnovaExtractor`](https://wolski.github.io/prolfqua/reference/AnovaExtractor.md),
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsFirth`](https://wolski.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsFirthFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsLMImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMImputeFacade.md),
[`ContrastsLMMissingFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMMissingFacade.md),
[`ContrastsLimma`](https://wolski.github.io/prolfqua/reference/ContrastsLimma.md),
[`ContrastsLimmaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaFacade.md),
[`ContrastsLmerFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLmerFacade.md),
[`ContrastsMissing`](https://wolski.github.io/prolfqua/reference/ContrastsMissing.md),
[`ContrastsModerated`](https://wolski.github.io/prolfqua/reference/ContrastsModerated.md),
[`ContrastsModeratedDEqMS`](https://wolski.github.io/prolfqua/reference/ContrastsModeratedDEqMS.md),
[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.md),
[`ContrastsRLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsRLMFacade.md),
[`ContrastsROPECA`](https://wolski.github.io/prolfqua/reference/ContrastsROPECA.md),
[`ContrastsROPECAFacade`](https://wolski.github.io/prolfqua/reference/ContrastsROPECAFacade.md),
[`ContrastsTable`](https://wolski.github.io/prolfqua/reference/ContrastsTable.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`LR_test()`](https://wolski.github.io/prolfqua/reference/LR_test.md),
[`Model`](https://wolski.github.io/prolfqua/reference/Model.md),
[`ModelFirth`](https://wolski.github.io/prolfqua/reference/ModelFirth.md),
[`ModelLimma`](https://wolski.github.io/prolfqua/reference/ModelLimma.md),
[`StrategyLM`](https://wolski.github.io/prolfqua/reference/StrategyLM.md),
[`StrategyLimma`](https://wolski.github.io/prolfqua/reference/StrategyLimma.md),
[`StrategyLmer`](https://wolski.github.io/prolfqua/reference/StrategyLmer.md),
[`StrategyLogistf`](https://wolski.github.io/prolfqua/reference/StrategyLogistf.md),
[`StrategyRLM`](https://wolski.github.io/prolfqua/reference/StrategyRLM.md),
[`build_contrast_analysis()`](https://wolski.github.io/prolfqua/reference/build_contrast_analysis.md),
[`build_model()`](https://wolski.github.io/prolfqua/reference/build_model.md),
[`build_model_glm_peptide()`](https://wolski.github.io/prolfqua/reference/build_model_glm_peptide.md),
[`build_model_glm_protein()`](https://wolski.github.io/prolfqua/reference/build_model_glm_protein.md),
[`build_model_impute()`](https://wolski.github.io/prolfqua/reference/build_model_impute.md),
[`build_model_limma()`](https://wolski.github.io/prolfqua/reference/build_model_limma.md),
[`build_model_logistf()`](https://wolski.github.io/prolfqua/reference/build_model_logistf.md),
[`compute_borrowed_variance()`](https://wolski.github.io/prolfqua/reference/compute_borrowed_variance.md),
[`compute_contrast()`](https://wolski.github.io/prolfqua/reference/compute_contrast.md),
[`compute_lmer_contrast()`](https://wolski.github.io/prolfqua/reference/compute_lmer_contrast.md),
[`contrasts_fisher_exact()`](https://wolski.github.io/prolfqua/reference/contrasts_fisher_exact.md),
[`get_anova_df()`](https://wolski.github.io/prolfqua/reference/get_anova_df.md),
[`get_complete_model_fit()`](https://wolski.github.io/prolfqua/reference/get_complete_model_fit.md),
[`get_p_values_pbeta()`](https://wolski.github.io/prolfqua/reference/get_p_values_pbeta.md),
[`group_label()`](https://wolski.github.io/prolfqua/reference/group_label.md),
[`impute_refit_singular()`](https://wolski.github.io/prolfqua/reference/impute_refit_singular.md),
[`isSingular_lm()`](https://wolski.github.io/prolfqua/reference/isSingular_lm.md),
[`linfct_all_possible_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_all_possible_contrasts.md),
[`linfct_factors_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_factors_contrasts.md),
[`linfct_from_model()`](https://wolski.github.io/prolfqua/reference/linfct_from_model.md),
[`linfct_matrix_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_matrix_contrasts.md),
[`merge_contrasts_results()`](https://wolski.github.io/prolfqua/reference/merge_contrasts_results.md),
[`model_analyse()`](https://wolski.github.io/prolfqua/reference/model_analyse.md),
[`model_summary()`](https://wolski.github.io/prolfqua/reference/model_summary.md),
[`moderated_p_deqms()`](https://wolski.github.io/prolfqua/reference/moderated_p_deqms.md),
[`moderated_p_deqms_long()`](https://wolski.github.io/prolfqua/reference/moderated_p_deqms_long.md),
[`moderated_p_limma()`](https://wolski.github.io/prolfqua/reference/moderated_p_limma.md),
[`moderated_p_limma_long()`](https://wolski.github.io/prolfqua/reference/moderated_p_limma_long.md),
[`new_lm_imputed()`](https://wolski.github.io/prolfqua/reference/new_lm_imputed.md),
[`pivot_model_contrasts_2_Wide()`](https://wolski.github.io/prolfqua/reference/pivot_model_contrasts_2_Wide.md),
[`plot_lmer_peptide_predictions()`](https://wolski.github.io/prolfqua/reference/plot_lmer_peptide_predictions.md),
[`sim_build_models_lm()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lm.md),
[`sim_build_models_lmer()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lmer.md),
[`sim_build_models_logistf()`](https://wolski.github.io/prolfqua/reference/sim_build_models_logistf.md),
[`sim_make_model_lm()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lm.md),
[`sim_make_model_lmer()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lmer.md),
[`strategy_limma()`](https://wolski.github.io/prolfqua/reference/strategy_limma.md),
[`summary_ROPECA_median_p.scaled()`](https://wolski.github.io/prolfqua/reference/summary_ROPECA_median_p.scaled.md)

## Examples

``` r
tmp <- strategy_logistf("bin_resp ~ condition", model_name = "parallel design")
tmp$model_fun(get_formula = TRUE)
#> bin_resp ~ condition
#> <environment: 0x558827f25ff0>

istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 10, with_missing = TRUE,
  weight_missing = 0.5, seed = 3)
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
istar$data <- encode_bin_resp(istar$data, istar$config)
#> completing cases
istar <- LFQData$new(istar$data, istar$config)
df <- istar$summarize_hierarchy()
df2 <- df[df[[ncol(df)]] > 1, ]
istar2 <- istar$get_subset(df2)
#> Joining with `by = join_by(protein_Id)`
istar2$data |>
  dplyr::group_by(protein_Id) |>
  tidyr::nest() -> nestProtein
modelFunction <- strategy_logistf("bin_resp ~ group_ + peptide_Id",
  model_name = "random_example")
modelFunction$model_fun(nestProtein$data[[1]])
#> logistf::logistf(formula = self$formula, data = DFT, weights = Freq)
#> Model fitted by Penalized ML
#> Confidence intervals and p-values by Profile Likelihood 
#> 
#> Coefficients:
#>        (Intercept)            group_B         group_Ctrl peptide_IdFLq7LKTq 
#>       2.101899e+00       6.773389e-01      -6.663876e-01       3.440276e-16 
#> peptide_IdJYhOpuPH peptide_IdLiw5EMKP peptide_IdVcatZJTa peptide_IdjrLUqOjg 
#>      -1.068335e+00       3.432936e-16      -1.068335e+00       1.197563e+00 
#> peptide_Idq2jTaC1y 
#>      -1.445323e+00 
#> 
#> Likelihood ratio test=10.07111 on 8 df, p=0.2600707, n=84
#> 
modelFunction$model_fun(nestProtein$data[[4]])
#> logistf::logistf(formula = self$formula, data = DFT, weights = Freq)
#> Model fitted by Penalized ML
#> Confidence intervals and p-values by Profile Likelihood 
#> 
#> Coefficients:
#>        (Intercept)            group_B         group_Ctrl peptide_IdWcAw5ozd 
#>       8.375544e-03       8.360400e-11       9.365126e-01      -3.147690e-01 
#> peptide_IdgdnXrza3 peptide_IdxvlVt88v 
#>      -2.035802e-09      -6.303721e-01 
#> 
#> Likelihood ratio test=3.276542 on 5 df, p=0.657435, n=48
#> 
modelFunction <- strategy_lmer("abundanceC ~ group_ + (1|peptide_Id)",
  model_name = "random_example")
modelFunction$model_fun(get_formula = TRUE)
#> abundanceC ~ group_ + (1 | peptide_Id)
#> <environment: 0x558823cf53e8>
tmp <- strategy_lm("Intensity ~ condition", model_name = "parallel design")
tmp$model_fun(get_formula = TRUE)
#> Intensity ~ condition
#> <environment: 0x558823c89c50>
tmp$weights
#> NULL
tmp <- strategy_rlm("Intensity ~ condition", model_name = "parallel design")
tmp$model_fun(get_formula = TRUE)
#> Intensity ~ condition
#> <environment: 0x558823be49e0>
```
