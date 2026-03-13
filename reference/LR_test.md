# Likelihood ratio test

Likelihood ratio test

## Usage

``` r
LR_test(
  modelProteinF,
  modelName,
  modelProteinF_Int,
  modelName_Int,
  subject_Id = "protein_Id",
  path = NULL
)
```

## Arguments

- modelProteinF:

  table with models (see build model)

- modelName:

  name of model

- modelProteinF_Int:

  reduced model

- modelName_Int:

  name of reduced model

- subject_Id:

  subject id typically Assession or protein_Id

- path:

  default NULL, set to a directory if you need to write diagnostic
  plots.

## See also

Other modelling:
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsFirth`](https://wolski.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsLMMissingFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMMissingFacade.md),
[`ContrastsLimma`](https://wolski.github.io/prolfqua/reference/ContrastsLimma.md),
[`ContrastsLimmaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaFacade.md),
[`ContrastsLmerFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLmerFacade.md),
[`ContrastsMissing`](https://wolski.github.io/prolfqua/reference/ContrastsMissing.md),
[`ContrastsModerated`](https://wolski.github.io/prolfqua/reference/ContrastsModerated.md),
[`ContrastsModeratedDEqMS`](https://wolski.github.io/prolfqua/reference/ContrastsModeratedDEqMS.md),
[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.md),
[`ContrastsProDA`](https://wolski.github.io/prolfqua/reference/ContrastsProDA.md),
[`ContrastsROPECA`](https://wolski.github.io/prolfqua/reference/ContrastsROPECA.md),
[`ContrastsROPECAFacade`](https://wolski.github.io/prolfqua/reference/ContrastsROPECAFacade.md),
[`ContrastsTable`](https://wolski.github.io/prolfqua/reference/ContrastsTable.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`Model`](https://wolski.github.io/prolfqua/reference/Model.md),
[`ModelFirth`](https://wolski.github.io/prolfqua/reference/ModelFirth.md),
[`ModelLimma`](https://wolski.github.io/prolfqua/reference/ModelLimma.md),
[`build_contrast_analysis()`](https://wolski.github.io/prolfqua/reference/build_contrast_analysis.md),
[`build_model()`](https://wolski.github.io/prolfqua/reference/build_model.md),
[`build_model_limma()`](https://wolski.github.io/prolfqua/reference/build_model_limma.md),
[`build_model_logistf()`](https://wolski.github.io/prolfqua/reference/build_model_logistf.md),
[`contrasts_fisher_exact()`](https://wolski.github.io/prolfqua/reference/contrasts_fisher_exact.md),
[`get_anova_df()`](https://wolski.github.io/prolfqua/reference/get_anova_df.md),
[`get_complete_model_fit()`](https://wolski.github.io/prolfqua/reference/get_complete_model_fit.md),
[`get_p_values_pbeta()`](https://wolski.github.io/prolfqua/reference/get_p_values_pbeta.md),
[`group_label()`](https://wolski.github.io/prolfqua/reference/group_label.md),
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
[`my_contest()`](https://wolski.github.io/prolfqua/reference/my_contest.md),
[`my_contrast()`](https://wolski.github.io/prolfqua/reference/my_contrast.md),
[`my_contrast_V1()`](https://wolski.github.io/prolfqua/reference/my_contrast_V1.md),
[`my_contrast_V2()`](https://wolski.github.io/prolfqua/reference/my_contrast_V2.md),
[`my_glht()`](https://wolski.github.io/prolfqua/reference/my_glht.md),
[`pivot_model_contrasts_2_Wide()`](https://wolski.github.io/prolfqua/reference/pivot_model_contrasts_2_Wide.md),
[`plot_lmer_peptide_predictions()`](https://wolski.github.io/prolfqua/reference/plot_lmer_peptide_predictions.md),
[`sim_build_models_lm()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lm.md),
[`sim_build_models_lmer()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lmer.md),
[`sim_build_models_logistf()`](https://wolski.github.io/prolfqua/reference/sim_build_models_logistf.md),
[`sim_make_model_lm()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lm.md),
[`sim_make_model_lmer()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lmer.md),
[`strategy_limma()`](https://wolski.github.io/prolfqua/reference/strategy_limma.md),
[`strategy_logistf()`](https://wolski.github.io/prolfqua/reference/strategy.md),
[`summary_ROPECA_median_p.scaled()`](https://wolski.github.io/prolfqua/reference/summary_ROPECA_median_p.scaled.md)

## Examples

``` r
data_2Factor <- prolfqua::sim_lfq_data_2Factor_config(
 Nprot = 200,
 with_missing = TRUE,
 weight_missing = 2)
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done

pMerged <- LFQData$new(data_2Factor$data, data_2Factor$config)

pMerged$config$get_response()
#> [1] "abundance"
pMerged$factors()
#> # A tibble: 16 × 4
#>    sample  sampleName Treatment Background
#>    <chr>   <chr>      <chr>     <chr>     
#>  1 A_V1    A_V1       A         X         
#>  2 A_V2    A_V2       A         X         
#>  3 A_V3    A_V3       A         X         
#>  4 A_V4    A_V4       A         X         
#>  5 B_V1    B_V1       B         X         
#>  6 B_V2    B_V2       B         X         
#>  7 B_V3    B_V3       B         X         
#>  8 B_V4    B_V4       B         X         
#>  9 C_V1    C_V1       B         Z         
#> 10 C_V2    C_V2       B         Z         
#> 11 C_V3    C_V3       B         Z         
#> 12 C_V4    C_V4       B         Z         
#> 13 Ctrl_V1 Ctrl_V1    A         Z         
#> 14 Ctrl_V2 Ctrl_V2    A         Z         
#> 15 Ctrl_V3 Ctrl_V3    A         Z         
#> 16 Ctrl_V4 Ctrl_V4    A         Z         

formula_condition_and_Batches <-
  prolfqua::strategy_lm("abundance ~ Treatment + Background")
modCB <- prolfqua::build_model(
  pMerged$data,
  formula_condition_and_Batches,
  subject_Id = pMerged$config$hierarchy_keys() )
#> Warning: There were 25 warnings in `dplyr::mutate()`.
#> The first warning was:
#> ℹ In argument: `linear_model = purrr::map(data, model_strategy$model_fun, pb =
#>   pb)`.
#> ℹ In group 33: `protein_Id = "Br6sVH~3679"`.
#> Caused by warning in `value[[3L]]()`:
#> ! WARN :Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]): contrasts can be applied only to factors with 2 or more levels
#> ℹ Run `dplyr::last_dplyr_warnings()` to see the 24 remaining warnings.
#> Joining with `by = join_by(protein_Id)`

formula_condition <-
  prolfqua::strategy_lm("abundance ~ Treatment")
modC <- prolfqua::build_model(
  pMerged$data,
  formula_condition,
  subject_Id = pMerged$config$hierarchy_keys() )
#> Warning: There were 19 warnings in `dplyr::mutate()`.
#> The first warning was:
#> ℹ In argument: `linear_model = purrr::map(data, model_strategy$model_fun, pb =
#>   pb)`.
#> ℹ In group 33: `protein_Id = "Br6sVH~3679"`.
#> Caused by warning in `value[[3L]]()`:
#> ! WARN :Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]): contrasts can be applied only to factors with 2 or more levels
#> ℹ Run `dplyr::last_dplyr_warnings()` to see the 18 remaining warnings.
#> Joining with `by = join_by(protein_Id)`

tmp <- LR_test(modCB$modelDF, "modCB", modC$modelDF, "modB")
hist(tmp$likelihood_ratio_test.pValue)

```
