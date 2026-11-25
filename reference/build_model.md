# Build protein models from data

Fits the model of `model_strategy` once per subject (by default per
protein) and returns the fits with their status. Failed fits are kept as
rows without a model.
[`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md)
calls it for the per-subject backends.

## Usage

``` r
build_model(
  data,
  model_strategy,
  subject_id = if (inherits(data, "LFQData")) {
     data$subject_id()
 } else {
    
    "protein_Id"
 },
  model_name = model_strategy$model_name
)
```

## Arguments

- data:

  data - a data frame or LFQData object

- model_strategy:

  model strategy object (e.g. from strategy_lmer or strategy_lm)

- subject_id:

  grouping variable

- model_name:

  model name

## Value

a object of class
[`Model`](https://fgcz.github.io/prolfqua/reference/Model.md)

## See also

[`strategy_lmer`](https://fgcz.github.io/prolfqua/reference/strategy.md),
[`strategy_lm`](https://fgcz.github.io/prolfqua/reference/strategy.md)

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

## Examples

``` r
D <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 20, weight_missing = 0.1)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 68 rows with NA response for feature and sample combinations absent from the input.
D$data$abundance |> is.na() |> sum()
#> [1] 68
D <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 20, weight_missing = 0.1, seed =3)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 59 rows with NA response for feature and sample combinations absent from the input.
D$data$abundance |> is.na() |> sum()
#> [1] 59
model_name <- "f_condtion_r_peptide"
formula_randomPeptide <-
  strategy_lmer("abundance  ~ group_ + (1 | peptide_Id) + (1 | sampleName)",
   model_name = model_name)


mod <- prolfqua::build_model(
 D$data,
 formula_randomPeptide,
 model_name = model_name,
 subject_id = D$config$hierarchy_keys_depth())
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> Warning: fit (f_condtion_r_peptide): 6 fits failed (e.g. 0GRprF~7339, AZPG26~9461, NNASGK~6205, quTD7H~9280, tCZCHm~5708, ...): Error: grouping factors must have > 1 sampled level
#> See $provenance().
#> fit (f_condtion_r_peptide): 14 singular fits (e.g. 4JK499~3111, 7IZdVV~6818, AoNKbb~3497, CibL2O~2149, EnbuYk~7811, ...): coefficients not estimable or fewer than 2 residual df
#> fit (f_condtion_r_peptide): 54 missing responses dropped (e.g. 4JK499~3111, 7IZdVV~6818, CibL2O~2149, EnbuYk~7811, JnEDsz~4232, ...): NA response (missing observation) left out of the fit
#> See $provenance().
aovtable <- mod$get_anova()

mod <- prolfqua::build_model(
 LFQData$new(D$data, D$config),
 formula_randomPeptide,
 model_name = model_name)
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> Warning: fit (f_condtion_r_peptide): 6 fits failed (e.g. 0GRprF~7339, AZPG26~9461, NNASGK~6205, quTD7H~9280, tCZCHm~5708, ...): Error: grouping factors must have > 1 sampled level
#> See $provenance().
#> fit (f_condtion_r_peptide): 14 singular fits (e.g. 4JK499~3111, 7IZdVV~6818, AoNKbb~3497, CibL2O~2149, EnbuYk~7811, ...): coefficients not estimable or fewer than 2 residual df
#> fit (f_condtion_r_peptide): 54 missing responses dropped (e.g. 4JK499~3111, 7IZdVV~6818, CibL2O~2149, EnbuYk~7811, JnEDsz~4232, ...): NA response (missing observation) left out of the fit
#> See $provenance().
head(mod$get_coefficients())
#> # A tibble: 6 × 9
#> # Groups:   protein_Id [2]
#>   protein_Id  factor      Estimate Std..Error    df t.value Pr...t.. isSingular
#>   <chr>       <chr>          <dbl>      <dbl> <dbl>   <dbl>    <dbl> <lgl>     
#> 1 4JK499~3111 (Intercept)   16.4        0.577 16.0   28.4   3.97e-15 TRUE      
#> 2 4JK499~3111 group_B        2.24       0.849 16.0    2.64  1.78e- 2 TRUE      
#> 3 4JK499~3111 group_Ctrl    -3.90       0.849 16.0   -4.59  3.03e- 4 TRUE      
#> 4 7IZdVV~6818 (Intercept)   24.0        1.69   7.95  14.2   6.19e- 7 TRUE      
#> 5 7IZdVV~6818 group_B        0.703      1.60  50.1    0.439 6.63e- 1 TRUE      
#> 6 7IZdVV~6818 group_Ctrl    -0.758      1.60  50.1   -0.473 6.38e- 1 TRUE      
#> # ℹ 1 more variable: nr_coef <int>
```
