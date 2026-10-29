# get linfct from model

get linfct from model

## Usage

``` r
linfct_from_model(m, as_list = TRUE)
```

## Arguments

- m:

  linear model

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
m <- sim_make_model_lm("factors")
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
linfct <- linfct_from_model(m, as_list = TRUE)

linfct$linfct_factors
#>             (Intercept) TreatmentB BackgroundZ
#> BackgroundX           1        0.5         0.0
#> BackgroundZ           1        0.5         1.0
#> TreatmentA            1        0.0         0.5
#> TreatmentB            1        1.0         0.5
linfct$linfct_interactions
#>                        (Intercept) TreatmentB BackgroundZ
#> TreatmentA:BackgroundX           1          0           0
#> TreatmentA:BackgroundZ           1          0           1
#> TreatmentB:BackgroundX           1          1           0
#> TreatmentB:BackgroundZ           1          1           1
lf <- matrix(
c(1, 1, 1, 1, 0.5, 0.5, 0, 1, 0, 1, 0.5, 0.5),
nrow = 4,
byrow = FALSE,
dimnames = list(c("BackgroundX", "BackgroundZ", "TreatmentA", "TreatmentB"),
                c("(Intercept)", "TreatmentB", "BackgroundZ"))
)
stopifnot(lf == linfct$linfct_factors)
m <- sim_make_model_lm("interaction")
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
linfct <- linfct_from_model(m)

m <- lm(Petal.Width ~ Species, data = iris)
linfct_from_model(m)
#> $linfct_factors
#>                   (Intercept) Speciesversicolor Speciesvirginica
#> Speciessetosa               1                 0                0
#> Speciesversicolor           1                 1                0
#> Speciesvirginica            1                 0                1
#> 
#> $linfct_interactions
#>                   (Intercept) Speciesversicolor Speciesvirginica
#> Speciessetosa               1                 0                0
#> Speciesversicolor           1                 1                0
#> Speciesvirginica            1                 0                1
#> 
xx <- data.frame( Y = 1:10 , Condition = c(rep("a",5), rep("b",5)) )
m <- lm(Y ~ Condition, data = xx)
linfct_from_model(m)
#> $linfct_factors
#>            (Intercept) Conditionb
#> Conditiona           1          0
#> Conditionb           1          1
#> 
#> $linfct_interactions
#>            (Intercept) Conditionb
#> Conditiona           1          0
#> Conditionb           1          1
#> 
xx <- data.frame( Y = 1:10 , Condition = c(rep("a",5), rep("b.b",5)) )
m <- lm(Y ~ Condition, data = xx)
linfct_from_model(m)
#> $linfct_factors
#>              (Intercept) Conditionb.b
#> Conditiona             1            0
#> Conditionb.b           1            1
#> 
#> $linfct_interactions
#>              (Intercept) Conditionb.b
#> Conditiona             1            0
#> Conditionb.b           1            1
#> 
xx <- data.frame( Y = 1:10 , Condition = c(rep("a",5), rep("ab",5)) )
m <- lm(Y ~ Condition, data = xx)
linfct_from_model(m)
#> $linfct_factors
#>             (Intercept) Conditionab
#> Conditiona            1           0
#> Conditionab           1           1
#> 
#> $linfct_interactions
#>             (Intercept) Conditionab
#> Conditiona            1           0
#> Conditionab           1           1
#> 
```
