# applies contrast computation using lmerTest::contest function

applies contrast computation using lmerTest::contest function

## Usage

``` r
compute_lmer_contrast(model, linfct, ddf = c("Satterthwaite", "Kenward-Roger"))
```

## Arguments

- model:

  mixed effects model

- linfct:

  linear function

- ddf:

  method to determine denominator degrees of freedom

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
[`ContrastsRfitImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsRfitImputeFacade.md),
[`ContrastsTable`](https://wolski.github.io/prolfqua/reference/ContrastsTable.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`LR_test()`](https://wolski.github.io/prolfqua/reference/LR_test.md),
[`Model`](https://wolski.github.io/prolfqua/reference/Model.md),
[`ModelFirth`](https://wolski.github.io/prolfqua/reference/ModelFirth.md),
[`ModelLimma`](https://wolski.github.io/prolfqua/reference/ModelLimma.md),
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
[`new_imputed_model()`](https://wolski.github.io/prolfqua/reference/new_imputed_model.md),
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

## Examples

``` r
mb <- sim_make_model_lmer("interaction")
#> Warning: Unknown or uninitialised column: `nr_peptides`.
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> Warning: There were 4 warnings in `dplyr::mutate()`.
#> The first warning was:
#> ℹ In argument: `linear_model = purrr::map(data, model_strategy$model_fun, pb =
#>   pb)`.
#> ℹ In group 2: `protein_Id = "7cbcrd~5725"`.
#> Caused by warning:
#> ! grouping factors must have > 1 sampled level
#> ℹ Run `dplyr::last_dplyr_warnings()` to see the 3 remaining warnings.
summary(mb)
#> Linear mixed model fit by REML. t-tests use Satterthwaite's method [
#> lmerModLmerTest]
#> Formula: self$formula
#>    Data: x
#> 
#> REML criterion at convergence: 193.5
#> 
#> Scaled residuals: 
#>      Min       1Q   Median       3Q      Max 
#> -1.69297 -0.80789  0.02554  0.77048  1.74487 
#> 
#> Random effects:
#>  Groups     Name        Variance Std.Dev.
#>  sampleName (Intercept) 0.000    0.000   
#>  peptide_Id (Intercept) 2.198    1.483   
#>  Residual               3.403    1.845   
#> Number of obs: 48, groups:  sampleName, 16; peptide_Id, 3
#> 
#> Fixed effects:
#>                        Estimate Std. Error      df t value Pr(>|t|)    
#> (Intercept)             16.7976     1.0082  3.1878  16.662  0.00033 ***
#> TreatmentB               5.9436     0.7531 42.0000   7.893 7.91e-10 ***
#> BackgroundZ              9.3826     0.7531 42.0000  12.459 1.07e-15 ***
#> TreatmentB:BackgroundZ  -9.3668     1.0650 42.0000  -8.795 4.48e-11 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Correlation of Fixed Effects:
#>             (Intr) TrtmnB BckgrZ
#> TreatmentB  -0.373              
#> BackgroundZ -0.373  0.500       
#> TrtmntB:BcZ  0.264 -0.707 -0.707
#> optimizer (nloptwrap) convergence code: 0 (OK)
#> boundary (singular) fit: see help('isSingular')
#> 

linfct <- linfct_from_model(mb)
names(linfct)
#> [1] "linfct_factors"      "linfct_interactions"
compute_lmer_contrast(mb, linfct$linfct_factors)
#> # A tibble: 4 × 9
#>   lhs        estimate std.error    df statistic conf.low conf.high p.value sigma
#>   <chr>         <dbl>     <dbl> <dbl>     <dbl>    <dbl>     <dbl>   <dbl> <dbl>
#> 1 Backgroun…     19.8     0.935  2.37      21.1     16.3      23.3 9.42e-4  1.84
#> 2 Backgroun…     24.5     0.935  2.37      26.2     21.0      27.9 5.69e-4  1.84
#> 3 TreatmentA     21.5     0.935  2.37      23.0     18.0      25.0 7.74e-4  1.84
#> 4 TreatmentB     22.7     0.935  2.37      24.3     19.3      26.2 6.76e-4  1.84
compute_lmer_contrast(mb, linfct$linfct_interactions)
#> # A tibble: 4 × 9
#>   lhs        estimate std.error    df statistic conf.low conf.high p.value sigma
#>   <chr>         <dbl>     <dbl> <dbl>     <dbl>    <dbl>     <dbl>   <dbl> <dbl>
#> 1 Treatment…     16.8      1.01  3.19      16.7     13.7      19.9 3.30e-4  1.84
#> 2 Treatment…     26.2      1.01  3.19      26.0     23.1      29.3 8.08e-5  1.84
#> 3 Treatment…     22.7      1.01  3.19      22.6     19.6      25.8 1.26e-4  1.84
#> 4 Treatment…     22.8      1.01  3.19      22.6     19.7      25.9 1.26e-4  1.84
length(lme4::fixef(mb))
#> [1] 4
lmerTest::contest(mb, c( 0 ,1 , 0 , 0),joint = FALSE)
#>   Estimate Std. Error df  t value    lower    upper     Pr(>|t|)
#> 1 5.943641  0.7530708 42 7.892539 4.423882 7.463399 7.909703e-10
summary(mb)
#> Linear mixed model fit by REML. t-tests use Satterthwaite's method [
#> lmerModLmerTest]
#> Formula: self$formula
#>    Data: x
#> 
#> REML criterion at convergence: 193.5
#> 
#> Scaled residuals: 
#>      Min       1Q   Median       3Q      Max 
#> -1.69297 -0.80789  0.02554  0.77048  1.74487 
#> 
#> Random effects:
#>  Groups     Name        Variance Std.Dev.
#>  sampleName (Intercept) 0.000    0.000   
#>  peptide_Id (Intercept) 2.198    1.483   
#>  Residual               3.403    1.845   
#> Number of obs: 48, groups:  sampleName, 16; peptide_Id, 3
#> 
#> Fixed effects:
#>                        Estimate Std. Error      df t value Pr(>|t|)    
#> (Intercept)             16.7976     1.0082  3.1878  16.662  0.00033 ***
#> TreatmentB               5.9436     0.7531 42.0000   7.893 7.91e-10 ***
#> BackgroundZ              9.3826     0.7531 42.0000  12.459 1.07e-15 ***
#> TreatmentB:BackgroundZ  -9.3668     1.0650 42.0000  -8.795 4.48e-11 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Correlation of Fixed Effects:
#>             (Intr) TrtmnB BckgrZ
#> TreatmentB  -0.373              
#> BackgroundZ -0.373  0.500       
#> TrtmntB:BcZ  0.264 -0.707 -0.707
#> optimizer (nloptwrap) convergence code: 0 (OK)
#> boundary (singular) fit: see help('isSingular')
#> 
```
