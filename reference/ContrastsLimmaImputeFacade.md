# Limma contrast analysis with LOD imputation facade

Limma contrast analysis with LOD imputation facade

Limma contrast analysis with LOD imputation facade

## Value

An R6 class generator.

## Details

Encapsulates the pipeline:
[`strategy_limma`](https://wolski.github.io/prolfqua/reference/strategy_limma.md)
-\>
[`build_model_limma_impute`](https://wolski.github.io/prolfqua/reference/build_model_limma_impute.md)
-\>
[`ContrastsLimma`](https://wolski.github.io/prolfqua/reference/ContrastsLimma.md).

Proteins whose limma fit produces NA coefficients (typically from entire
missing groups) are recovered by imputing missing values with the limit
of detection (LOD) and refitting. The variance is borrowed from
successful fits and degrees of freedom are corrected so that inference
is not artificially precise from the constant imputation.

## See also

Other modelling:
[`AnovaExtractor`](https://wolski.github.io/prolfqua/reference/AnovaExtractor.md),
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsDEqMSVoomFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSVoomFacade.md),
[`ContrastsFacadeBase`](https://wolski.github.io/prolfqua/reference/ContrastsFacadeBase.md),
[`ContrastsFirth`](https://wolski.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsFirthFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsFirthNestedFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthNestedFacade.md),
[`ContrastsLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsLMImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMImputeFacade.md),
[`ContrastsLMMissingFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMMissingFacade.md),
[`ContrastsLimma`](https://wolski.github.io/prolfqua/reference/ContrastsLimma.md),
[`ContrastsLimmaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaFacade.md),
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

## Super classes

[`prolfqua::ContrastsInterface`](https://wolski.github.io/prolfqua/reference/ContrastsInterface.md)
-\>
[`prolfqua::ContrastsFacadeBase`](https://wolski.github.io/prolfqua/reference/ContrastsFacadeBase.md)
-\> `ContrastsLimmaImputeFacade`

## Public fields

- `model`:

  ModelLimma object (with imputed proteins)

- `contrast`:

  ContrastsLimma object

- `.lfqdata`:

  stored reference to input LFQData

- `.contrast_names`:

  names of the requested contrasts

## Methods

### Public methods

- [`ContrastsLimmaImputeFacade$new()`](#method-ContrastsLimmaImputeFacade-new)

- [`ContrastsLimmaImputeFacade$clone()`](#method-ContrastsLimmaImputeFacade-clone)

Inherited methods

- [`prolfqua::ContrastsInterface$column_description()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-column_description)
- [`prolfqua::ContrastsInterface$contrast_summary_table()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-contrast_summary_table)
- [`prolfqua::ContrastsInterface$extra_artifacts()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-extra_artifacts)
- [`prolfqua::ContrastsInterface$filter_significant()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-filter_significant)
- [`prolfqua::ContrastsInterface$get_config()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_config)
- [`prolfqua::ContrastsInterface$get_contrast_sides()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_contrast_sides)
- [`prolfqua::ContrastsInterface$get_ora()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_ora)
- [`prolfqua::ContrastsInterface$get_rank()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_rank)
- [`prolfqua::ContrastsFacadeBase$get_Plotter()`](https://wolski.github.io/prolfqua/html/ContrastsFacadeBase.html#method-ContrastsFacadeBase-get_Plotter)
- [`prolfqua::ContrastsFacadeBase$get_contrasts()`](https://wolski.github.io/prolfqua/html/ContrastsFacadeBase.html#method-ContrastsFacadeBase-get_contrasts)
- [`prolfqua::ContrastsFacadeBase$get_missing()`](https://wolski.github.io/prolfqua/html/ContrastsFacadeBase.html#method-ContrastsFacadeBase-get_missing)
- [`prolfqua::ContrastsFacadeBase$to_wide()`](https://wolski.github.io/prolfqua/html/ContrastsFacadeBase.html#method-ContrastsFacadeBase-to_wide)

------------------------------------------------------------------------

### Method `new()`

initialize

#### Usage

    ContrastsLimmaImputeFacade$new(
      lfqdata,
      modelstr,
      contrasts,
      lod = NULL,
      df_method = c("observed", "borrowed"),
      weights = lfqdata$nr_children_col(),
      ...
    )

#### Arguments

- `lfqdata`:

  LFQData object (aggregated to protein level)

- `modelstr`:

  model formula string (e.g. "~ group\_")

- `contrasts`:

  named character vector of contrasts

- `lod`:

  numeric limit of detection; if NULL, auto-computed from data

- `df_method`:

  "observed" uses max(n_observed - p, 1); "borrowed" uses median df from
  successful fits

- `weights`:

  column name for per-observation weights (default:
  `lfqdata$nr_children_col()`). Pass `NULL` for unweighted.

- `...`:

  passed to
  [`strategy_limma`](https://wolski.github.io/prolfqua/reference/strategy_limma.md)
  (e.g. trend, robust)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsLimmaImputeFacade$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- sim_lfq_data_protein_config(Nprot = 30, weight_missing = 0.5)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfqdata <- LFQData$new(istar$data, istar$config)
lfqdata$rename_response("transformedIntensity")
contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
fa <- ContrastsLimmaImputeFacade$new(lfqdata, "~ group_", contrasts)
#> Warning: Partial NA coefficients for 4 probe(s)
head(fa$get_contrasts())
#> # A tibble: 6 × 14
#>   modelName  estimate_type protein_Id contrast    diff   FDR std.error statistic
#>   <chr>      <chr>         <chr>      <chr>      <dbl> <dbl>     <dbl>     <dbl>
#> 1 limma_imp… observed      0EfVhX~29… A_vs_Ct…  1.24   0.427     0.666    1.86  
#> 2 limma_imp… observed      0m5WN4~67… A_vs_Ct… -0.0361 0.971     0.964   -0.0374
#> 3 limma_imp… observed      7QuTub~61… A_vs_Ct…  0.909  0.773     0.908    1.00  
#> 4 limma_imp… observed      7cbcrd~26… A_vs_Ct…  0.612  0.855     1.04     0.588 
#> 5 limma_imp… observed      9VUkAq~34… A_vs_Ct…  0.768  0.855     1.16     0.664 
#> 6 limma_imp… observed      At886V~77… A_vs_Ct… -1.86   0.212     0.744   -2.50  
#> # ℹ 6 more variables: p.value <dbl>, sigma <dbl>, df <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, avgAbd <dbl>
fa$to_wide()
#> # A tibble: 30 × 5
#>    protein_Id diff.A_vs_Ctrl p.value.A_vs_Ctrl FDR.A_vs_Ctrl statistic.A_vs_Ctrl
#>    <chr>               <dbl>             <dbl>         <dbl>               <dbl>
#>  1 0EfVhX~29…         1.24             0.0853          0.427              1.86  
#>  2 0m5WN4~67…        -0.0361           0.971           0.971             -0.0374
#>  3 7QuTub~61…         0.909            0.352           0.773              1.00  
#>  4 7cbcrd~26…         0.612            0.570           0.855              0.588 
#>  5 9VUkAq~34…         0.768            0.519           0.855              0.664 
#>  6 At886V~77…        -1.86             0.0258          0.212             -2.50  
#>  7 BEJI92~27…        -1.30             0.182           0.545             -1.41  
#>  8 CGzoYe~08…         0.196            0.850           0.971              0.194 
#>  9 CtOJ9t~91…        -0.0717           0.929           0.971             -0.0912
#> 10 DoWup2~28…        -1.82             0.00858         0.129             -3.07  
#> # ℹ 20 more rows
```
