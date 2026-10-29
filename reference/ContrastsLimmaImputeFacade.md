# Limma contrast analysis with LOD imputation facade

Limma contrast analysis with LOD imputation facade

Limma contrast analysis with LOD imputation facade

## Value

An R6 class generator.

## Details

Encapsulates the pipeline:
[`strategy_limma`](https://fgcz.github.io/prolfqua/reference/strategy_limma.md)
-\>
[`build_model_limma_impute`](https://fgcz.github.io/prolfqua/reference/build_model_limma_impute.md)
-\>
[`ContrastsLimma`](https://fgcz.github.io/prolfqua/reference/ContrastsLimma.md).

Proteins whose limma fit produces NA coefficients (typically from entire
missing groups) are recovered by imputing missing values with the limit
of detection (LOD) and refitting. The variance is borrowed from
successful fits and degrees of freedom are corrected so that inference
is not artificially precise from the constant imputation.

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

## Super classes

[`prolfqua::ContrastsInterface`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.md)
-\>
[`prolfqua::ContrastsFacadeBase`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.md)
-\> `ContrastsLimmaImputeFacade`

## Methods

### Public methods

- [`ContrastsLimmaImputeFacade$new()`](#method-ContrastsLimmaImputeFacade-new)

- [`ContrastsLimmaImputeFacade$clone()`](#method-ContrastsLimmaImputeFacade-clone)

Inherited methods

- [`prolfqua::ContrastsInterface$column_description()`](https://fgcz.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-column_description)
- [`prolfqua::ContrastsInterface$contrast_summary_table()`](https://fgcz.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-contrast_summary_table)
- [`prolfqua::ContrastsInterface$extra_artifacts()`](https://fgcz.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-extra_artifacts)
- [`prolfqua::ContrastsInterface$filter_significant()`](https://fgcz.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-filter_significant)
- [`prolfqua::ContrastsInterface$get_contrast_sides()`](https://fgcz.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_contrast_sides)
- [`prolfqua::ContrastsInterface$get_ora()`](https://fgcz.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_ora)
- [`prolfqua::ContrastsInterface$get_rank()`](https://fgcz.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_rank)
- [`prolfqua::ContrastsFacadeBase$get_Plotter()`](https://fgcz.github.io/prolfqua/html/ContrastsFacadeBase.html#method-ContrastsFacadeBase-get_Plotter)
- [`prolfqua::ContrastsFacadeBase$get_config()`](https://fgcz.github.io/prolfqua/html/ContrastsFacadeBase.html#method-ContrastsFacadeBase-get_config)
- [`prolfqua::ContrastsFacadeBase$get_contrasts()`](https://fgcz.github.io/prolfqua/html/ContrastsFacadeBase.html#method-ContrastsFacadeBase-get_contrasts)
- [`prolfqua::ContrastsFacadeBase$get_missing()`](https://fgcz.github.io/prolfqua/html/ContrastsFacadeBase.html#method-ContrastsFacadeBase-get_missing)
- [`prolfqua::ContrastsFacadeBase$to_wide()`](https://fgcz.github.io/prolfqua/html/ContrastsFacadeBase.html#method-ContrastsFacadeBase-to_wide)

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
  [`strategy_limma`](https://fgcz.github.io/prolfqua/reference/strategy_limma.md)
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
#> # A tibble: 6 × 16
#>   modelName    estimate_type protein_Id  contrast     diff   FDR
#>   <chr>        <chr>         <chr>       <chr>       <dbl> <dbl>
#> 1 limma_impute observed      0EfVhX~2956 A_vs_Ctrl  1.24   0.427
#> 2 limma_impute observed      0m5WN4~6730 A_vs_Ctrl -0.0361 0.971
#> 3 limma_impute observed      7QuTub~6175 A_vs_Ctrl  0.909  0.773
#> 4 limma_impute observed      7cbcrd~2687 A_vs_Ctrl  0.612  0.855
#> 5 limma_impute observed      9VUkAq~3402 A_vs_Ctrl  0.768  0.855
#> 6 limma_impute observed      At886V~7759 A_vs_Ctrl -1.86   0.212
#> # ℹ 10 more variables: std.error.unmoderated <dbl>, df.unmoderated <dbl>,
#> #   std.error <dbl>, statistic <dbl>, p.value <dbl>, sigma <dbl>, df <dbl>,
#> #   conf.low <dbl>, conf.high <dbl>, avgAbd <dbl>
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
