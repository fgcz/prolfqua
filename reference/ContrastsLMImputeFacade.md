# LM contrast analysis with LOD imputation facade

LM contrast analysis with LOD imputation facade

LM contrast analysis with LOD imputation facade

## Value

An R6 class generator.

## Details

Encapsulates the pipeline:
[`strategy_lm`](https://wolski.github.io/prolfqua/reference/strategy.md)
-\>
[`build_model`](https://wolski.github.io/prolfqua/reference/build_model.md)
(with `impute = TRUE`) -\>
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md)
-\>
[`ContrastsModerated`](https://wolski.github.io/prolfqua/reference/ContrastsModerated.md).

Proteins whose initial lm fit fails or produces NA coefficients are
re-fitted after imputing missing values with the limit of detection
(LOD). The covariance matrix is borrowed from successful fits so that
the variance is not underestimated by the constant imputation.

## See also

Other modelling:
[`AnovaExtractor`](https://wolski.github.io/prolfqua/reference/AnovaExtractor.md),
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsDEqMSVoomFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSVoomFacade.md),
[`ContrastsFirth`](https://wolski.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsFirthFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsLMMissingFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMMissingFacade.md),
[`ContrastsLimma`](https://wolski.github.io/prolfqua/reference/ContrastsLimma.md),
[`ContrastsLimmaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaFacade.md),
[`ContrastsLimmaImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaImputeFacade.md),
[`ContrastsLimmaVoomFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaVoomFacade.md),
[`ContrastsLimmaVoomImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaVoomImputeFacade.md),
[`ContrastsLimpaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimpaFacade.md),
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
[`StrategyLimpa`](https://wolski.github.io/prolfqua/reference/StrategyLimpa.md),
[`StrategyLmer`](https://wolski.github.io/prolfqua/reference/StrategyLmer.md),
[`StrategyLogistf`](https://wolski.github.io/prolfqua/reference/StrategyLogistf.md),
[`StrategyRLM`](https://wolski.github.io/prolfqua/reference/StrategyRLM.md),
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
[`new_lm_imputed()`](https://wolski.github.io/prolfqua/reference/new_lm_imputed.md),
[`pivot_model_contrasts_to_wide()`](https://wolski.github.io/prolfqua/reference/pivot_model_contrasts_to_wide.md),
[`plot_lmer_peptide_predictions()`](https://wolski.github.io/prolfqua/reference/plot_lmer_peptide_predictions.md),
[`register_facade()`](https://wolski.github.io/prolfqua/reference/register_facade.md),
[`sim_build_models_lm()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lm.md),
[`sim_build_models_lmer()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lmer.md),
[`sim_build_models_logistf()`](https://wolski.github.io/prolfqua/reference/sim_build_models_logistf.md),
[`sim_make_model_lm()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lm.md),
[`sim_make_model_lmer()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lmer.md),
[`strategy_limma()`](https://wolski.github.io/prolfqua/reference/strategy_limma.md),
[`strategy_limpa()`](https://wolski.github.io/prolfqua/reference/strategy_limpa.md),
[`strategy_logistf()`](https://wolski.github.io/prolfqua/reference/strategy.md),
[`summary_ROPECA_median_p.scaled()`](https://wolski.github.io/prolfqua/reference/summary_ROPECA_median_p.scaled.md),
[`unregister_facade()`](https://wolski.github.io/prolfqua/reference/unregister_facade.md)

## Super class

[`prolfqua::ContrastsInterface`](https://wolski.github.io/prolfqua/reference/ContrastsInterface.md)
-\> `ContrastsLMImputeFacade`

## Public fields

- `model`:

  Model object (with imputed proteins)

- `contrast`:

  ContrastsModerated object

- `.lfqdata`:

  stored reference to input LFQData

- `.contrast_names`:

  names of the requested contrasts

## Methods

### Public methods

- [`ContrastsLMImputeFacade$new()`](#method-ContrastsLMImputeFacade-new)

- [`ContrastsLMImputeFacade$get_contrasts()`](#method-ContrastsLMImputeFacade-get_contrasts)

- [`ContrastsLMImputeFacade$get_missing()`](#method-ContrastsLMImputeFacade-get_missing)

- [`ContrastsLMImputeFacade$get_Plotter()`](#method-ContrastsLMImputeFacade-get_Plotter)

- [`ContrastsLMImputeFacade$to_wide()`](#method-ContrastsLMImputeFacade-to_wide)

- [`ContrastsLMImputeFacade$clone()`](#method-ContrastsLMImputeFacade-clone)

Inherited methods

- [`prolfqua::ContrastsInterface$column_description()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-column_description)
- [`prolfqua::ContrastsInterface$contrast_summary_table()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-contrast_summary_table)
- [`prolfqua::ContrastsInterface$extra_artifacts()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-extra_artifacts)
- [`prolfqua::ContrastsInterface$filter_significant()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-filter_significant)
- [`prolfqua::ContrastsInterface$get_config()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_config)
- [`prolfqua::ContrastsInterface$get_contrast_sides()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_contrast_sides)
- [`prolfqua::ContrastsInterface$get_ora()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_ora)
- [`prolfqua::ContrastsInterface$get_rank()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_rank)

------------------------------------------------------------------------

### Method `new()`

initialize

#### Usage

    ContrastsLMImputeFacade$new(
      lfqdata,
      modelstr,
      contrasts,
      lod = NULL,
      borrow_method = c("sigma", "vcov"),
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

- `borrow_method`:

  "sigma" borrows scalar sigma and uses per-protein (X'X)^-1; "vcov"
  borrows element-wise median of full vcov matrices

- `df_method`:

  "observed" uses max(n_observed - p, 1); "borrowed" uses median df from
  successful fits

- `weights`:

  column name for per-observation weights (default:
  `lfqdata$nr_children_col()`). Pass `NULL` for unweighted.

- `...`:

  passed to
  [`strategy_lm`](https://wolski.github.io/prolfqua/reference/strategy.md)

------------------------------------------------------------------------

### Method `get_contrasts()`

get contrast results

#### Usage

    ContrastsLMImputeFacade$get_contrasts(...)

#### Arguments

- `...`:

  passed to ContrastsModerated\$get_contrasts

------------------------------------------------------------------------

### Method `get_missing()`

get protein x contrast pairs that could not be estimated

#### Usage

    ContrastsLMImputeFacade$get_missing()

------------------------------------------------------------------------

### Method `get_Plotter()`

get ContrastsPlotter

#### Usage

    ContrastsLMImputeFacade$get_Plotter(...)

#### Arguments

- `...`:

  passed to ContrastsModerated\$get_Plotter

------------------------------------------------------------------------

### Method `to_wide()`

convert results to wide format

#### Usage

    ContrastsLMImputeFacade$to_wide(...)

#### Arguments

- `...`:

  passed to ContrastsModerated\$to_wide

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsLMImputeFacade$clone(deep = FALSE)

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
fa <- ContrastsLMImputeFacade$new(lfqdata, "~ group_", contrasts)
head(fa$get_contrasts())
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`
#> # A tibble: 6 × 14
#>   facade  modelName protein_Id contrast    diff std.error avgAbd statistic    df
#>   <chr>   <chr>     <chr>      <chr>      <dbl>     <dbl>  <dbl>     <dbl> <dbl>
#> 1 lm_imp… WaldTest… 0EfVhX~29… A_vs_Ct…  1.24       0.731   22.6    1.88   15.7 
#> 2 lm_imp… WaldTest… 0m5WN4~67… A_vs_Ct… -0.0361     0.614   20.8   -0.0334 13.7 
#> 3 lm_imp… WaldTest… 7QuTub~61… A_vs_Ct…  0.909      0.943   16.7    0.961   8.74
#> 4 lm_imp… WaldTest… 7cbcrd~26… A_vs_Ct…  0.612      1.08    21.9    0.574  11.7 
#> 5 lm_imp… WaldTest… 9VUkAq~34… A_vs_Ct…  0.768      1.42    20.0    0.710  14.7 
#> 6 lm_imp… WaldTest… At886V~77… A_vs_Ct… -1.86       0.706   29.1   -2.40   15.7 
#> # ℹ 5 more variables: p.value <dbl>, conf.low <dbl>, conf.high <dbl>,
#> #   sigma <dbl>, FDR <dbl>
fa$to_wide()
#> # A tibble: 30 × 5
#>    protein_Id diff.A_vs_Ctrl p.value.A_vs_Ctrl FDR.A_vs_Ctrl statistic.A_vs_Ctrl
#>    <chr>               <dbl>             <dbl>         <dbl>               <dbl>
#>  1 0EfVhX~29…         1.24             0.0791          0.474              1.88  
#>  2 0m5WN4~67…        -0.0361           0.974           0.974             -0.0334
#>  3 7QuTub~61…         0.909            0.362           0.776              0.961 
#>  4 7cbcrd~26…         0.612            0.577           0.883              0.574 
#>  5 9VUkAq~34…         0.768            0.489           0.883              0.710 
#>  6 At886V~77…        -1.86             0.0294          0.220             -2.40  
#>  7 BEJI92~27…        -1.30             0.221           0.603             -1.28  
#>  8 CGzoYe~08…         0.196            0.860           0.974              0.180 
#>  9 CtOJ9t~91…        -0.0717           0.926           0.974             -0.0938
#> 10 DoWup2~28…        -1.82             0.00762         0.114             -3.06  
#> # ℹ 20 more rows
```
