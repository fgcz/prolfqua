# Rfit rank-regression contrast analysis with LOD imputation facade

Rfit rank-regression contrast analysis with LOD imputation facade

Rfit rank-regression contrast analysis with LOD imputation facade

## Value

An R6 class generator.

## Details

The rank-based analogue of
[`ContrastsLMImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMImputeFacade.md).
Encapsulates the pipeline:
[`strategy_rfit`](https://wolski.github.io/prolfqua/reference/strategy.md)
-\>
[`build_model_impute`](https://wolski.github.io/prolfqua/reference/build_model_impute.md)
-\>
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md)
-\>
[`ContrastsModerated`](https://wolski.github.io/prolfqua/reference/ContrastsModerated.md),
using [`rfit`](https://rdrr.io/pkg/Rfit/man/rfit.html).

Proteins whose initial `rfit` fit fails, is singular, or has an
incomplete design are re-fitted after imputing missing values with the
limit of detection (LOD). Uncertainty is borrowed from successful donor
fits so it is not underestimated by the constant imputation.

Differences from
[`ContrastsLMImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMImputeFacade.md):

- Covariance borrowing is **vcov-only** (element-wise median of the
  donors' full named covariance matrices). `rfit` exposes no `lm`-style
  `cov.unscaled`, so the scalar-sigma borrowing mode is not available
  and is not exposed.

- Borrowing is **fail-hard** (`on_misalign = "fail"`): if the donor
  covariance matrices cannot be aligned by coefficient name, the rescue
  is skipped and those proteins remain in `get_missing()` rather than
  carrying an `lm`-specific approximation.

- No observation weights (`rfit` has no `weights` argument), matching
  plain `rfit`.

Rescued rows surface as `modelName = "WaldTest_moderated_imputed"` in
the output, the same tag used by
[`ContrastsLMImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMImputeFacade.md).

## See also

[`ContrastsLMImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMImputeFacade.md),
[`ContrastsRfitFacade`](https://wolski.github.io/prolfqua/reference/ContrastsRfitFacade.md),
[`build_model_impute`](https://wolski.github.io/prolfqua/reference/build_model_impute.md)

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

## Super class

[`prolfqua::ContrastsInterface`](https://wolski.github.io/prolfqua/reference/ContrastsInterface.md)
-\> `ContrastsRfitImputeFacade`

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

- [`ContrastsRfitImputeFacade$new()`](#method-ContrastsRfitImputeFacade-new)

- [`ContrastsRfitImputeFacade$get_contrasts()`](#method-ContrastsRfitImputeFacade-get_contrasts)

- [`ContrastsRfitImputeFacade$get_missing()`](#method-ContrastsRfitImputeFacade-get_missing)

- [`ContrastsRfitImputeFacade$get_Plotter()`](#method-ContrastsRfitImputeFacade-get_Plotter)

- [`ContrastsRfitImputeFacade$to_wide()`](#method-ContrastsRfitImputeFacade-to_wide)

- [`ContrastsRfitImputeFacade$clone()`](#method-ContrastsRfitImputeFacade-clone)

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

    ContrastsRfitImputeFacade$new(
      lfqdata,
      modelstr,
      contrasts,
      lod = NULL,
      df_method = c("observed", "borrowed"),
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

- `...`:

  passed to
  [`strategy_rfit`](https://wolski.github.io/prolfqua/reference/strategy.md)

------------------------------------------------------------------------

### Method `get_contrasts()`

get contrast results

#### Usage

    ContrastsRfitImputeFacade$get_contrasts(...)

#### Arguments

- `...`:

  passed to ContrastsModerated\$get_contrasts

------------------------------------------------------------------------

### Method `get_missing()`

get protein x contrast pairs that could not be estimated

#### Usage

    ContrastsRfitImputeFacade$get_missing()

------------------------------------------------------------------------

### Method `get_Plotter()`

get ContrastsPlotter

#### Usage

    ContrastsRfitImputeFacade$get_Plotter(...)

#### Arguments

- `...`:

  passed to ContrastsModerated\$get_Plotter

------------------------------------------------------------------------

### Method `to_wide()`

convert results to wide format

#### Usage

    ContrastsRfitImputeFacade$to_wide(...)

#### Arguments

- `...`:

  passed to ContrastsModerated\$to_wide

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsRfitImputeFacade$clone(deep = FALSE)

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
fa <- ContrastsRfitImputeFacade$new(lfqdata, "~ group_", contrasts)
#> Warning: There were 2 warnings in `dplyr::mutate()`.
#> The first warning was:
#> ℹ In argument: `linear_model = purrr::map(data, model_strategy$model_fun, pb =
#>   pb)`.
#> ℹ In group 20: `protein_Id = "LlQY7I~3486"`.
#> Caused by warning in `sqrt()`:
#> ! NaNs produced
#> ℹ Run `dplyr::last_dplyr_warnings()` to see the 1 remaining warning.
head(fa$get_contrasts())
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`
#> # A tibble: 6 × 14
#>   facade  modelName protein_Id contrast    diff std.error avgAbd statistic    df
#>   <chr>   <chr>     <chr>      <chr>      <dbl>     <dbl>  <dbl>     <dbl> <dbl>
#> 1 rfit_i… WaldTest… 0EfVhX~29… A_vs_Ct…  1.40       0.999   22.7    1.55    27.3
#> 2 rfit_i… WaldTest… 0m5WN4~67… A_vs_Ct…  0.0318     0.857   20.7    0.0304  25.3
#> 3 rfit_i… WaldTest… 7QuTub~61… A_vs_Ct…  0.909      1.15    16.7    0.580   20.3
#> 4 rfit_i… WaldTest… 7cbcrd~26… A_vs_Ct…  0.612      1.64    21.8    0.553   23.3
#> 5 rfit_i… WaldTest… 9VUkAq~34… A_vs_Ct…  0.768      2.04    19.9    0.787   26.3
#> 6 rfit_i… WaldTest… At886V~77… A_vs_Ct… -1.91       1.03    29.2   -2.11    27.3
#> # ℹ 5 more variables: p.value <dbl>, conf.low <dbl>, conf.high <dbl>,
#> #   sigma <dbl>, FDR <dbl>
fa$to_wide()
#> # A tibble: 30 × 5
#>    protein_Id diff.A_vs_Ctrl p.value.A_vs_Ctrl FDR.A_vs_Ctrl statistic.A_vs_Ctrl
#>    <chr>               <dbl>             <dbl>         <dbl>               <dbl>
#>  1 0EfVhX~29…         1.40              0.133          0.666              1.55  
#>  2 0m5WN4~67…         0.0318            0.976          1                  0.0304
#>  3 7QuTub~61…         0.909             0.568          1                  0.580 
#>  4 7cbcrd~26…         0.612             0.586          1                  0.553 
#>  5 9VUkAq~34…         0.768             0.439          1                  0.787 
#>  6 At886V~77…        -1.91              0.0443         0.443             -2.11  
#>  7 BEJI92~27…        -1.30              0.195          0.711             -1.33  
#>  8 CGzoYe~08…         0.867             0.383          1                  0.888 
#>  9 CtOJ9t~91…         0.107             0.906          1                  0.119 
#> 10 DoWup2~28…        -1.50              0.108          0.646             -1.66  
#> # ℹ 20 more rows
```
