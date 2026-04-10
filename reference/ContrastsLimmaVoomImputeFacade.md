# Limma-voom contrast analysis with LOD imputation facade

Limma-voom contrast analysis with LOD imputation facade

Limma-voom contrast analysis with LOD imputation facade

## Details

Encapsulates the pipeline:
[`strategy_limma`](https://wolski.github.io/prolfqua/reference/strategy_limma.md)
-\>
[`build_model_limma_voom_impute`](https://wolski.github.io/prolfqua/reference/build_model_limma_voom_impute.md)
-\>
[`ContrastsLimma`](https://wolski.github.io/prolfqua/reference/ContrastsLimma.md).

Combines vooma precision weights with LOD imputation for proteins whose
fit produces NA coefficients (typically from entire missing groups).

## See also

Other modelling:
[`AnovaExtractor`](https://wolski.github.io/prolfqua/reference/AnovaExtractor.md),
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsDEqMSVoomFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSVoomFacade.md),
[`ContrastsFirth`](https://wolski.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsFirthFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsLMImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMImputeFacade.md),
[`ContrastsLMMissingFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMMissingFacade.md),
[`ContrastsLimma`](https://wolski.github.io/prolfqua/reference/ContrastsLimma.md),
[`ContrastsLimmaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaFacade.md),
[`ContrastsLimmaImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaImputeFacade.md),
[`ContrastsLimmaVoomFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaVoomFacade.md),
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
[`sim_build_models_lm()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lm.md),
[`sim_build_models_lmer()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lmer.md),
[`sim_build_models_logistf()`](https://wolski.github.io/prolfqua/reference/sim_build_models_logistf.md),
[`sim_make_model_lm()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lm.md),
[`sim_make_model_lmer()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lmer.md),
[`strategy_limma()`](https://wolski.github.io/prolfqua/reference/strategy_limma.md),
[`strategy_limpa()`](https://wolski.github.io/prolfqua/reference/strategy_limpa.md),
[`strategy_logistf()`](https://wolski.github.io/prolfqua/reference/strategy.md),
[`summary_ROPECA_median_p.scaled()`](https://wolski.github.io/prolfqua/reference/summary_ROPECA_median_p.scaled.md)

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

- [`ContrastsLimmaVoomImputeFacade$new()`](#method-ContrastsLimmaVoomImputeFacade-new)

- [`ContrastsLimmaVoomImputeFacade$get_contrasts()`](#method-ContrastsLimmaVoomImputeFacade-get_contrasts)

- [`ContrastsLimmaVoomImputeFacade$get_missing()`](#method-ContrastsLimmaVoomImputeFacade-get_missing)

- [`ContrastsLimmaVoomImputeFacade$get_Plotter()`](#method-ContrastsLimmaVoomImputeFacade-get_Plotter)

- [`ContrastsLimmaVoomImputeFacade$to_wide()`](#method-ContrastsLimmaVoomImputeFacade-to_wide)

- [`ContrastsLimmaVoomImputeFacade$clone()`](#method-ContrastsLimmaVoomImputeFacade-clone)

------------------------------------------------------------------------

### Method `new()`

initialize

#### Usage

    ContrastsLimmaVoomImputeFacade$new(
      lfqdata,
      modelstr,
      contrasts,
      lod = NULL,
      df_method = c("observed", "borrowed"),
      weights = lfqdata$config$nr_children,
      span = 0.5,
      plot = FALSE,
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
  `lfqdata$config$nr_children`). Pass `NULL` for unweighted.

- `span`:

  lowess smoother span for vooma trend (default 0.5)

- `plot`:

  logical; if TRUE, plot the mean-variance trend

- `...`:

  passed to
  [`strategy_limma`](https://wolski.github.io/prolfqua/reference/strategy_limma.md)
  (e.g. trend, robust)

------------------------------------------------------------------------

### Method `get_contrasts()`

get contrast results (rows with NA diff are filtered out)

#### Usage

    ContrastsLimmaVoomImputeFacade$get_contrasts(...)

#### Arguments

- `...`:

  passed to ContrastsLimma\$get_contrasts

------------------------------------------------------------------------

### Method `get_missing()`

get protein x contrast pairs that could not be estimated

#### Usage

    ContrastsLimmaVoomImputeFacade$get_missing()

------------------------------------------------------------------------

### Method `get_Plotter()`

get ContrastsPlotter

#### Usage

    ContrastsLimmaVoomImputeFacade$get_Plotter(...)

#### Arguments

- `...`:

  passed to ContrastsLimma\$get_Plotter

------------------------------------------------------------------------

### Method `to_wide()`

convert results to wide format

#### Usage

    ContrastsLimmaVoomImputeFacade$to_wide(...)

#### Arguments

- `...`:

  passed to ContrastsLimma\$to_wide

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsLimmaVoomImputeFacade$clone(deep = FALSE)

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
fa <- ContrastsLimmaVoomImputeFacade$new(lfqdata, "~ group_", contrasts)
#> Warning: Partial NA coefficients for 4 probe(s)
head(fa$get_contrasts())
#> # A tibble: 6 × 14
#>   facade modelName protein_Id contrast    diff   FDR std.error statistic p.value
#>   <chr>  <chr>     <chr>      <chr>      <dbl> <dbl>     <dbl>     <dbl>   <dbl>
#> 1 limma… limma     0EfVhX~29… A_vs_Ct…  1.24   0.325     0.617    2.00    0.0650
#> 2 limma… limma     0m5WN4~67… A_vs_Ct… -0.0361 0.973     1.03    -0.0352  0.973 
#> 3 limma… limma     7QuTub~61… A_vs_Ct…  0.909  0.694     0.813    1.12    0.301 
#> 4 limma… limma     7cbcrd~26… A_vs_Ct…  0.612  0.855     0.968    0.633   0.541 
#> 5 limma… limma     9VUkAq~34… A_vs_Ct…  0.768  0.855     1.15     0.671   0.514 
#> 6 limma… limma     At886V~77… A_vs_Ct… -1.86   0.240     0.822   -2.26    0.0400
#> # ℹ 5 more variables: sigma <dbl>, df <dbl>, conf.low <dbl>, conf.high <dbl>,
#> #   avgAbd <dbl>
fa$to_wide()
#> # A tibble: 30 × 5
#>    protein_Id diff.A_vs_Ctrl p.value.A_vs_Ctrl FDR.A_vs_Ctrl statistic.A_vs_Ctrl
#>    <chr>               <dbl>             <dbl>         <dbl>               <dbl>
#>  1 0EfVhX~29…         1.24             0.0650         0.325               2.00  
#>  2 0m5WN4~67…        -0.0361           0.973          0.973              -0.0352
#>  3 7QuTub~61…         0.909            0.301          0.694               1.12  
#>  4 7cbcrd~26…         0.612            0.541          0.855               0.633 
#>  5 9VUkAq~34…         0.768            0.514          0.855               0.671 
#>  6 At886V~77…        -1.86             0.0400         0.240              -2.26  
#>  7 BEJI92~27…        -1.30             0.154          0.462              -1.51  
#>  8 CGzoYe~08…         0.196            0.849          0.951               0.194 
#>  9 CtOJ9t~91…        -0.0717           0.931          0.973              -0.0882
#> 10 DoWup2~28…        -1.82             0.00796        0.0865             -3.09  
#> # ℹ 20 more rows
```
