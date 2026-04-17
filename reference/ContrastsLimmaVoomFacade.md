# Limma-voom contrast analysis facade

Limma-voom contrast analysis facade

Limma-voom contrast analysis facade

## Details

Encapsulates the pipeline:
[`strategy_limma`](https://wolski.github.io/prolfqua/reference/strategy_limma.md)
-\>
[`build_model_limma_voom`](https://wolski.github.io/prolfqua/reference/build_model_limma_voom.md)
-\>
[`ContrastsLimma`](https://wolski.github.io/prolfqua/reference/ContrastsLimma.md).

Uses vooma-style precision weights derived from a mean-variance trend,
optionally combined with external weights (e.g. peptide/precursor
counts).

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

  ModelLimma object

- `contrast`:

  ContrastsLimma object

- `.lfqdata`:

  stored reference to input LFQData

- `.contrast_names`:

  names of the requested contrasts

## Methods

### Public methods

- [`ContrastsLimmaVoomFacade$new()`](#method-ContrastsLimmaVoomFacade-new)

- [`ContrastsLimmaVoomFacade$get_contrasts()`](#method-ContrastsLimmaVoomFacade-get_contrasts)

- [`ContrastsLimmaVoomFacade$get_missing()`](#method-ContrastsLimmaVoomFacade-get_missing)

- [`ContrastsLimmaVoomFacade$get_Plotter()`](#method-ContrastsLimmaVoomFacade-get_Plotter)

- [`ContrastsLimmaVoomFacade$to_wide()`](#method-ContrastsLimmaVoomFacade-to_wide)

- [`ContrastsLimmaVoomFacade$clone()`](#method-ContrastsLimmaVoomFacade-clone)

------------------------------------------------------------------------

### Method `new()`

initialize

#### Usage

    ContrastsLimmaVoomFacade$new(
      lfqdata,
      modelstr,
      contrasts,
      weights = lfqdata$nr_children_col(),
      span = 0.5,
      plot = FALSE,
      ...
    )

#### Arguments

- `lfqdata`:

  LFQData object

- `modelstr`:

  model formula string (e.g. "~ group\_")

- `contrasts`:

  named character vector of contrasts

- `weights`:

  column name for per-observation weights (default:
  `lfqdata$nr_children_col()`). Pass `NULL` for unweighted.

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

    ContrastsLimmaVoomFacade$get_contrasts(...)

#### Arguments

- `...`:

  passed to ContrastsLimma\$get_contrasts

------------------------------------------------------------------------

### Method `get_missing()`

get protein x contrast pairs that could not be estimated

#### Usage

    ContrastsLimmaVoomFacade$get_missing()

------------------------------------------------------------------------

### Method `get_Plotter()`

get ContrastsPlotter

#### Usage

    ContrastsLimmaVoomFacade$get_Plotter(...)

#### Arguments

- `...`:

  passed to ContrastsLimma\$get_Plotter

------------------------------------------------------------------------

### Method `to_wide()`

convert results to wide format

#### Usage

    ContrastsLimmaVoomFacade$to_wide(...)

#### Arguments

- `...`:

  passed to ContrastsLimma\$to_wide

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsLimmaVoomFacade$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- sim_lfq_data_protein_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfqdata <- LFQData$new(istar$data, istar$config)
lfqdata$rename_response("transformedIntensity")
contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
fa <- ContrastsLimmaVoomFacade$new(lfqdata, "~ group_", contrasts)
#> Warning: Partial NA coefficients for 1 probe(s)
head(fa$get_contrasts())
#> # A tibble: 6 × 14
#>   facade     modelName protein_Id  contrast     diff     FDR std.error statistic
#>   <chr>      <chr>     <chr>       <chr>       <dbl>   <dbl>     <dbl>     <dbl>
#> 1 limma_voom limma     0EfVhX~0087 A_vs_Ctrl -2.62   0.00198     0.482   -5.42  
#> 2 limma_voom limma     7cbcrd~5725 A_vs_Ctrl  2.80   0.0590      0.985    2.84  
#> 3 limma_voom limma     9VUkAq~4703 A_vs_Ctrl  1.67   0.0590      0.565    2.95  
#> 4 limma_voom limma     BEJI92~5282 A_vs_Ctrl  0.424  0.922       1.81     0.234 
#> 5 limma_voom limma     CGzoYe~2147 A_vs_Ctrl -0.598  0.794       1.17    -0.512 
#> 6 limma_voom limma     Fl4JiV~8625 A_vs_Ctrl -0.0494 0.955       0.851   -0.0581
#> # ℹ 6 more variables: p.value <dbl>, sigma <dbl>, df <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, avgAbd <dbl>
fa$to_wide()
#> # A tibble: 10 × 5
#>    protein_Id diff.A_vs_Ctrl p.value.A_vs_Ctrl FDR.A_vs_Ctrl statistic.A_vs_Ctrl
#>    <chr>               <dbl>             <dbl>         <dbl>               <dbl>
#>  1 0EfVhX~00…        -2.62            0.000220       0.00198             -5.42  
#>  2 7cbcrd~57…         2.80            0.0197         0.0590               2.84  
#>  3 9VUkAq~47…         1.67            0.0148         0.0590               2.95  
#>  4 BEJI92~52…         0.424           0.819          0.922                0.234 
#>  5 CGzoYe~21…        -0.598           0.618          0.794               -0.512 
#>  6 DoWup2~58…        NA              NA             NA                   NA     
#>  7 Fl4JiV~86…        -0.0494          0.955          0.955               -0.0581
#>  8 HvIpHG~90…        -0.809           0.389          0.713               -0.897 
#>  9 JcKVfU~96…         0.642           0.466          0.713                0.754 
#> 10 SGIVBl~57…        -0.494           0.475          0.713               -0.740 
```
