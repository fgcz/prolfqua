# LM + missing-value imputation contrast analysis facade

LM + missing-value imputation contrast analysis facade

LM + missing-value imputation contrast analysis facade

## Details

Encapsulates the pipeline:
[`strategy_lm`](https://wolski.github.io/prolfqua/reference/strategy.md)
-\>
[`build_model`](https://wolski.github.io/prolfqua/reference/build_model.md)
-\>
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md)
-\> merge with
[`ContrastsMissing`](https://wolski.github.io/prolfqua/reference/ContrastsMissing.md)
-\>
[`ContrastsModerated`](https://wolski.github.io/prolfqua/reference/ContrastsModerated.md).

Proteins without a fitted model get their contrasts filled in from the
group-mean imputation method
([`ContrastsMissing`](https://wolski.github.io/prolfqua/reference/ContrastsMissing.md)).

## See also

Other modelling:
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsFirth`](https://wolski.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsFirthFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMFacade.md),
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
[`LR_test()`](https://wolski.github.io/prolfqua/reference/LR_test.md),
[`Model`](https://wolski.github.io/prolfqua/reference/Model.md),
[`ModelFirth`](https://wolski.github.io/prolfqua/reference/ModelFirth.md),
[`ModelLimma`](https://wolski.github.io/prolfqua/reference/ModelLimma.md),
[`build_contrast_analysis()`](https://wolski.github.io/prolfqua/reference/build_contrast_analysis.md),
[`build_model()`](https://wolski.github.io/prolfqua/reference/build_model.md),
[`build_model_glm_peptide()`](https://wolski.github.io/prolfqua/reference/build_model_glm_peptide.md),
[`build_model_glm_protein()`](https://wolski.github.io/prolfqua/reference/build_model_glm_protein.md),
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

## Public fields

- `model`:

  Model object

- `contrast`:

  ContrastsModerated object (merged with ContrastsMissing)

- `missing_contrast`:

  ContrastsMissing object

- `merged`:

  merged contrast result list from merge_contrasts_results

## Methods

### Public methods

- [`ContrastsLMMissingFacade$new()`](#method-ContrastsLMMissingFacade-new)

- [`ContrastsLMMissingFacade$get_contrasts()`](#method-ContrastsLMMissingFacade-get_contrasts)

- [`ContrastsLMMissingFacade$get_Plotter()`](#method-ContrastsLMMissingFacade-get_Plotter)

- [`ContrastsLMMissingFacade$to_wide()`](#method-ContrastsLMMissingFacade-to_wide)

- [`ContrastsLMMissingFacade$clone()`](#method-ContrastsLMMissingFacade-clone)

------------------------------------------------------------------------

### Method `new()`

initialize

#### Usage

    ContrastsLMMissingFacade$new(lfqdata, modelstr, contrasts, ...)

#### Arguments

- `lfqdata`:

  LFQData object

- `modelstr`:

  model formula string (e.g. "~ group\_")

- `contrasts`:

  named character vector of contrasts

- `...`:

  passed to
  [`strategy_lm`](https://wolski.github.io/prolfqua/reference/strategy.md)

------------------------------------------------------------------------

### Method `get_contrasts()`

get contrast results

#### Usage

    ContrastsLMMissingFacade$get_contrasts(...)

#### Arguments

- `...`:

  passed to ContrastsTable\$get_contrasts

------------------------------------------------------------------------

### Method `get_Plotter()`

get ContrastsPlotter

#### Usage

    ContrastsLMMissingFacade$get_Plotter(...)

#### Arguments

- `...`:

  passed to ContrastsTable\$get_Plotter

------------------------------------------------------------------------

### Method `to_wide()`

convert results to wide format

#### Usage

    ContrastsLMMissingFacade$to_wide(...)

#### Arguments

- `...`:

  passed to ContrastsTable\$to_wide

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsLMMissingFacade$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# ContrastsMissing requires protein-level data (hierarchyDepth == len(hierarchy_keys()))
istar <- sim_lfq_data_protein_config(Nprot = 30)
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
lfqdata <- LFQData$new(istar$data, istar$config)
lfqdata$rename_response("transformedIntensity")
contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
fa <- ContrastsLMMissingFacade$new(lfqdata, "~ group_", contrasts)
#> Joining with `by = join_by(protein_Id)`
#> determine linear functions:
#> Warning: linfct_matrix_contrasts: computed 0/2 contrasts; failed 2: A_vs_Ctrl, avg_A_vs_Ctrl. ℹ In argument: `A_vs_Ctrl = group_A - group_Ctrl`.
#> Caused by error:
#> ! object 'group_A' not found; ℹ In argument: `avg_A_vs_Ctrl = (group_A + group_Ctrl)/2`.
#> Caused by error:
#> ! object 'group_A' not found
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`
#> completing cases
#> A_vs_Ctrl=group_A - group_Ctrl
#> A_vs_Ctrl=group_A - group_Ctrl
#> A_vs_Ctrl=group_A - group_Ctrl
#> Joining with `by = join_by(protein_Id, contrast)`
#> Joining with `by = join_by(protein_Id, contrast)`
head(fa$get_contrasts())
#> # A tibble: 6 × 14
#>   facade  modelName protein_Id contrast    diff std.error avgAbd statistic    df
#>   <chr>   <fct>     <chr>      <chr>      <dbl>     <dbl>  <dbl>     <dbl> <dbl>
#> 1 lm_mis… WaldTest… 0EfVhX~29… A_vs_Ct…  1.24       0.731   22.6    1.63    31.8
#> 2 lm_mis… WaldTest… 0m5WN4~67… A_vs_Ct… -0.0361     0.614   20.8   -0.0412  29.8
#> 3 lm_mis… WaldTest… 7QuTub~61… A_vs_Ct… -0.680      0.806   16.6   -0.831   29.8
#> 4 lm_mis… WaldTest… 7cbcrd~26… A_vs_Ct…  0.704      0.718   22.0    0.929   31.8
#> 5 lm_mis… WaldTest… 9VUkAq~34… A_vs_Ct…  0.768      1.42    20.0    0.939   30.8
#> 6 lm_mis… WaldTest… At886V~77… A_vs_Ct… -1.86       0.706   29.1   -2.46    31.8
#> # ℹ 5 more variables: p.value <dbl>, conf.low <dbl>, conf.high <dbl>,
#> #   sigma <dbl>, FDR <dbl>
fa$to_wide()
#> # A tibble: 30 × 5
#>    protein_Id diff.A_vs_Ctrl p.value.A_vs_Ctrl FDR.A_vs_Ctrl statistic.A_vs_Ctrl
#>    <chr>               <dbl>             <dbl>         <dbl>               <dbl>
#>  1 0EfVhX~29…         1.24              0.113          0.509              1.63  
#>  2 0m5WN4~67…        -0.0361            0.967          0.967             -0.0412
#>  3 7QuTub~61…        -0.680             0.413          0.630             -0.831 
#>  4 7cbcrd~26…         0.704             0.360          0.630              0.929 
#>  5 9VUkAq~34…         0.768             0.355          0.630              0.939 
#>  6 At886V~77…        -1.86              0.0196         0.164             -2.46  
#>  7 BEJI92~27…        -0.721             0.348          0.630             -0.952 
#>  8 CGzoYe~08…        -0.389             0.611          0.770             -0.514 
#>  9 CtOJ9t~91…        -0.0717            0.925          0.958             -0.0946
#> 10 DoWup2~28…        -1.82              0.0226         0.164             -2.40  
#> # ℹ 20 more rows
```
