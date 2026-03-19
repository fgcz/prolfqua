# Lmer contrast analysis facade

Lmer contrast analysis facade

Lmer contrast analysis facade

## Details

Encapsulates the pipeline:
[`strategy_lmer`](https://wolski.github.io/prolfqua/reference/strategy.md)
-\>
[`build_model`](https://wolski.github.io/prolfqua/reference/build_model.md)
-\>
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md)
-\>
[`ContrastsModerated`](https://wolski.github.io/prolfqua/reference/ContrastsModerated.md).

This facade requires data with hierarchy below the analysis subject, for
example peptide-level measurements nested within proteins.

## See also

Other modelling:
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsFirth`](https://wolski.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsFirthFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsLMMissingFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMMissingFacade.md),
[`ContrastsLimma`](https://wolski.github.io/prolfqua/reference/ContrastsLimma.md),
[`ContrastsLimmaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaFacade.md),
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

  ContrastsModerated object

## Methods

### Public methods

- [`ContrastsLmerFacade$new()`](#method-ContrastsLmerFacade-new)

- [`ContrastsLmerFacade$get_contrasts()`](#method-ContrastsLmerFacade-get_contrasts)

- [`ContrastsLmerFacade$get_Plotter()`](#method-ContrastsLmerFacade-get_Plotter)

- [`ContrastsLmerFacade$to_wide()`](#method-ContrastsLmerFacade-to_wide)

- [`ContrastsLmerFacade$clone()`](#method-ContrastsLmerFacade-clone)

------------------------------------------------------------------------

### Method `new()`

initialize

#### Usage

    ContrastsLmerFacade$new(lfqdata, modelstr, contrasts, ...)

#### Arguments

- `lfqdata`:

  LFQData object

- `modelstr`:

  model formula string (e.g. "~ group\_ + (1 \| peptide_Id)")

- `contrasts`:

  named character vector of contrasts

- `...`:

  passed to
  [`strategy_lmer`](https://wolski.github.io/prolfqua/reference/strategy.md)

------------------------------------------------------------------------

### Method `get_contrasts()`

get contrast results

#### Usage

    ContrastsLmerFacade$get_contrasts(...)

#### Arguments

- `...`:

  passed to ContrastsModerated\$get_contrasts

------------------------------------------------------------------------

### Method `get_Plotter()`

get ContrastsPlotter

#### Usage

    ContrastsLmerFacade$get_Plotter(...)

#### Arguments

- `...`:

  passed to ContrastsModerated\$get_Plotter

------------------------------------------------------------------------

### Method `to_wide()`

convert results to wide format

#### Usage

    ContrastsLmerFacade$to_wide(...)

#### Arguments

- `...`:

  passed to ContrastsModerated\$to_wide

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsLmerFacade$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
istar$config <- old2new(istar$config)
lfqdata <- LFQData$new(istar$data, istar$config)
lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#> Column added : log2_abundance
contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
fa <- ContrastsLmerFacade$new(
  lfqdata,
  "~ group_ + (1 | peptide_Id) + (1 | sampleName)",
  contrasts
)
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> Warning: There were 4 warnings in `dplyr::mutate()`.
#> The first warning was:
#> ℹ In argument: `linear_model = purrr::map(data,
#>   model_strategy$model_fun, pb = pb)`.
#> ℹ In group 2: `protein_Id = "7cbcrd~5725"`.
#> Caused by warning in `value[[3L]]()`:
#> ! WARN :Error: grouping factors must have > 1 sampled level
#> ℹ Run `dplyr::last_dplyr_warnings()` to see the 3 remaining
#>   warnings.
#> Joining with `by = join_by(protein_Id)`
head(fa$get_contrasts())
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`
#> # A tibble: 6 × 14
#>   facade modelName protein_Id contrast     diff std.error avgAbd statistic    df
#>   <chr>  <chr>     <chr>      <chr>       <dbl>     <dbl>  <dbl>     <dbl> <dbl>
#> 1 lmer   WaldTest… 0EfVhX~00… A_vs_Ct… -8.32e-4    0.0730   4.34   -0.0115  28.9
#> 2 lmer   WaldTest… BEJI92~52… A_vs_Ct…  3.22e-1    0.0832   4.22    2.81    11.6
#> 3 lmer   WaldTest… Fl4JiV~86… A_vs_Ct… -4.13e-2    0.0850   4.38   -0.503   39.5
#> 4 lmer   WaldTest… HvIpHG~90… A_vs_Ct… -3.72e-1    0.0616   4.40   -5.65    21.8
#> 5 lmer   WaldTest… JcKVfU~96… A_vs_Ct… -1.07e-1    0.0577   5.05   -1.88    79.8
#> 6 lmer   WaldTest… SGIVBl~57… A_vs_Ct…  3.07e-2    0.0695   4.68    0.452   61.0
#> # ℹ 5 more variables: p.value <dbl>, conf.low <dbl>, conf.high <dbl>,
#> #   sigma <dbl>, FDR <dbl>
fa$to_wide()
#> # A tibble: 6 × 5
#>   protein_Id  diff.A_vs_Ctrl p.value.A_vs_Ctrl FDR.A_vs_Ctrl statistic.A_vs_Ctrl
#>   <chr>                <dbl>             <dbl>         <dbl>               <dbl>
#> 1 0EfVhX~0087      -0.000832         0.991         0.991                 -0.0115
#> 2 BEJI92~5282       0.322            0.0161        0.0484                 2.81  
#> 3 Fl4JiV~8625      -0.0413           0.618         0.783                 -0.503 
#> 4 HvIpHG~9079      -0.372            0.0000113     0.0000680             -5.65  
#> 5 JcKVfU~9653      -0.107            0.0637        0.127                 -1.88  
#> 6 SGIVBl~5782       0.0307           0.653         0.783                  0.452 
```
