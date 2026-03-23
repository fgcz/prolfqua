# ROPECA contrast analysis facade

ROPECA contrast analysis facade

ROPECA contrast analysis facade

## Details

Encapsulates the pipeline:
[`strategy_lm`](https://wolski.github.io/prolfqua/reference/strategy.md)
-\>
[`build_model`](https://wolski.github.io/prolfqua/reference/build_model.md)
-\>
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md)
-\>
[`ContrastsROPECA`](https://wolski.github.io/prolfqua/reference/ContrastsROPECA.md).

ROPECA operates on peptide-level data and aggregates peptide-level
p-values to the protein level. The `lfqdata` object must contain
peptide-level data (i.e. `hierarchyDepth >= 2`).

## See also

Other modelling:
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsFirth`](https://wolski.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsFirthFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsLMImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMImputeFacade.md),
[`ContrastsLMMissingFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMMissingFacade.md),
[`ContrastsLimma`](https://wolski.github.io/prolfqua/reference/ContrastsLimma.md),
[`ContrastsLimmaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaFacade.md),
[`ContrastsLmerFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLmerFacade.md),
[`ContrastsMissing`](https://wolski.github.io/prolfqua/reference/ContrastsMissing.md),
[`ContrastsModerated`](https://wolski.github.io/prolfqua/reference/ContrastsModerated.md),
[`ContrastsModeratedDEqMS`](https://wolski.github.io/prolfqua/reference/ContrastsModeratedDEqMS.md),
[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.md),
[`ContrastsRLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsRLMFacade.md),
[`ContrastsROPECA`](https://wolski.github.io/prolfqua/reference/ContrastsROPECA.md),
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
[`build_model_impute()`](https://wolski.github.io/prolfqua/reference/build_model_impute.md),
[`build_model_limma()`](https://wolski.github.io/prolfqua/reference/build_model_limma.md),
[`build_model_logistf()`](https://wolski.github.io/prolfqua/reference/build_model_logistf.md),
[`compute_borrowed_variance()`](https://wolski.github.io/prolfqua/reference/compute_borrowed_variance.md),
[`contrasts_fisher_exact()`](https://wolski.github.io/prolfqua/reference/contrasts_fisher_exact.md),
[`get_anova_df()`](https://wolski.github.io/prolfqua/reference/get_anova_df.md),
[`get_complete_model_fit()`](https://wolski.github.io/prolfqua/reference/get_complete_model_fit.md),
[`get_p_values_pbeta()`](https://wolski.github.io/prolfqua/reference/get_p_values_pbeta.md),
[`group_label()`](https://wolski.github.io/prolfqua/reference/group_label.md),
[`impute_refit_singular()`](https://wolski.github.io/prolfqua/reference/impute_refit_singular.md),
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
[`my_contrast_V2()`](https://wolski.github.io/prolfqua/reference/my_contrast_V2.md),
[`my_glht()`](https://wolski.github.io/prolfqua/reference/my_glht.md),
[`new_lm_imputed()`](https://wolski.github.io/prolfqua/reference/new_lm_imputed.md),
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

  Model object (peptide-level)

- `contrast`:

  ContrastsROPECA object

- `.lfqdata`:

  stored reference to input LFQData

- `.contrast_names`:

  names of the requested contrasts

## Methods

### Public methods

- [`ContrastsROPECAFacade$new()`](#method-ContrastsROPECAFacade-new)

- [`ContrastsROPECAFacade$get_contrasts()`](#method-ContrastsROPECAFacade-get_contrasts)

- [`ContrastsROPECAFacade$get_missing()`](#method-ContrastsROPECAFacade-get_missing)

- [`ContrastsROPECAFacade$get_Plotter()`](#method-ContrastsROPECAFacade-get_Plotter)

- [`ContrastsROPECAFacade$to_wide()`](#method-ContrastsROPECAFacade-to_wide)

- [`ContrastsROPECAFacade$clone()`](#method-ContrastsROPECAFacade-clone)

------------------------------------------------------------------------

### Method `new()`

initialize

#### Usage

    ContrastsROPECAFacade$new(lfqdata, modelstr, contrasts, ...)

#### Arguments

- `lfqdata`:

  LFQData object with peptide-level data (hierarchyDepth \>= 2)

- `modelstr`:

  model formula string (e.g. "~ group\_")

- `contrasts`:

  named character vector of contrasts

- `...`:

  passed to
  [`strategy_lm`](https://wolski.github.io/prolfqua/reference/strategy.md)

------------------------------------------------------------------------

### Method `get_contrasts()`

Get contrast results with standardized column names.

ROPECA's `beta.based.significance` is mapped to `p.value` and
`FDR.beta.based.significance` to `FDR`.

Columns not directly produced by ROPECA are derived heuristically:

- `std.error = diff / statistic` (algebraic: t = estimate / SE)

- `sigma = mad.estimate` (MAD of peptide fold changes)

- `df = n_not_na` (number of contributing peptides)

- `conf.low/conf.high` via `diff ± qt(0.975, df) * |std.error|`

#### Usage

    ContrastsROPECAFacade$get_contrasts()

------------------------------------------------------------------------

### Method `get_missing()`

get protein × contrast pairs that could not be estimated

#### Usage

    ContrastsROPECAFacade$get_missing()

------------------------------------------------------------------------

### Method `get_Plotter()`

get ContrastsPlotter (uses standardized column names)

#### Usage

    ContrastsROPECAFacade$get_Plotter(FCthreshold = 2, FDRthreshold = 0.1)

#### Arguments

- `FCthreshold`:

  fold change threshold

- `FDRthreshold`:

  FDR threshold

------------------------------------------------------------------------

### Method `to_wide()`

convert results to wide format

#### Usage

    ContrastsROPECAFacade$to_wide(columns = c("p.value", "FDR", "statistic"))

#### Arguments

- `columns`:

  value columns to pivot, default `c("p.value", "FDR", "statistic")`

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsROPECAFacade$clone(deep = FALSE)

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
fa <- ContrastsROPECAFacade$new(lfqdata, "~ group_", contrasts)
#> Joining with `by = join_by(protein_Id, peptide_Id)`
head(fa$get_contrasts())
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, peptide_Id, contrast)`
#> # A tibble: 6 × 14
#> # Groups:   contrast [1]
#>   facade protein_Id  modelName contrast  avgAbd    diff        FDR statistic
#>   <chr>  <chr>       <chr>     <chr>      <dbl>   <dbl>      <dbl>     <dbl>
#> 1 ropeca 0EfVhX~0087 ROPECA    A_vs_Ctrl   4.27 -0.0742 0.0528         -1.75
#> 2 ropeca 7cbcrd~5725 ROPECA    A_vs_Ctrl   4.51  0.741  0.0000991       8.79
#> 3 ropeca 9VUkAq~4703 ROPECA    A_vs_Ctrl   4.47 -0.598  0.00000691    -12.7 
#> 4 ropeca BEJI92~5282 ROPECA    A_vs_Ctrl   4.23  0.277  0.00187         3.94
#> 5 ropeca CGzoYe~2147 ROPECA    A_vs_Ctrl   4.76 -0.310  0.0000374      -9.26
#> 6 ropeca DoWup2~5896 ROPECA    A_vs_Ctrl   4.43  0.295  0.00000138     14.7 
#> # ℹ 6 more variables: std.error <dbl>, df <int>, p.value <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, sigma <dbl>
fa$to_wide()
#> # A tibble: 10 × 5
#>    protein_Id diff.A_vs_Ctrl p.value.A_vs_Ctrl FDR.A_vs_Ctrl statistic.A_vs_Ctrl
#>    <chr>               <dbl>             <dbl>         <dbl>               <dbl>
#>  1 0EfVhX~00…        -0.0742       0.0423         0.0528                   -1.75
#>  2 7cbcrd~57…         0.741        0.0000496      0.0000991                 8.79
#>  3 9VUkAq~47…        -0.598        0.00000138     0.00000691              -12.7 
#>  4 BEJI92~52…         0.277        0.00112        0.00187                   3.94
#>  5 CGzoYe~21…        -0.310        0.0000150      0.0000374                -9.26
#>  6 DoWup2~58…         0.295        0.000000138    0.00000138               14.7 
#>  7 Fl4JiV~86…         0.0700       0.740          0.823                     2.08
#>  8 HvIpHG~90…        -0.384        0.00000389     0.0000130                -7.28
#>  9 JcKVfU~96…        -0.0634       0.00218        0.00312                  -1.87
#> 10 SGIVBl~57…        -0.122        0.993          0.993                    -3.88
```
