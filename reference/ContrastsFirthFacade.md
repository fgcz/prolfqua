# Firth logistic missingness contrast analysis facade

Firth logistic missingness contrast analysis facade

Firth logistic missingness contrast analysis facade

## Details

Encapsulates the pipeline: encode missingness -\>
[`build_model_glm_protein`](https://wolski.github.io/prolfqua/reference/build_model_glm_protein.md)
or
[`build_model_glm_peptide`](https://wolski.github.io/prolfqua/reference/build_model_glm_peptide.md)
-\>
[`ContrastsFirth`](https://wolski.github.io/prolfqua/reference/ContrastsFirth.md).

The input may be aggregated protein-level data or nested peptide-level
data. The correct builder is chosen from the `LFQData` hierarchy
automatically.

Supports `options(prolfqua.vectorize = TRUE)` for faster contrast
computation. See
[`build_contrast_analysis`](https://wolski.github.io/prolfqua/reference/build_contrast_analysis.md)
for details.

## See also

Other modelling:
[`AnovaExtractor`](https://wolski.github.io/prolfqua/reference/AnovaExtractor.md),
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsDEqMSVoomFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSVoomFacade.md),
[`ContrastsFirth`](https://wolski.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsLMImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMImputeFacade.md),
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
[`new_lm_imputed()`](https://wolski.github.io/prolfqua/reference/new_lm_imputed.md),
[`pivot_model_contrasts_2_Wide()`](https://wolski.github.io/prolfqua/reference/pivot_model_contrasts_2_Wide.md),
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

  ModelFirth object

- `contrast`:

  ContrastsFirth object

- `.lfqdata`:

  stored reference to input LFQData

- `.contrast_names`:

  names of the requested contrasts

## Methods

### Public methods

- [`ContrastsFirthFacade$new()`](#method-ContrastsFirthFacade-new)

- [`ContrastsFirthFacade$get_contrasts()`](#method-ContrastsFirthFacade-get_contrasts)

- [`ContrastsFirthFacade$get_missing()`](#method-ContrastsFirthFacade-get_missing)

- [`ContrastsFirthFacade$get_Plotter()`](#method-ContrastsFirthFacade-get_Plotter)

- [`ContrastsFirthFacade$to_wide()`](#method-ContrastsFirthFacade-to_wide)

- [`ContrastsFirthFacade$clone()`](#method-ContrastsFirthFacade-clone)

------------------------------------------------------------------------

### Method `new()`

initialize

#### Usage

    ContrastsFirthFacade$new(lfqdata, modelstr, contrasts)

#### Arguments

- `lfqdata`:

  LFQData object

- `modelstr`:

  model formula string (e.g. "~ group\_")

- `contrasts`:

  named character vector of contrasts

------------------------------------------------------------------------

### Method `get_contrasts()`

get contrast results

#### Usage

    ContrastsFirthFacade$get_contrasts(...)

#### Arguments

- `...`:

  passed to ContrastsFirth\$get_contrasts

------------------------------------------------------------------------

### Method `get_missing()`

get protein × contrast pairs that could not be estimated

#### Usage

    ContrastsFirthFacade$get_missing()

------------------------------------------------------------------------

### Method `get_Plotter()`

get ContrastsPlotter

#### Usage

    ContrastsFirthFacade$get_Plotter(...)

#### Arguments

- `...`:

  passed to ContrastsFirth\$get_Plotter

------------------------------------------------------------------------

### Method `to_wide()`

convert results to wide format

#### Usage

    ContrastsFirthFacade$to_wide(...)

#### Arguments

- `...`:

  passed to ContrastsFirth\$to_wide

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsFirthFacade$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- sim_lfq_data_protein_config(Nprot = 20, weight_missing = 0.5)
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
lfqdata <- LFQData$new(istar$data, istar$config)
contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
fa <- ContrastsFirthFacade$new(lfqdata, "~ group_", contrasts)
#> completing cases
#> Joining with `by = join_by(protein_Id)`
head(fa$get_contrasts())
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct_firth
#> Joining with `by = join_by(protein_Id, contrast)`
#> # A tibble: 6 × 14
#> # Groups:   contrast [1]
#>   facade modelName     protein_Id contrast sigma    df      diff   FDR std.error
#>   <chr>  <chr>         <chr>      <chr>    <dbl> <int>     <dbl> <dbl>     <dbl>
#> 1 firth  WaldTestFirth 0EfVhX~59… A_vs_Ct…     1     9  1.07e-15 1          2.11
#> 2 firth  WaldTestFirth 0m5WN4~14… A_vs_Ct…     1     9  8.47e- 1 0.978      1.32
#> 3 firth  WaldTestFirth 7cbcrd~83… A_vs_Ct…     1     9  1.07e-15 1          2.11
#> 4 firth  WaldTestFirth 9VUkAq~45… A_vs_Ct…     1     9 -1.35e+ 0 0.978      1.78
#> 5 firth  WaldTestFirth At886V~32… A_vs_Ct…     1     9 -8.47e- 1 0.978      1.32
#> 6 firth  WaldTestFirth BEJI92~91… A_vs_Ct…     1     9 -1.35e+ 0 0.978      1.78
#> # ℹ 5 more variables: statistic <dbl>, p.value <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, avgAbd <dbl>
fa$to_wide()
#> # A tibble: 20 × 5
#>    protein_Id diff.A_vs_Ctrl p.value.A_vs_Ctrl FDR.A_vs_Ctrl statistic.A_vs_Ctrl
#>    <chr>               <dbl>             <dbl>         <dbl>               <dbl>
#>  1 0EfVhX~59…       1.07e-15             1             1                5.08e-16
#>  2 0m5WN4~14…       8.47e- 1             0.538         0.978            6.40e- 1
#>  3 7cbcrd~83…       1.07e-15             1             1                5.08e-16
#>  4 9VUkAq~45…      -1.35e+ 0             0.468         0.978           -7.58e- 1
#>  5 At886V~32…      -8.47e- 1             0.538         0.978           -6.40e- 1
#>  6 BEJI92~91…      -1.35e+ 0             0.468         0.978           -7.58e- 1
#>  7 CGzoYe~28…      -4.13e-16             1             1               -1.96e-16
#>  8 CtOJ9t~28…       1.35e+ 0             0.468         0.978            7.58e- 1
#>  9 DoWup2~29…       2.20e+ 0             0.238         0.978            1.26e+ 0
#> 10 DuwH7n~34…       8.47e- 1             0.538         0.978            6.40e- 1
#> 11 Fl4JiV~75…      -2.85e-17             1             1               -2.26e-17
#> 12 HC8K98~49…       8.47e- 1             0.538         0.978            6.40e- 1
#> 13 HvIpHG~40…       1.07e-15             1             1                5.08e-16
#> 14 I1Jk2Z~08…      -8.47e- 1             0.538         0.978           -6.40e- 1
#> 15 JV3Z7t~29…       1.07e-15             1             1                5.08e-16
#> 16 JcKVfU~08…      -1.35e+ 0             0.468         0.978           -7.58e- 1
#> 17 JfvT8X~27…      -2.20e+ 0             0.238         0.978           -1.26e+ 0
#> 18 R2i6w7~02…       6.65e-17             1             1                5.26e-17
#> 19 SGIVBl~95…       1.07e-15             1             1                5.08e-16
#> 20 r2J0Eh~26…      -4.13e-16             1             1               -1.96e-16
```
