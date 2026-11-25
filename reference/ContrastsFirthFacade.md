# Firth logistic missingness contrast analysis facade

Encapsulates the pipeline: encode missingness -\>
[`build_model_glm_protein`](https://fgcz.github.io/prolfqua/reference/build_model_glm_protein.md)
-\>
[`ContrastsFirth`](https://fgcz.github.io/prolfqua/reference/ContrastsFirth.md).

Requires aggregated (protein-level) LFQData. For nested peptide-level
input use
[`ContrastsBinomialFirthNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsBinomialFirthNestedFacade.md).

Contrasts are Wald tests on the log odds of detection with the normal
reference distribution (z, `df = Inf`), as in the coefficient table,
[`confint()`](https://rdrr.io/r/stats/confint.html) and
[`logistf`](https://rdrr.io/pkg/logistf/man/logistf.html) itself. With
few samples the test is conservative: in a null simulation with 10
samples per group and detection probability 0.5, 4.1% of the p-values
fell below 0.05.

## Value

An R6 class generator.

## See also

Other modelling:
[`AnovaExtractor`](https://fgcz.github.io/prolfqua/reference/AnovaExtractor.md),
[`Contrasts`](https://fgcz.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsComparison`](https://fgcz.github.io/prolfqua/reference/ContrastsComparison.md),
[`ContrastsDEqMSFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsDEqMSVoomFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsDEqMSVoomFacade.md),
[`ContrastsFacadeBase`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.md),
[`ContrastsFirth`](https://fgcz.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsHurdleFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsHurdleFacade.md),
[`ContrastsHurdleNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsHurdleNestedFacade.md),
[`ContrastsLMFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsLMImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLMImputeFacade.md),
[`ContrastsLimma`](https://fgcz.github.io/prolfqua/reference/ContrastsLimma.md),
[`ContrastsLimmaFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimmaFacade.md),
[`ContrastsLimmaImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimmaImputeFacade.md),
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
[`ContrastsROPECAV2`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAV2.md),
[`ContrastsROPECAV2NestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAV2NestedFacade.md),
[`ContrastsRfitFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsRfitFacade.md),
[`ContrastsRfitImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsRfitImputeFacade.md),
[`ContrastsTable`](https://fgcz.github.io/prolfqua/reference/ContrastsTable.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
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
[`build_model_glm_protein()`](https://fgcz.github.io/prolfqua/reference/build_model_glm_protein.md),
[`build_model_impute()`](https://fgcz.github.io/prolfqua/reference/build_model_impute.md),
[`build_model_limma()`](https://fgcz.github.io/prolfqua/reference/build_model_limma.md),
[`build_model_limma_impute()`](https://fgcz.github.io/prolfqua/reference/build_model_limma_impute.md),
[`build_model_limma_voom()`](https://fgcz.github.io/prolfqua/reference/build_model_limma_voom.md),
[`build_model_limma_voom_impute()`](https://fgcz.github.io/prolfqua/reference/build_model_limma_voom_impute.md),
[`build_model_limpa()`](https://fgcz.github.io/prolfqua/reference/build_model_limpa.md),
[`compute_borrowed_variance()`](https://fgcz.github.io/prolfqua/reference/compute_borrowed_variance.md),
[`compute_borrowed_variance_limma()`](https://fgcz.github.io/prolfqua/reference/compute_borrowed_variance_limma.md),
[`df.residual.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/df.residual.rfit_prolfqua.md),
[`group_label()`](https://fgcz.github.io/prolfqua/reference/group_label.md),
[`impute_from_model()`](https://fgcz.github.io/prolfqua/reference/impute_from_model.md),
[`impute_refit_singular()`](https://fgcz.github.io/prolfqua/reference/impute_refit_singular.md),
[`list_facades()`](https://fgcz.github.io/prolfqua/reference/list_facades.md),
[`lookup_facade()`](https://fgcz.github.io/prolfqua/reference/lookup_facade.md),
[`lr_test()`](https://fgcz.github.io/prolfqua/reference/lr_test.md),
[`merge_contrasts_results()`](https://fgcz.github.io/prolfqua/reference/merge_contrasts_results.md),
[`new_imputed_model()`](https://fgcz.github.io/prolfqua/reference/new_imputed_model.md),
[`nobs.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/nobs.rfit_prolfqua.md),
[`pivot_model_contrasts_to_wide()`](https://fgcz.github.io/prolfqua/reference/pivot_model_contrasts_to_wide.md),
[`prolfqua_analysis_methods`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md),
[`register_facade()`](https://fgcz.github.io/prolfqua/reference/register_facade.md),
[`sigma.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/sigma.rfit_prolfqua.md),
[`strategy_limma()`](https://fgcz.github.io/prolfqua/reference/strategy_limma.md),
[`strategy_limpa()`](https://fgcz.github.io/prolfqua/reference/strategy_limpa.md),
[`strategy_logistf()`](https://fgcz.github.io/prolfqua/reference/strategy.md),
[`vcov.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/vcov.rfit_prolfqua.md)

## Super classes

[`ContrastsInterface`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.md)
-\>
[`ContrastsFacadeBase`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.md)
-\> `ContrastsFirthFacade`

## Methods

### Public methods

- [`ContrastsFirthFacade$new()`](#method-ContrastsFirthFacade-initialize)

- [`ContrastsFirthFacade$clone()`](#method-ContrastsFirthFacade-clone)

Inherited methods

- [`ContrastsInterface$column_description()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-column_description)
- [`ContrastsInterface$contrast_summary_table()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-contrast_summary_table)
- [`ContrastsInterface$extra_artifacts()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-extra_artifacts)
- [`ContrastsInterface$filter_significant()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-filter_significant)
- [`ContrastsInterface$get_contrast_sides()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_contrast_sides)
- [`ContrastsInterface$get_ora()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_ora)
- [`ContrastsInterface$get_rank()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_rank)
- [`ContrastsFacadeBase$.recorded_provenance()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-.recorded_provenance)
- [`ContrastsFacadeBase$get_Plotter()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-get_Plotter)
- [`ContrastsFacadeBase$get_config()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-get_config)
- [`ContrastsFacadeBase$get_contrasts()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-get_contrasts)
- [`ContrastsFacadeBase$get_missing()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-get_missing)
- [`ContrastsFacadeBase$provenance()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-provenance)
- [`ContrastsFacadeBase$to_wide()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-to_wide)

------------------------------------------------------------------------

### `ContrastsFirthFacade$new()`

initialize

#### Usage

    ContrastsFirthFacade$new(lfqdata, modelstr, contrasts, control = list())

#### Arguments

- `lfqdata`:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

- `modelstr`:

  model formula string (e.g. "~ group\_")

- `contrasts`:

  named character vector of contrasts

- `control`:

  convergence controls of
  [`logistf`](https://rdrr.io/pkg/logistf/man/logistf.html), see
  [`strategy_logistf`](https://fgcz.github.io/prolfqua/reference/strategy.md)

------------------------------------------------------------------------

### `ContrastsFirthFacade$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsFirthFacade$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- sim_lfq_data_protein_config(n_proteins = 20, weight_missing = 0.5)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 5 unconfigured columns: average_prot_abundance, sd, mean, N, Replicate.
#> setup_analysis: added 51 rows with NA response for feature and sample combinations absent from the input.
lfqdata <- LFQData$new(istar$data, istar$config)
contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
fa <- ContrastsFirthFacade$new(lfqdata, "~ group_", contrasts)
head(fa$get_contrasts())
#> # A tibble: 6 × 16
#>   modelName estimate_type protein_Id  contrast  sigma    df      diff   FDR
#>   <chr>     <chr>         <chr>       <chr>     <dbl> <dbl>     <dbl> <dbl>
#> 1 firth     observed      0EfVhX~5954 A_vs_Ctrl     1   Inf  1.07e-15 1    
#> 2 firth     observed      0m5WN4~1448 A_vs_Ctrl     1   Inf  8.47e- 1 0.949
#> 3 firth     observed      7cbcrd~8305 A_vs_Ctrl     1   Inf  1.07e-15 1    
#> 4 firth     observed      9VUkAq~4562 A_vs_Ctrl     1   Inf -1.35e+ 0 0.949
#> 5 firth     observed      At886V~3296 A_vs_Ctrl     1   Inf -8.47e- 1 0.949
#> 6 firth     observed      BEJI92~9143 A_vs_Ctrl     1   Inf -1.35e+ 0 0.949
#> # ℹ 8 more variables: std.error <dbl>, statistic <dbl>, p.value <dbl>,
#> #   conf.low <dbl>, conf.high <dbl>, avgAbd <dbl>, std.error.unmoderated <dbl>,
#> #   df.unmoderated <dbl>
fa$to_wide()
#> # A tibble: 20 × 5
#>    protein_Id diff.A_vs_Ctrl p.value.A_vs_Ctrl FDR.A_vs_Ctrl statistic.A_vs_Ctrl
#>    <chr>               <dbl>             <dbl>         <dbl>               <dbl>
#>  1 0EfVhX~59…       1.07e-15             1             1                5.08e-16
#>  2 0m5WN4~14…       8.47e- 1             0.522         0.949            6.40e- 1
#>  3 7cbcrd~83…       1.07e-15             1             1                5.08e-16
#>  4 9VUkAq~45…      -1.35e+ 0             0.449         0.949           -7.58e- 1
#>  5 At886V~32…      -8.47e- 1             0.522         0.949           -6.40e- 1
#>  6 BEJI92~91…      -1.35e+ 0             0.449         0.949           -7.58e- 1
#>  7 CGzoYe~28…      -4.13e-16             1             1               -1.96e-16
#>  8 CtOJ9t~28…       1.35e+ 0             0.449         0.949            7.58e- 1
#>  9 DoWup2~29…       2.20e+ 0             0.206         0.949            1.26e+ 0
#> 10 DuwH7n~34…       8.47e- 1             0.522         0.949            6.40e- 1
#> 11 Fl4JiV~75…      -2.94e-17             1             1               -2.32e-17
#> 12 HC8K98~49…       8.47e- 1             0.522         0.949            6.40e- 1
#> 13 HvIpHG~40…       1.07e-15             1             1                5.08e-16
#> 14 I1Jk2Z~08…      -8.47e- 1             0.522         0.949           -6.40e- 1
#> 15 JV3Z7t~29…       1.07e-15             1             1                5.08e-16
#> 16 JcKVfU~08…      -1.35e+ 0             0.449         0.949           -7.58e- 1
#> 17 JfvT8X~27…      -2.20e+ 0             0.206         0.949           -1.26e+ 0
#> 18 R2i6w7~02…       4.22e-16             1             1                3.33e-16
#> 19 SGIVBl~95…       1.07e-15             1             1                5.08e-16
#> 20 r2J0Eh~26…      -4.13e-16             1             1               -1.96e-16
```
