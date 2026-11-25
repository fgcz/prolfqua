# ROPECA V2 on LOD-imputed peptide models, for nested input

Encapsulates the pipeline:
[`strategy_lm`](https://fgcz.github.io/prolfqua/reference/strategy.md)
-\>
[`build_model_impute`](https://fgcz.github.io/prolfqua/reference/build_model_impute.md)
(per peptide) -\>
[`Contrasts`](https://fgcz.github.io/prolfqua/reference/Contrasts.md)
-\>
[`ContrastsModerated`](https://fgcz.github.io/prolfqua/reference/ContrastsModerated.md)
-\>
[`ContrastsROPECAV2`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAV2.md).

Like
[`ContrastsROPECAV2NestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAV2NestedFacade.md)
(`ropecaV2_nested`), but the peptide models come from
[`build_model_impute`](https://fgcz.github.io/prolfqua/reference/build_model_impute.md),
as in `lm_impute`: a peptide with a group without observations is refit
at the limit of detection with a borrowed variance instead of being left
out. The correlation of the peptides of a protein is estimated from the
observed peptide models only, because refit peptides have near-constant
responses. Results that use a refit peptide have
`estimate_type = "lod_imputed"`. Input must be nested LFQData.

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
[`ContrastsFirthFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsFirthFacade.md),
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
[`ContrastsROPECAImputeNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAImputeNestedFacade.md),
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
-\>
[`ContrastsROPECAV2NestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAV2NestedFacade.md)
-\> `ContrastsROPECAV2ImputeNestedFacade`

## Methods

### Public methods

- [`ContrastsROPECAV2ImputeNestedFacade$new()`](#method-ContrastsROPECAV2ImputeNestedFacade-initialize)

- [`ContrastsROPECAV2ImputeNestedFacade$clone()`](#method-ContrastsROPECAV2ImputeNestedFacade-clone)

Inherited methods

- [`ContrastsInterface$column_description()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-column_description)
- [`ContrastsInterface$contrast_summary_table()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-contrast_summary_table)
- [`ContrastsInterface$extra_artifacts()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-extra_artifacts)
- [`ContrastsInterface$filter_significant()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-filter_significant)
- [`ContrastsInterface$get_contrast_sides()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_contrast_sides)
- [`ContrastsInterface$get_ora()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_ora)
- [`ContrastsInterface$get_rank()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_rank)
- [`ContrastsFacadeBase$.recorded_provenance()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-.recorded_provenance)
- [`ContrastsFacadeBase$get_config()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-get_config)
- [`ContrastsFacadeBase$get_missing()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-get_missing)
- [`ContrastsFacadeBase$provenance()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-provenance)
- [`ContrastsROPECAV2NestedFacade$get_Plotter()`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAV2NestedFacade.html#method-get_Plotter)
- [`ContrastsROPECAV2NestedFacade$get_contrasts()`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAV2NestedFacade.html#method-get_contrasts)
- [`ContrastsROPECAV2NestedFacade$to_wide()`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAV2NestedFacade.html#method-to_wide)

------------------------------------------------------------------------

### `ContrastsROPECAV2ImputeNestedFacade$new()`

initialize

#### Usage

    ContrastsROPECAV2ImputeNestedFacade$new(
      lfqdata,
      modelstr,
      contrasts,
      lod = NULL,
      borrow_method = c("sigma", "vcov"),
      df_method = c("observed", "borrowed")
    )

#### Arguments

- `lfqdata`:

  nested LFQData (peptide-level)

- `modelstr`:

  model formula string (e.g. "~ group\_")

- `contrasts`:

  named character vector of contrasts

- `lod`:

  limit of detection; NULL estimates it from the data (see
  [`build_model_impute`](https://fgcz.github.io/prolfqua/reference/build_model_impute.md))

- `borrow_method, df_method`:

  passed to
  [`build_model_impute`](https://fgcz.github.io/prolfqua/reference/build_model_impute.md)

------------------------------------------------------------------------

### `ContrastsROPECAV2ImputeNestedFacade$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsROPECAV2ImputeNestedFacade$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- sim_lfq_data_peptide_config(weight_missing = 0.3)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 52 rows with NA response for feature and sample combinations absent from the input.
lfqdata <- LFQData$new(istar$data, istar$config)
lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#> Transformed response column added: log2_abundance
contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
fa <- ContrastsROPECAV2ImputeNestedFacade$new(lfqdata, "~ group_", contrasts)
#> fit (ModelImputed): 52 missing responses dropped (e.g. 0EfVhX~0087~ITLb4x1q, 0EfVhX~0087~ahQLlQY7, 0EfVhX~0087~dJkdz7so, 7cbcrd~5725~D5dQ4nKk, 9VUkAq~4703~eIC06D7g, ...): NA response (missing observation) left out of the fit
#> See $provenance().
head(fa$get_contrasts())
#> # A tibble: 6 × 14
#>   modelName  estimate_type protein_Id contrast avgAbd    diff      FDR statistic
#>   <chr>      <chr>         <chr>      <chr>     <dbl>   <dbl>    <dbl>     <dbl>
#> 1 ropecaV2_… observed      0EfVhX~00… A_vs_Ct…   4.27 -0.0742 5.26e- 3     -2.82
#> 2 ropecaV2_… observed      7cbcrd~57… A_vs_Ct…   4.57  0.632  4.36e-10      6.49
#> 3 ropecaV2_… observed      9VUkAq~47… A_vs_Ct…   4.48 -0.573  8.70e-10     -6.32
#> 4 ropecaV2_… observed      BEJI92~52… A_vs_Ct…   4.25  0.245  3.48e- 5      4.26
#> 5 ropecaV2_… observed      CGzoYe~21… A_vs_Ct…   4.76 -0.310  8.39e- 8     -5.48
#> 6 ropecaV2_… observed      DoWup2~58… A_vs_Ct…   4.43  0.290  8.39e- 8      5.51
#> # ℹ 6 more variables: std.error <dbl>, df <dbl>, p.value <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, sigma <dbl>
fa$contrast$rho
#> [1] 0
```
