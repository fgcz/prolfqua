# ROPECA reproducibility-optimization method

ROPECA optimizes the reproducibility of statistical testing on
peptide-level and aggregates the peptide-level changes to determine
differential protein-level expression.

## Value

An R6 class generator.

## See also

Other modelling:
[`AnovaExtractor`](https://fgcz.github.io/prolfqua/reference/AnovaExtractor.md),
[`Contrasts`](https://fgcz.github.io/prolfqua/reference/Contrasts.md),
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

## Super class

[`ContrastsInterface`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.md)
-\> `ContrastsROPECA`

## Public fields

- `Contrast`:

  Contrast

- `contrast_result`:

  contrast result

- `model_name`:

  model name

- `subject_id`:

  columns with protein ID's

- `p.adjust`:

  method to use for p.value adjustment

## Methods

### Public methods

- [`ContrastsROPECA$new()`](#method-ContrastsROPECA-initialize)

- [`ContrastsROPECA$get_contrast_sides()`](#method-ContrastsROPECA-get_contrast_sides)

- [`ContrastsROPECA$get_linfct()`](#method-ContrastsROPECA-get_linfct)

- [`ContrastsROPECA$get_contrasts()`](#method-ContrastsROPECA-get_contrasts)

- [`ContrastsROPECA$get_Plotter()`](#method-ContrastsROPECA-get_Plotter)

- [`ContrastsROPECA$to_wide()`](#method-ContrastsROPECA-to_wide)

- [`ContrastsROPECA$clone()`](#method-ContrastsROPECA-clone)

Inherited methods

- [`ContrastsInterface$column_description()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-column_description)
- [`ContrastsInterface$contrast_summary_table()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-contrast_summary_table)
- [`ContrastsInterface$extra_artifacts()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-extra_artifacts)
- [`ContrastsInterface$filter_significant()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-filter_significant)
- [`ContrastsInterface$get_config()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_config)
- [`ContrastsInterface$get_missing()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_missing)
- [`ContrastsInterface$get_ora()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_ora)
- [`ContrastsInterface$get_rank()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_rank)

------------------------------------------------------------------------

### `ContrastsROPECA$new()`

initialize

#### Usage

    ContrastsROPECA$new(
      Contrast,
      model_name = "ROPECA",
      p.adjust = prolfqua::adjust_p_values
    )

#### Arguments

- `Contrast`:

  e.g. instance of Contrasts class, or ContrastsModerated

- `model_name`:

  default ROPECA

- `p.adjust`:

  function to use for p.value adjustment

------------------------------------------------------------------------

### `ContrastsROPECA$get_contrast_sides()`

show names of contrasts

#### Usage

    ContrastsROPECA$get_contrast_sides()

#### Returns

data.frame

------------------------------------------------------------------------

### `ContrastsROPECA$get_linfct()`

get linear function used to determine contrasts

#### Usage

    ContrastsROPECA$get_linfct()

#### Returns

data.frame

------------------------------------------------------------------------

### `ContrastsROPECA$get_contrasts()`

get contrasts

#### Usage

    ContrastsROPECA$get_contrasts(all = FALSE)

#### Arguments

- `all`:

  should all columns be returned (default FALSE)

- `global`:

  use a global linear function (determined by get_linfct)

#### Returns

data.frame

------------------------------------------------------------------------

### `ContrastsROPECA$get_Plotter()`

get
[`ContrastsPlotter`](https://fgcz.github.io/prolfqua/reference/ContrastsPlotter.md)

#### Usage

    ContrastsROPECA$get_Plotter(fdr_threshold = 0.1, fc_threshold = 2)

#### Arguments

- `fdr_threshold`:

  FDR threshold

- `fc_threshold`:

  FC threshold

#### Returns

[`ContrastsPlotter`](https://fgcz.github.io/prolfqua/reference/ContrastsPlotter.md)

------------------------------------------------------------------------

### `ContrastsROPECA$to_wide()`

convert to wide format

#### Usage

    ContrastsROPECA$to_wide(
      columns = c("beta.based.significance", "FDR.beta.based.significance")
    )

#### Arguments

- `columns`:

  value column default beta.based.significance

#### Returns

data.frame

------------------------------------------------------------------------

### `ContrastsROPECA$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsROPECA$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- prolfqua::sim_lfq_data_peptide_config(n_proteins=20)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 110 rows with NA response for feature and sample combinations absent from the input.
istar_data <- istar$data
modelFunction <-
  strategy_lm("abundance  ~ group_")
pepIntensity <- istar_data
config <- istar$config$clone(deep = TRUE)
config$hierarchy_depth <- 2
config$hierarchy_keys_depth()
#> [1] "protein_Id" "peptide_Id"

mod <- build_model(
 pepIntensity,
 modelFunction,
 subject_id = config$hierarchy_keys_depth())
#> fit (Model): 110 missing responses dropped (e.g. 0m5WN4~1448~UwJbtbRY, 0m5WN4~1448~wpHqGLYB, 7cbcrd~8305~O9LMQBev, 9VUkAq~4562~6mb5eKpF, 9VUkAq~4562~AYQOqCVO, ...): NA response (missing observation) left out of the fit
#> See $provenance().

 Contr <- c("AvsCtrl" = "group_A - group_Ctrl")


 contr <- prolfqua::Contrasts$new(mod, Contr)
 dim(contr$get_contrasts())
#> [1] 78 17
 contrM <- prolfqua::ContrastsModerated$new(contr)
 dim(contrM$get_contrasts())
#> [1] 78 17
 contrast <- prolfqua::ContrastsROPECA$new(contrM)
 contrast$get_contrasts()
#> # A tibble: 20 × 10
#>    modelName estimate_type protein_Id  contrast     n   diff statistic avgAbd
#>    <chr>     <chr>         <chr>       <chr>    <int>  <dbl>     <dbl>  <dbl>
#>  1 ROPECA    observed      0EfVhX~5954 AvsCtrl      1  9.23     13.8     25.1
#>  2 ROPECA    observed      0m5WN4~1448 AvsCtrl      2 -1.24     -1.66    18.5
#>  3 ROPECA    observed      7cbcrd~8305 AvsCtrl      1  4.65      6.43    25.7
#>  4 ROPECA    observed      9VUkAq~4562 AvsCtrl     16  0.980     1.38    18.1
#>  5 ROPECA    observed      At886V~3296 AvsCtrl      5 -2.67     -4.00    18.1
#>  6 ROPECA    observed      BEJI92~9143 AvsCtrl      4  1.66      2.49    24.6
#>  7 ROPECA    observed      CGzoYe~2857 AvsCtrl      1 -1.23     -1.43    17.6
#>  8 ROPECA    observed      CtOJ9t~2837 AvsCtrl      5  7.02     10.5     23.4
#>  9 ROPECA    observed      DoWup2~2934 AvsCtrl      7  0.831     1.24    21.2
#> 10 ROPECA    observed      DuwH7n~3402 AvsCtrl      3  0.655     0.980   19.6
#> 11 ROPECA    observed      Fl4JiV~7526 AvsCtrl      1  1.27      1.89    23.4
#> 12 ROPECA    observed      HC8K98~4958 AvsCtrl      2 -1.56     -1.37    15.3
#> 13 ROPECA    observed      HvIpHG~4015 AvsCtrl      2  1.99      2.72    16.4
#> 14 ROPECA    observed      I1Jk2Z~0821 AvsCtrl     10 -3.88     -4.13    15.5
#> 15 ROPECA    observed      JV3Z7t~2956 AvsCtrl      1 -5.06     -7.56    27.3
#> 16 ROPECA    observed      JcKVfU~0815 AvsCtrl      1 -1.51     -2.26    22.2
#> 17 ROPECA    observed      JfvT8X~2727 AvsCtrl     11  3.67      5.07    21.1
#> 18 ROPECA    observed      R2i6w7~0288 AvsCtrl      2  6.27      8.83    22.1
#> 19 ROPECA    observed      SGIVBl~9558 AvsCtrl      2  7.27     10.9     31.2
#> 20 ROPECA    observed      r2J0Eh~2687 AvsCtrl      1 -0.157    -0.235   22.3
#> # ℹ 2 more variables: beta.based.significance <dbl>,
#> #   FDR.beta.based.significance <dbl>
 contrast <- prolfqua::ContrastsROPECA$new(contr)
 tmp <- contrast$get_contrasts()
 dim(tmp)
#> [1] 20 10
 pl <- contrast$get_Plotter()
 contrast$to_wide()
#> # A tibble: 20 × 4
#>    protein_Id  diff.AvsCtrl beta.based.significance.Avs…¹ FDR.beta.based.signi…²
#>    <chr>              <dbl>                         <dbl>                  <dbl>
#>  1 0EfVhX~5954        9.23                      6.43 e- 7              2.57 e- 6
#>  2 0m5WN4~1448       -1.24                      9.64 e- 1              1.000e+ 0
#>  3 7cbcrd~8305        4.65                      1.15 e- 3              2.55 e- 3
#>  4 9VUkAq~4562        0.980                     1.01 e- 1              1.68 e- 1
#>  5 At886V~3296       -2.67                      1.11 e- 6              3.70 e- 6
#>  6 BEJI92~9143        1.66                      1.77 e- 3              3.54 e- 3
#>  7 CGzoYe~2857       -1.23                      2.05 e- 1              2.92 e- 1
#>  8 CtOJ9t~2837        7.02                      3.30 e-16              6.61 e-15
#>  9 DoWup2~2934        0.831                     1.75 e- 1              2.69 e- 1
#> 10 DuwH7n~3402        0.655                     2.25 e- 1              3.00 e- 1
#> 11 Fl4JiV~7526        1.27                      2.50 e- 1              3.12 e- 1
#> 12 HC8K98~4958       -1.56                      1.000e+ 0              1.000e+ 0
#> 13 HvIpHG~4015        1.99                      8.99 e- 1              9.99 e- 1
#> 14 I1Jk2Z~0821       -3.88                      2.32 e- 8              1.16 e- 7
#> 15 JV3Z7t~2956       -5.06                      1.18 e- 4              2.96 e- 4
#> 16 JcKVfU~0815       -1.51                      2.69 e- 2              4.90 e- 2
#> 17 JfvT8X~2727        3.67                      8.43 e-16              8.43 e-15
#> 18 R2i6w7~0288        6.27                      2.10 e- 5              6.01 e- 5
#> 19 SGIVBl~9558        7.27                      4.51 e- 9              3.01 e- 8
#> 20 r2J0Eh~2687       -0.157                     8.48 e- 1              9.98 e- 1
#> # ℹ abbreviated names: ¹​beta.based.significance.AvsCtrl,
#> #   ²​FDR.beta.based.significance.AvsCtrl
 contrast$get_linfct()
#>             (Intercept) group_B group_Ctrl
#> AvsCtrl               0       0       -1.0
#> avg_AvsCtrl           1       0        0.5
 contrast$get_contrast_sides()
#> # A tibble: 1 × 3
#>   contrast group_1 group_2   
#>   <chr>    <chr>   <chr>     
#> 1 AvsCtrl  group_A group_Ctrl
 pl$histogram()
#> $beta.based.significance

#> 
#> $FDR.beta.based.significance

#> 
 pl$ma_plot()

```
