# ROPECA V2: median peptide z-score with a correlation-aware null

Summarises peptide-level contrasts per protein like
[`ContrastsROPECA`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECA.md),
but tests the median peptide z-score against its exact null
distribution: the median of n standard normal scores with pairwise
correlation `rho`. The protein p-values are uniform under the null for
any number of peptides, and correlated peptides do not make proteins
with many peptides more significant.
[`ContrastsROPECAV2NestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAV2NestedFacade.md)
estimates `rho` from the peptide model residuals.

Each peptide t-statistic is converted to the standard normal score with
the same tail probability. `diff` is the median peptide fold change,
`statistic` the normal score with the protein's two-sided p-value,
`std.error` = `|diff / statistic|`, `df` = `Inf`, and `conf.low`,
`conf.high` the matching normal 95\\ `sigma` is the MAD of the peptide
fold changes.

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
[`ContrastsROPECA`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECA.md),
[`ContrastsROPECANestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECANestedFacade.md),
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
[`build_model_logistf()`](https://fgcz.github.io/prolfqua/reference/build_model_logistf.md),
[`compute_borrowed_variance()`](https://fgcz.github.io/prolfqua/reference/compute_borrowed_variance.md),
[`compute_borrowed_variance_limma()`](https://fgcz.github.io/prolfqua/reference/compute_borrowed_variance_limma.md),
[`compute_contrast()`](https://fgcz.github.io/prolfqua/reference/compute_contrast.md),
[`compute_lmer_contrast()`](https://fgcz.github.io/prolfqua/reference/compute_lmer_contrast.md),
[`df.residual.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/df.residual.rfit_prolfqua.md),
[`get_anova_df()`](https://fgcz.github.io/prolfqua/reference/get_anova_df.md),
[`get_complete_model_fit()`](https://fgcz.github.io/prolfqua/reference/get_complete_model_fit.md),
[`get_p_values_pbeta()`](https://fgcz.github.io/prolfqua/reference/get_p_values_pbeta.md),
[`group_label()`](https://fgcz.github.io/prolfqua/reference/group_label.md),
[`impute_from_model()`](https://fgcz.github.io/prolfqua/reference/impute_from_model.md),
[`impute_refit_singular()`](https://fgcz.github.io/prolfqua/reference/impute_refit_singular.md),
[`linfct_from_model()`](https://fgcz.github.io/prolfqua/reference/linfct_from_model.md),
[`linfct_matrix_contrasts()`](https://fgcz.github.io/prolfqua/reference/linfct_matrix_contrasts.md),
[`list_facades()`](https://fgcz.github.io/prolfqua/reference/list_facades.md),
[`lookup_facade()`](https://fgcz.github.io/prolfqua/reference/lookup_facade.md),
[`lr_test()`](https://fgcz.github.io/prolfqua/reference/lr_test.md),
[`merge_contrasts_results()`](https://fgcz.github.io/prolfqua/reference/merge_contrasts_results.md),
[`model_analyse()`](https://fgcz.github.io/prolfqua/reference/model_analyse.md),
[`model_summary()`](https://fgcz.github.io/prolfqua/reference/model_summary.md),
[`moderated_p_deqms()`](https://fgcz.github.io/prolfqua/reference/moderated_p_deqms.md),
[`moderated_p_deqms_long()`](https://fgcz.github.io/prolfqua/reference/moderated_p_deqms_long.md),
[`moderated_p_limma()`](https://fgcz.github.io/prolfqua/reference/moderated_p_limma.md),
[`new_imputed_model()`](https://fgcz.github.io/prolfqua/reference/new_imputed_model.md),
[`nobs.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/nobs.rfit_prolfqua.md),
[`pivot_model_contrasts_to_wide()`](https://fgcz.github.io/prolfqua/reference/pivot_model_contrasts_to_wide.md),
[`prolfqua_analysis_methods`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md),
[`register_facade()`](https://fgcz.github.io/prolfqua/reference/register_facade.md),
[`sigma.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/sigma.rfit_prolfqua.md),
[`sim_build_models_lm()`](https://fgcz.github.io/prolfqua/reference/sim_build_models_lm.md),
[`sim_build_models_lmer()`](https://fgcz.github.io/prolfqua/reference/sim_build_models_lmer.md),
[`sim_build_models_logistf()`](https://fgcz.github.io/prolfqua/reference/sim_build_models_logistf.md),
[`sim_make_model_lm()`](https://fgcz.github.io/prolfqua/reference/sim_make_model_lm.md),
[`sim_make_model_lmer()`](https://fgcz.github.io/prolfqua/reference/sim_make_model_lmer.md),
[`strategy_limma()`](https://fgcz.github.io/prolfqua/reference/strategy_limma.md),
[`strategy_limpa()`](https://fgcz.github.io/prolfqua/reference/strategy_limpa.md),
[`strategy_logistf()`](https://fgcz.github.io/prolfqua/reference/strategy.md),
[`summary_ropeca_median_p_scaled()`](https://fgcz.github.io/prolfqua/reference/summary_ropeca_median_p_scaled.md),
[`vcov.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/vcov.rfit_prolfqua.md)

## Super class

[`ContrastsInterface`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.md)
-\> `ContrastsROPECAV2`

## Public fields

- `Contrast`:

  peptide-level contrasts, e.g. `ContrastsModerated`

- `rho`:

  correlation of the peptides of a protein

- `contrast_result`:

  contrast result

- `model_name`:

  model name

- `subject_id`:

  peptide-level subject columns; the protein columns are all but the
  last

- `p.adjust`:

  method to use for p.value adjustment

## Methods

### Public methods

- [`ContrastsROPECAV2$new()`](#method-ContrastsROPECAV2-initialize)

- [`ContrastsROPECAV2$get_contrast_sides()`](#method-ContrastsROPECAV2-get_contrast_sides)

- [`ContrastsROPECAV2$get_linfct()`](#method-ContrastsROPECAV2-get_linfct)

- [`ContrastsROPECAV2$get_contrasts()`](#method-ContrastsROPECAV2-get_contrasts)

- [`ContrastsROPECAV2$get_Plotter()`](#method-ContrastsROPECAV2-get_Plotter)

- [`ContrastsROPECAV2$to_wide()`](#method-ContrastsROPECAV2-to_wide)

- [`ContrastsROPECAV2$clone()`](#method-ContrastsROPECAV2-clone)

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

### `ContrastsROPECAV2$new()`

initialize

#### Usage

    ContrastsROPECAV2$new(
      Contrast,
      rho,
      model_name = "ROPECAV2",
      p.adjust = prolfqua::adjust_p_values
    )

#### Arguments

- `Contrast`:

  peptide-level contrasts with columns `diff`, `statistic` and `df`

- `rho`:

  correlation of the peptides of a protein, in \[0, 1)

- `model_name`:

  default ROPECAV2

- `p.adjust`:

  function to use for p.value adjustment

------------------------------------------------------------------------

### `ContrastsROPECAV2$get_contrast_sides()`

show names of contrasts

#### Usage

    ContrastsROPECAV2$get_contrast_sides()

#### Returns

data.frame

------------------------------------------------------------------------

### `ContrastsROPECAV2$get_linfct()`

get linear function used to determine contrasts

#### Usage

    ContrastsROPECAV2$get_linfct()

#### Returns

data.frame

------------------------------------------------------------------------

### `ContrastsROPECAV2$get_contrasts()`

get protein-level contrasts

#### Usage

    ContrastsROPECAV2$get_contrasts(all = FALSE)

#### Arguments

- `all`:

  `FALSE` (default) returns the standard columns; `TRUE` appends
  `n_peptides`, `median_z` and `rho`

#### Returns

data.frame

------------------------------------------------------------------------

### `ContrastsROPECAV2$get_Plotter()`

get
[`ContrastsPlotter`](https://fgcz.github.io/prolfqua/reference/ContrastsPlotter.md)

#### Usage

    ContrastsROPECAV2$get_Plotter(fdr_threshold = 0.1, fc_threshold = 2)

#### Arguments

- `fdr_threshold`:

  FDR threshold

- `fc_threshold`:

  FC threshold

#### Returns

[`ContrastsPlotter`](https://fgcz.github.io/prolfqua/reference/ContrastsPlotter.md)

------------------------------------------------------------------------

### `ContrastsROPECAV2$to_wide()`

convert to wide format

#### Usage

    ContrastsROPECAV2$to_wide(columns = c("p.value", "FDR"))

#### Arguments

- `columns`:

  value columns, default `c("p.value", "FDR")`

#### Returns

data.frame

------------------------------------------------------------------------

### `ContrastsROPECAV2$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsROPECAV2$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- sim_lfq_data_peptide_config(n_proteins = 20)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 110 rows with NA response for feature and sample combinations absent from the input.
config <- istar$config$clone(deep = TRUE)
config$hierarchy_depth <- 2
mod <- build_model(
  istar$data,
  strategy_lm("abundance ~ group_"),
  subject_id = config$hierarchy_keys_depth()
)
#> fit (Model): 110 missing responses dropped (e.g. 0m5WN4~1448~UwJbtbRY, 0m5WN4~1448~wpHqGLYB, 7cbcrd~8305~O9LMQBev, 9VUkAq~4562~6mb5eKpF, 9VUkAq~4562~AYQOqCVO, ...): NA response (missing observation) left out of the fit
#> See $provenance().
contrasts <- c("AvsCtrl" = "group_A - group_Ctrl")
peptide_contrasts <- ContrastsModerated$new(Contrasts$new(mod, contrasts))
ropeca2 <- ContrastsROPECAV2$new(peptide_contrasts, rho = 0.2)
head(ropeca2$get_contrasts())
#> # A tibble: 6 × 14
#>   protein_Id  modelName estimate_type contrast avgAbd   diff      FDR statistic
#>   <chr>       <chr>     <chr>         <chr>     <dbl>  <dbl>    <dbl>     <dbl>
#> 1 0EfVhX~5954 ROPECAV2  observed      AvsCtrl    25.1  9.23  8.17e-42     13.7 
#> 2 0m5WN4~1448 ROPECAV2  observed      AvsCtrl    18.5 -1.24  4.56e- 2     -2.15
#> 3 7cbcrd~8305 ROPECAV2  observed      AvsCtrl    25.7  4.65  3.28e-10      6.43
#> 4 9VUkAq~4562 ROPECAV2  observed      AvsCtrl    18.1  0.980 1.35e- 2      2.65
#> 5 At886V~3296 ROPECAV2  observed      AvsCtrl    18.1 -2.67  3.66e- 9     -6.03
#> 6 BEJI92~9143 ROPECAV2  observed      AvsCtrl    24.6  1.66  3.62e- 4      3.74
#> # ℹ 6 more variables: std.error <dbl>, df <dbl>, p.value <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, sigma <dbl>
ropeca2$to_wide()
#> # A tibble: 20 × 4
#>    protein_Id  diff.AvsCtrl p.value.AvsCtrl FDR.AvsCtrl
#>    <chr>              <dbl>           <dbl>       <dbl>
#>  1 0EfVhX~5954        9.23         8.17e-43    8.17e-42
#>  2 0m5WN4~1448       -1.24         3.19e- 2    4.56e- 2
#>  3 7cbcrd~8305        4.65         1.31e-10    3.28e-10
#>  4 9VUkAq~4562        0.980        8.12e- 3    1.35e- 2
#>  5 At886V~3296       -2.67         1.65e- 9    3.66e- 9
#>  6 BEJI92~9143        1.66         1.81e- 4    3.62e- 4
#>  7 CGzoYe~2857       -1.23         1.53e- 1    1.70e- 1
#>  8 CtOJ9t~2837        7.02         3.52e-54    7.03e-53
#>  9 DoWup2~2934        0.831        4.08e- 2    5.44e- 2
#> 10 DuwH7n~3402        0.655        1.90e- 1    2.00e- 1
#> 11 Fl4JiV~7526        1.27         5.85e- 2    7.31e- 2
#> 12 HC8K98~4958       -1.56         7.72e- 2    9.08e- 2
#> 13 HvIpHG~4015        1.99         4.56e- 4    8.29e- 4
#> 14 I1Jk2Z~0821       -3.88         1.75e-13    5.01e-13
#> 15 JV3Z7t~2956       -5.06         4.48e-14    1.49e-13
#> 16 JcKVfU~0815       -1.51         2.37e- 2    3.65e- 2
#> 17 JfvT8X~2727        3.67         1.54e-19    6.14e-19
#> 18 R2i6w7~0288        6.27         1.85e-29    9.27e-29
#> 19 SGIVBl~9558        7.27         1.12e-29    7.48e-29
#> 20 r2J0Eh~2687       -0.157        8.15e- 1    8.15e- 1
pl <- ropeca2$get_Plotter()
pl$volcano()
#> $FDR

#> 
```
