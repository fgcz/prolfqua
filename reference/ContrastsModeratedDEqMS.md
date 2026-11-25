# DEqMS count-dependent moderated contrasts

Decorator that wraps any Contrasts object and applies count-dependent
empirical Bayes variance shrinkage. Similar to
[`ContrastsModerated`](https://fgcz.github.io/prolfqua/reference/ContrastsModerated.md)
but the prior variance depends on the number of quantified peptides/PSMs
per protein: proteins with many peptides get less shrinkage, proteins
with few peptides get more.

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

## Super class

[`ContrastsInterface`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.md)
-\> `ContrastsModeratedDEqMS`

## Public fields

- `Contrast`:

  Class implementing the Contrast interface

- `count_df`:

  data.frame with subject_id + count column

- `count_column`:

  name of the count column in count_df

- `loess_span`:

  span parameter for LOESS fit

- `model_name`:

  name of model

- `subject_id`:

  columns with subject_id (proteinID)

- `p.adjust`:

  function to adjust p-values

## Methods

### Public methods

- [`ContrastsModeratedDEqMS$new()`](#method-ContrastsModeratedDEqMS-initialize)

- [`ContrastsModeratedDEqMS$get_contrast_sides()`](#method-ContrastsModeratedDEqMS-get_contrast_sides)

- [`ContrastsModeratedDEqMS$get_linfct()`](#method-ContrastsModeratedDEqMS-get_linfct)

- [`ContrastsModeratedDEqMS$get_contrasts()`](#method-ContrastsModeratedDEqMS-get_contrasts)

- [`ContrastsModeratedDEqMS$clone()`](#method-ContrastsModeratedDEqMS-clone)

Inherited methods

- [`ContrastsInterface$column_description()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-column_description)
- [`ContrastsInterface$contrast_summary_table()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-contrast_summary_table)
- [`ContrastsInterface$extra_artifacts()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-extra_artifacts)
- [`ContrastsInterface$filter_significant()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-filter_significant)
- [`ContrastsInterface$get_Plotter()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_Plotter)
- [`ContrastsInterface$get_config()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_config)
- [`ContrastsInterface$get_missing()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_missing)
- [`ContrastsInterface$get_ora()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_ora)
- [`ContrastsInterface$get_rank()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_rank)
- [`ContrastsInterface$to_wide()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-to_wide)

------------------------------------------------------------------------

### `ContrastsModeratedDEqMS$new()`

initialize

#### Usage

    ContrastsModeratedDEqMS$new(
      Contrast,
      count_df,
      count_column,
      loess_span = 0.75,
      model_name = paste0(Contrast$model_name, "_DEqMS"),
      p.adjust = prolfqua::adjust_p_values
    )

#### Arguments

- `Contrast`:

  class implementing the ContrastInterface

- `count_df`:

  data.frame with subject_id columns and a count column

- `count_column`:

  name of the count column in count_df

- `loess_span`:

  span for LOESS variance fit (default 0.75)

- `model_name`:

  name of the model

- `p.adjust`:

  function to adjust p-values - default BH

------------------------------------------------------------------------

### `ContrastsModeratedDEqMS$get_contrast_sides()`

get both sides of contrasts

#### Usage

    ContrastsModeratedDEqMS$get_contrast_sides()

------------------------------------------------------------------------

### `ContrastsModeratedDEqMS$get_linfct()`

get linear functions from contrasts

#### Usage

    ContrastsModeratedDEqMS$get_linfct(global = TRUE)

#### Arguments

- `global`:

  ignored; kept for compatibility with
  `Contrasts$get_linfct(global = )`. The linear functions come from the
  wrapped contrast object

------------------------------------------------------------------------

### `ContrastsModeratedDEqMS$get_contrasts()`

applies DEqMS-style count-dependent moderation

#### Usage

    ContrastsModeratedDEqMS$get_contrasts(all = FALSE)

#### Arguments

- `all`:

  should all columns be returned (default FALSE)

------------------------------------------------------------------------

### `ContrastsModeratedDEqMS$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsModeratedDEqMS$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- sim_lfq_data_protein_config(n_proteins = 50)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 5 unconfigured columns: average_prot_abundance, sd, mean, N, Replicate.
#> setup_analysis: added 82 rows with NA response for feature and sample combinations absent from the input.
protIntensity <- istar$data
config <- istar$config

lProt <- LFQData$new(protIntensity, config)
lProt$rename_response("transformedIntensity")
modelFunction <-
  strategy_lm("transformedIntensity ~ group_")
mod <- build_model(
 lProt,
 modelFunction)
#> fit (Model): 1 singular fits (e.g. VUnZ6o~1614): coefficients not estimable or fewer than 2 residual df
#> fit (Model): 82 missing responses dropped (e.g. 0EfVhX~7161, 0m5WN4~3543, 7QuTub~5556, 7cbcrd~0495, 9VUkAq~8655, ...): NA response (missing observation) left out of the fit
#> See $provenance().

Contr <- c("dil.b_vs_a" = "group_A - group_Ctrl")
contrast <- prolfqua::Contrasts$new(mod, Contr)

# Build count_df from config
count_df <- dplyr::select(protIntensity,
  dplyr::all_of(c(config$hierarchy_keys_depth(), "nrPeptides"))) |>
  dplyr::distinct()

deqms <- ContrastsModeratedDEqMS$new(contrast,
  count_df = count_df,
  count_column = "nrPeptides")

bb <- deqms$get_contrasts()
#> Warning: moderated_p_deqms_long: condition messages in 1/1 groups. contrast=dil.b_vs_a (pseudoinverse used at 1; neighborhood radius 1; reciprocal condition number  7.3328e-17)
stopifnot(all(c("diff", "p.value", "FDR", "sigma") %in% colnames(bb)))

# Merge with LOD-imputed refits for proteins the plain lm cannot estimate
imp <- ContrastsLMImputeFacade$new(lProt, "~ group_", Contr)$contrast
#> fit (ModelImputed): 1 refit at the limit of detection (e.g. VUnZ6o~1614): fit failed, singular or missing a group; missing responses imputed at the LOD
#> fit (ModelImputed): 73 missing responses dropped (e.g. 0EfVhX~7161, 0m5WN4~3543, 7QuTub~5556, 7cbcrd~0495, 9VUkAq~8655, ...): NA response (missing observation) left out of the fit
#> See $provenance().
merged <- merge_contrasts_results(deqms, imp)
#> Warning: moderated_p_deqms_long: condition messages in 1/1 groups. contrast=dil.b_vs_a (pseudoinverse used at 1; neighborhood radius 1; reciprocal condition number  7.3328e-17)

cs <- deqms$get_contrast_sides()
cslf <- deqms$get_linfct()
ctrwide <- deqms$to_wide()
#> Warning: moderated_p_deqms_long: condition messages in 1/1 groups. contrast=dil.b_vs_a (pseudoinverse used at 1; neighborhood radius 1; reciprocal condition number  7.3328e-17)
cp <- deqms$get_Plotter()
#> Warning: moderated_p_deqms_long: condition messages in 1/1 groups. contrast=dil.b_vs_a (pseudoinverse used at 1; neighborhood radius 1; reciprocal condition number  7.3328e-17)
cp$volcano()
#> $FDR

#> 
```
