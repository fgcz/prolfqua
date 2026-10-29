# DEqMS count-dependent moderated contrasts

DEqMS count-dependent moderated contrasts

DEqMS count-dependent moderated contrasts

## Value

An R6 class generator.

## Details

Decorator that wraps any Contrasts object and applies count-dependent
empirical Bayes variance shrinkage. Similar to
[`ContrastsModerated`](https://fgcz.github.io/prolfqua/reference/ContrastsModerated.md)
but the prior variance depends on the number of quantified peptides/PSMs
per protein: proteins with many peptides get less shrinkage, proteins
with few peptides get more.

## See also

[`moderated_p_deqms_long`](https://fgcz.github.io/prolfqua/reference/moderated_p_deqms_long.md)

Other modelling:
[`AnovaExtractor`](https://fgcz.github.io/prolfqua/reference/AnovaExtractor.md),
[`Contrasts`](https://fgcz.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsDEqMSVoomFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsDEqMSVoomFacade.md),
[`ContrastsFacadeBase`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.md),
[`ContrastsFirth`](https://fgcz.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsFirthFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsFirthNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsFirthNestedFacade.md),
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
[`ContrastsRfitFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsRfitFacade.md),
[`ContrastsRfitImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsRfitImputeFacade.md),
[`ContrastsTable`](https://fgcz.github.io/prolfqua/reference/ContrastsTable.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`LR_test()`](https://fgcz.github.io/prolfqua/reference/LR_test.md),
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
[`build_model_glm_peptide()`](https://fgcz.github.io/prolfqua/reference/build_model_glm_peptide.md),
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
[`merge_contrasts_results()`](https://fgcz.github.io/prolfqua/reference/merge_contrasts_results.md),
[`model_analyse()`](https://fgcz.github.io/prolfqua/reference/model_analyse.md),
[`model_summary()`](https://fgcz.github.io/prolfqua/reference/model_summary.md),
[`moderated_p_deqms()`](https://fgcz.github.io/prolfqua/reference/moderated_p_deqms.md),
[`moderated_p_deqms_long()`](https://fgcz.github.io/prolfqua/reference/moderated_p_deqms_long.md),
[`moderated_p_limma()`](https://fgcz.github.io/prolfqua/reference/moderated_p_limma.md),
[`new_imputed_model()`](https://fgcz.github.io/prolfqua/reference/new_imputed_model.md),
[`pivot_model_contrasts_to_wide()`](https://fgcz.github.io/prolfqua/reference/pivot_model_contrasts_to_wide.md),
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
[`summary_ROPECA_median_p_scaled()`](https://fgcz.github.io/prolfqua/reference/summary_ROPECA_median_p_scaled.md),
[`vcov.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/vcov.rfit_prolfqua.md)

## Super class

[`prolfqua::ContrastsInterface`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.md)
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

- [`ContrastsModeratedDEqMS$new()`](#method-ContrastsModeratedDEqMS-new)

- [`ContrastsModeratedDEqMS$get_contrast_sides()`](#method-ContrastsModeratedDEqMS-get_contrast_sides)

- [`ContrastsModeratedDEqMS$get_linfct()`](#method-ContrastsModeratedDEqMS-get_linfct)

- [`ContrastsModeratedDEqMS$get_contrasts()`](#method-ContrastsModeratedDEqMS-get_contrasts)

- [`ContrastsModeratedDEqMS$clone()`](#method-ContrastsModeratedDEqMS-clone)

Inherited methods

- [`prolfqua::ContrastsInterface$column_description()`](https://fgcz.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-column_description)
- [`prolfqua::ContrastsInterface$contrast_summary_table()`](https://fgcz.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-contrast_summary_table)
- [`prolfqua::ContrastsInterface$extra_artifacts()`](https://fgcz.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-extra_artifacts)
- [`prolfqua::ContrastsInterface$filter_significant()`](https://fgcz.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-filter_significant)
- [`prolfqua::ContrastsInterface$get_Plotter()`](https://fgcz.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_Plotter)
- [`prolfqua::ContrastsInterface$get_config()`](https://fgcz.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_config)
- [`prolfqua::ContrastsInterface$get_missing()`](https://fgcz.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_missing)
- [`prolfqua::ContrastsInterface$get_ora()`](https://fgcz.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_ora)
- [`prolfqua::ContrastsInterface$get_rank()`](https://fgcz.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_rank)
- [`prolfqua::ContrastsInterface$to_wide()`](https://fgcz.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-to_wide)

------------------------------------------------------------------------

### Method `new()`

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

### Method `get_contrast_sides()`

get both sides of contrasts

#### Usage

    ContrastsModeratedDEqMS$get_contrast_sides()

------------------------------------------------------------------------

### Method `get_linfct()`

get linear functions from contrasts

#### Usage

    ContrastsModeratedDEqMS$get_linfct(global = TRUE)

#### Arguments

- `global`:

  logical TRUE - get linear functions for all models

------------------------------------------------------------------------

### Method `get_contrasts()`

applies DEqMS-style count-dependent moderation

#### Usage

    ContrastsModeratedDEqMS$get_contrasts(all = FALSE)

#### Arguments

- `all`:

  should all columns be returned (default FALSE)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsModeratedDEqMS$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- sim_lfq_data_protein_config(Nprot = 50)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
protIntensity <- istar$data
config <- istar$config

lProt <- LFQData$new(protIntensity, config)
lProt$rename_response("transformedIntensity")
modelFunction <-
  strategy_lm("transformedIntensity ~ group_")
mod <- build_model(
 lProt,
 modelFunction)

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
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`
#> Warning: moderated_p_deqms_long: condition messages in 1/1 groups. contrast=dil.b_vs_a (pseudoinverse used at 1; neighborhood radius 1; reciprocal condition number  2.362e-17)
stopifnot(all(c("diff", "p.value", "FDR", "sigma") %in% colnames(bb)))

# Merge with LOD-imputed refits for proteins the plain lm cannot estimate
imp <- ContrastsLMImputeFacade$new(lProt, "~ group_", Contr)$contrast
#> Warning: essentially perfect fit: summary may be unreliable
merged <- merge_contrasts_results(deqms, imp)
#> Warning: moderated_p_deqms_long: condition messages in 1/1 groups. contrast=dil.b_vs_a (pseudoinverse used at 1; neighborhood radius 1; reciprocal condition number  2.362e-17)
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`
#> Joining with `by = join_by(protein_Id, contrast)`
#> Joining with `by = join_by(protein_Id, contrast)`

cs <- deqms$get_contrast_sides()
cslf <- deqms$get_linfct()
ctrwide <- deqms$to_wide()
#> Warning: moderated_p_deqms_long: condition messages in 1/1 groups. contrast=dil.b_vs_a (pseudoinverse used at 1; neighborhood radius 1; reciprocal condition number  2.362e-17)
cp <- deqms$get_Plotter()
#> Warning: moderated_p_deqms_long: condition messages in 1/1 groups. contrast=dil.b_vs_a (pseudoinverse used at 1; neighborhood radius 1; reciprocal condition number  2.362e-17)
cp$volcano()
#> $FDR

#> 
```
