# DEqMS count-dependent moderated contrasts

DEqMS count-dependent moderated contrasts

DEqMS count-dependent moderated contrasts

## Details

Decorator that wraps any Contrasts object and applies count-dependent
empirical Bayes variance shrinkage. Similar to
[`ContrastsModerated`](https://wolski.github.io/prolfqua/reference/ContrastsModerated.md)
but the prior variance depends on the number of quantified peptides/PSMs
per protein: proteins with many peptides get less shrinkage, proteins
with few peptides get more.

## See also

[`moderated_p_deqms_long`](https://wolski.github.io/prolfqua/reference/moderated_p_deqms_long.md)

Other modelling:
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsFirth`](https://wolski.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsFirthFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsLMMissingFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMMissingFacade.md),
[`ContrastsLimma`](https://wolski.github.io/prolfqua/reference/ContrastsLimma.md),
[`ContrastsLimmaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaFacade.md),
[`ContrastsLmerFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLmerFacade.md),
[`ContrastsMissing`](https://wolski.github.io/prolfqua/reference/ContrastsMissing.md),
[`ContrastsModerated`](https://wolski.github.io/prolfqua/reference/ContrastsModerated.md),
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

## Super class

[`prolfqua::ContrastsInterface`](https://wolski.github.io/prolfqua/reference/ContrastsInterface.md)
-\> `ContrastsModeratedDEqMS`

## Public fields

- `Contrast`:

  Class implementing the Contrast interface

- `count_df`:

  data.frame with subject_Id + count column

- `count_column`:

  name of the count column in count_df

- `loess_span`:

  span parameter for LOESS fit

- `modelName`:

  name of model

- `subject_Id`:

  columns with subject_Id (proteinID)

- `p.adjust`:

  function to adjust p-values

## Methods

### Public methods

- [`ContrastsModeratedDEqMS$new()`](#method-ContrastsModeratedDEqMS-new)

- [`ContrastsModeratedDEqMS$get_contrast_sides()`](#method-ContrastsModeratedDEqMS-get_contrast_sides)

- [`ContrastsModeratedDEqMS$get_linfct()`](#method-ContrastsModeratedDEqMS-get_linfct)

- [`ContrastsModeratedDEqMS$get_contrasts()`](#method-ContrastsModeratedDEqMS-get_contrasts)

- [`ContrastsModeratedDEqMS$get_Plotter()`](#method-ContrastsModeratedDEqMS-get_Plotter)

- [`ContrastsModeratedDEqMS$to_wide()`](#method-ContrastsModeratedDEqMS-to_wide)

- [`ContrastsModeratedDEqMS$clone()`](#method-ContrastsModeratedDEqMS-clone)

Inherited methods

- [`prolfqua::ContrastsInterface$column_description()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-column_description)

------------------------------------------------------------------------

### Method `new()`

initialize

#### Usage

    ContrastsModeratedDEqMS$new(
      Contrast,
      count_df,
      count_column,
      loess_span = 0.75,
      modelName = paste0(Contrast$modelName, "_DEqMS"),
      p.adjust = prolfqua::adjust_p_values
    )

#### Arguments

- `Contrast`:

  class implementing the ContrastInterface

- `count_df`:

  data.frame with subject_Id columns and a count column

- `count_column`:

  name of the count column in count_df

- `loess_span`:

  span for LOESS variance fit (default 0.75)

- `modelName`:

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

### Method `get_Plotter()`

get
[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.md)

#### Usage

    ContrastsModeratedDEqMS$get_Plotter(FCthreshold = 1, FDRthreshold = 0.1)

#### Arguments

- `FCthreshold`:

  fold change threshold to show in plots

- `FDRthreshold`:

  FDR threshold to show in plots

------------------------------------------------------------------------

### Method `to_wide()`

convert to wide format

#### Usage

    ContrastsModeratedDEqMS$to_wide(columns = c("p.value", "FDR", "statistic"))

#### Arguments

- `columns`:

  value column default p.value, FDR, statistic

#### Returns

data.frame

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
#> creating sampleName from fileName column
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
#> Joining with `by = join_by(protein_Id)`

Contr <- c("dil.b_vs_a" = "group_A - group_Ctrl")
contrast <- prolfqua::Contrasts$new(mod, Contr)

# Build count_df from config
count_df <- dplyr::select(protIntensity,
  dplyr::all_of(c(config$hierarchy_keys_depth(), "nr_peptides"))) |>
  dplyr::distinct()

deqms <- ContrastsModeratedDEqMS$new(contrast,
  count_df = count_df,
  count_column = "nr_peptides")

bb <- deqms$get_contrasts()
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`
#> Warning: moderated_p_deqms_long: warnings in 1/1 groups. contrast=dil.b_vs_a (pseudoinverse used at 1; neighborhood radius 1; reciprocal condition number  2.362e-17)
stopifnot(all(c("diff", "p.value", "FDR", "sigma") %in% colnames(bb)))

# Merge with ContrastsMissing
csi <- ContrastsMissing$new(lProt, contrasts = Contr)
merged <- merge_contrasts_results(deqms, csi)
#> Warning: moderated_p_deqms_long: warnings in 1/1 groups. contrast=dil.b_vs_a (pseudoinverse used at 1; neighborhood radius 1; reciprocal condition number  2.362e-17)
#> completing cases
#> dil.b_vs_a=group_A - group_Ctrl
#> dil.b_vs_a=group_A - group_Ctrl
#> dil.b_vs_a=group_A - group_Ctrl
#> Joining with `by = join_by(protein_Id, contrast)`
#> Joining with `by = join_by(protein_Id, contrast)`

cs <- deqms$get_contrast_sides()
cslf <- deqms$get_linfct()
ctrwide <- deqms$to_wide()
#> Warning: moderated_p_deqms_long: warnings in 1/1 groups. contrast=dil.b_vs_a (pseudoinverse used at 1; neighborhood radius 1; reciprocal condition number  2.362e-17)
cp <- deqms$get_Plotter()
#> Warning: moderated_p_deqms_long: warnings in 1/1 groups. contrast=dil.b_vs_a (pseudoinverse used at 1; neighborhood radius 1; reciprocal condition number  2.362e-17)
cp$volcano()
#> $FDR

#> 
```
