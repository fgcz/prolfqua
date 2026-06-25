# Estimate contrasts using Wald Test

Estimate contrasts using Wald Test

Estimate contrasts using Wald Test

## Value

An R6 class generator.

## Details

The per-protein contrast computation uses
[`compute_contrast`](https://wolski.github.io/prolfqua/reference/compute_contrast.md)
and
[`linfct_matrix_contrasts`](https://wolski.github.io/prolfqua/reference/linfct_matrix_contrasts.md)
internally. Both support a vectorized code path activated by
`options(prolfqua.vectorize = TRUE)`. This can give a significant
speed-up for large datasets. The vectorized path produces numerically
identical results.

## See also

Other modelling:
[`AnovaExtractor`](https://wolski.github.io/prolfqua/reference/AnovaExtractor.md),
[`ContrastsDEqMSFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsDEqMSVoomFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSVoomFacade.md),
[`ContrastsFacadeBase`](https://wolski.github.io/prolfqua/reference/ContrastsFacadeBase.md),
[`ContrastsFirth`](https://wolski.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsFirthFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsFirthNestedFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthNestedFacade.md),
[`ContrastsLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsLMImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMImputeFacade.md),
[`ContrastsLMMissingFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMMissingFacade.md),
[`ContrastsLimma`](https://wolski.github.io/prolfqua/reference/ContrastsLimma.md),
[`ContrastsLimmaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaFacade.md),
[`ContrastsLimmaImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaImputeFacade.md),
[`ContrastsLimmaVoomFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaVoomFacade.md),
[`ContrastsLimmaVoomImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaVoomImputeFacade.md),
[`ContrastsLimpaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimpaFacade.md),
[`ContrastsLimpaNestedFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimpaNestedFacade.md),
[`ContrastsLmerNestedFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLmerNestedFacade.md),
[`ContrastsMissing`](https://wolski.github.io/prolfqua/reference/ContrastsMissing.md),
[`ContrastsModerated`](https://wolski.github.io/prolfqua/reference/ContrastsModerated.md),
[`ContrastsModeratedDEqMS`](https://wolski.github.io/prolfqua/reference/ContrastsModeratedDEqMS.md),
[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.md),
[`ContrastsRLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsRLMFacade.md),
[`ContrastsROPECA`](https://wolski.github.io/prolfqua/reference/ContrastsROPECA.md),
[`ContrastsROPECANestedFacade`](https://wolski.github.io/prolfqua/reference/ContrastsROPECANestedFacade.md),
[`ContrastsRfitFacade`](https://wolski.github.io/prolfqua/reference/ContrastsRfitFacade.md),
[`ContrastsRfitImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsRfitImputeFacade.md),
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
[`StrategyRfit`](https://wolski.github.io/prolfqua/reference/StrategyRfit.md),
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
[`df.residual.rfit_prolfqua()`](https://wolski.github.io/prolfqua/reference/df.residual.rfit_prolfqua.md),
[`get_anova_df()`](https://wolski.github.io/prolfqua/reference/get_anova_df.md),
[`get_complete_model_fit()`](https://wolski.github.io/prolfqua/reference/get_complete_model_fit.md),
[`get_p_values_pbeta()`](https://wolski.github.io/prolfqua/reference/get_p_values_pbeta.md),
[`group_label()`](https://wolski.github.io/prolfqua/reference/group_label.md),
[`impute_refit_singular()`](https://wolski.github.io/prolfqua/reference/impute_refit_singular.md),
[`is_singular_lm()`](https://wolski.github.io/prolfqua/reference/is_singular_lm.md),
[`linfct_all_possible_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_all_possible_contrasts.md),
[`linfct_factors_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_factors_contrasts.md),
[`linfct_from_model()`](https://wolski.github.io/prolfqua/reference/linfct_from_model.md),
[`linfct_matrix_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_matrix_contrasts.md),
[`list_facades()`](https://wolski.github.io/prolfqua/reference/list_facades.md),
[`lookup_facade()`](https://wolski.github.io/prolfqua/reference/lookup_facade.md),
[`merge_contrasts_results()`](https://wolski.github.io/prolfqua/reference/merge_contrasts_results.md),
[`model_analyse()`](https://wolski.github.io/prolfqua/reference/model_analyse.md),
[`model_summary()`](https://wolski.github.io/prolfqua/reference/model_summary.md),
[`moderated_p_deqms()`](https://wolski.github.io/prolfqua/reference/moderated_p_deqms.md),
[`moderated_p_deqms_long()`](https://wolski.github.io/prolfqua/reference/moderated_p_deqms_long.md),
[`moderated_p_limma()`](https://wolski.github.io/prolfqua/reference/moderated_p_limma.md),
[`moderated_p_limma_long()`](https://wolski.github.io/prolfqua/reference/moderated_p_limma_long.md),
[`new_imputed_model()`](https://wolski.github.io/prolfqua/reference/new_imputed_model.md),
[`pivot_model_contrasts_to_wide()`](https://wolski.github.io/prolfqua/reference/pivot_model_contrasts_to_wide.md),
[`plot_lmer_peptide_predictions()`](https://wolski.github.io/prolfqua/reference/plot_lmer_peptide_predictions.md),
[`register_facade()`](https://wolski.github.io/prolfqua/reference/register_facade.md),
[`sigma.rfit_prolfqua()`](https://wolski.github.io/prolfqua/reference/sigma.rfit_prolfqua.md),
[`sim_build_models_lm()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lm.md),
[`sim_build_models_lmer()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lmer.md),
[`sim_build_models_logistf()`](https://wolski.github.io/prolfqua/reference/sim_build_models_logistf.md),
[`sim_make_model_lm()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lm.md),
[`sim_make_model_lmer()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lmer.md),
[`strategy_limma()`](https://wolski.github.io/prolfqua/reference/strategy_limma.md),
[`strategy_limpa()`](https://wolski.github.io/prolfqua/reference/strategy_limpa.md),
[`strategy_logistf()`](https://wolski.github.io/prolfqua/reference/strategy.md),
[`summary_ROPECA_median_p.scaled()`](https://wolski.github.io/prolfqua/reference/summary_ROPECA_median_p.scaled.md),
[`unregister_facade()`](https://wolski.github.io/prolfqua/reference/unregister_facade.md),
[`vcov.rfit_prolfqua()`](https://wolski.github.io/prolfqua/reference/vcov.rfit_prolfqua.md)

## Super class

[`prolfqua::ContrastsInterface`](https://wolski.github.io/prolfqua/reference/ContrastsInterface.md)
-\> `Contrast`

## Public fields

- `models`:

  Model

- `contrasts`:

  character with contrasts

- `contrastfun`:

  function to compute contrasts

- `model_name`:

  model name

- `subject_id`:

  name of column containing e.g., protein Id's

- `p.adjust`:

  function to adjust p-values (default prolfqua::adjust_p_values)

- `contrast_result`:

  data frame containing results of contrast computation

- `global`:

  use a global linear function (determined by get_linfct)

- `protein_annot`:

  holds protein annotation

## Methods

### Public methods

- [`Contrasts$new()`](#method-Contrast-new)

- [`Contrasts$get_contrast_sides()`](#method-Contrast-get_contrast_sides)

- [`Contrasts$get_linfct()`](#method-Contrast-get_linfct)

- [`Contrasts$get_contrasts()`](#method-Contrast-get_contrasts)

- [`Contrasts$get_Plotter()`](#method-Contrast-get_Plotter)

- [`Contrasts$to_wide()`](#method-Contrast-to_wide)

- [`Contrasts$clone()`](#method-Contrast-clone)

Inherited methods

- [`prolfqua::ContrastsInterface$column_description()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-column_description)
- [`prolfqua::ContrastsInterface$contrast_summary_table()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-contrast_summary_table)
- [`prolfqua::ContrastsInterface$extra_artifacts()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-extra_artifacts)
- [`prolfqua::ContrastsInterface$filter_significant()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-filter_significant)
- [`prolfqua::ContrastsInterface$get_config()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_config)
- [`prolfqua::ContrastsInterface$get_missing()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_missing)
- [`prolfqua::ContrastsInterface$get_ora()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_ora)
- [`prolfqua::ContrastsInterface$get_rank()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_rank)

------------------------------------------------------------------------

### Method `new()`

initialize create Contrast

#### Usage

    Contrasts$new(
      model,
      contrasts,
      p.adjust = prolfqua::adjust_p_values,
      global = FALSE,
      model_name = "WaldTest"
    )

#### Arguments

- `model`:

  a dataframe with a structure similar to that generated by
  [`build_model`](https://wolski.github.io/prolfqua/reference/build_model.md)

- `contrasts`:

  a character vector with contrast specificiation

- `p.adjust`:

  function to adjust the p-values

- `global`:

  development/internal argument (if FALSE determine linfct for each
  model.)

- `model_name`:

  name of contrast method, default WaldTest

------------------------------------------------------------------------

### Method `get_contrast_sides()`

get both sides of contrasts

#### Usage

    Contrasts$get_contrast_sides()

------------------------------------------------------------------------

### Method `get_linfct()`

get linear functions from contrasts

#### Usage

    Contrasts$get_linfct(global = TRUE, avg = TRUE)

#### Arguments

- `global`:

  logical TRUE - get the a linear functions for all models, FALSE -
  linear function for each model

- `avg`:

  logical TRUE - get also linfct for averages

------------------------------------------------------------------------

### Method `get_contrasts()`

get table with contrast estimates

#### Usage

    Contrasts$get_contrasts(all = FALSE)

#### Arguments

- `all`:

  should all columns be returned (default FALSE)

#### Returns

data.frame with contrasts

------------------------------------------------------------------------

### Method `get_Plotter()`

return
[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.md)
creates Contrast_Plotter

#### Usage

    Contrasts$get_Plotter(fc_threshold = 1, fdr_threshold = 0.1)

#### Arguments

- `fc_threshold`:

  fold change threshold to show in plots

- `fdr_threshold`:

  FDR threshold to show in plots

#### Returns

[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.md)

------------------------------------------------------------------------

### Method `to_wide()`

convert to wide format

#### Usage

    Contrasts$to_wide(columns = c("p.value", "FDR", "statistic"))

#### Arguments

- `columns`:

  value column default p.value

#### Returns

data.frame

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Contrasts$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# Fitting mixed effects model to peptide data
istar <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done

modelFunction <-
  strategy_lmer("abundance  ~ group_ + (1 | peptide_Id) + (1 | sample)")

config <- istar$config
config$hierarchy_keys_depth()
#> [1] "protein_Id"

mod <- build_model(
  istar$data,
  modelFunction,
  subject_id = config$hierarchy_keys_depth()
)
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> Warning: There were 4 warnings in `dplyr::mutate()`.
#> The first warning was:
#> ℹ In argument: `linear_model = purrr::map(data, model_strategy$model_fun, pb =
#>   pb)`.
#> ℹ In group 2: `protein_Id = "7cbcrd~5725"`.
#> Caused by warning:
#> ! grouping factors must have > 1 sampled level
#> ℹ Run `dplyr::last_dplyr_warnings()` to see the 3 remaining warnings.

ref_lfc <- data.frame(
  `(Intercept)` = c(0, 0, 0),
  group_B = c(0, 1, 0),
  group_Ctrl = c(-1, -1, 1),
  row.names = c("groupA_vs_Ctrl", "dil.e_vs_b", "dil.ctrl_vs_b")
)
prolfqua::model_summary(mod)
#> $exists
#> 
#> FALSE  TRUE 
#>     4     6 
#> 
#> $isSingular
#> 
#> FALSE  TRUE 
#>     1     5 
#> 
Contr <- c(
  "groupA_vs_Ctrl" = "group_A - group_Ctrl",
  "dil.e_vs_b" = "group_B - group_Ctrl",
  "dil.ctrl_vs_b" = "group_Ctrl - group_A"
)
contrastX <- prolfqua::Contrasts$new(mod, Contr)
y <- contrastX$get_linfct(avg = FALSE)
stopifnot(all(ref_lfc == y))
t <- contrastX$get_linfct(global = FALSE)

x <- contrastX$get_contrasts()
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`
stopifnot(all(x$p.value < 1 & x$p.value > 0))
stopifnot(all(x$avgAbd > 0))
stopifnot(all(x$FDR > 0 & x$FDR < 1))

x <- contrastX$get_contrast_sides()
xd <- contrastX$column_description()
modelFunction <-
  strategy_lm("abundance  ~ group_")
mod <- build_model(
  istar$data,
  modelFunction,
  subject_id = config$hierarchy_keys_depth()
)
contrastX <- prolfqua::Contrasts$new(mod, Contr)
y <- contrastX$get_linfct(avg = FALSE)
stopifnot(all(ref_lfc == y))
```
