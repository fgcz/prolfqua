# Linear mixed-effects model strategy (R6 class)

Linear mixed-effects model strategy (R6 class)

Linear mixed-effects model strategy (R6 class)

## Details

Encapsulates everything needed to fit per-protein linear mixed-effects
models via [`lmer`](https://rdrr.io/pkg/lmerTest/man/lmer.html) and
extract contrasts.

## See also

Other modelling:
[`AnovaExtractor`](https://wolski.github.io/prolfqua/reference/AnovaExtractor.md),
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
[`ContrastsROPECAFacade`](https://wolski.github.io/prolfqua/reference/ContrastsROPECAFacade.md),
[`ContrastsTable`](https://wolski.github.io/prolfqua/reference/ContrastsTable.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`LR_test()`](https://wolski.github.io/prolfqua/reference/LR_test.md),
[`Model`](https://wolski.github.io/prolfqua/reference/Model.md),
[`ModelFirth`](https://wolski.github.io/prolfqua/reference/ModelFirth.md),
[`ModelLimma`](https://wolski.github.io/prolfqua/reference/ModelLimma.md),
[`StrategyLM`](https://wolski.github.io/prolfqua/reference/StrategyLM.md),
[`StrategyLimma`](https://wolski.github.io/prolfqua/reference/StrategyLimma.md),
[`StrategyLogistf`](https://wolski.github.io/prolfqua/reference/StrategyLogistf.md),
[`StrategyRLM`](https://wolski.github.io/prolfqua/reference/StrategyRLM.md),
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

- `formula`:

  model formula

- `model_name`:

  name of model

- `report_columns`:

  columns to report

- `is_mixed`:

  always TRUE for lmer

- `anova_df`:

  list with anova function and column names

## Methods

### Public methods

- [`StrategyLmer$new()`](#method-StrategyLmer-new)

- [`StrategyLmer$model_fun()`](#method-StrategyLmer-model_fun)

- [`StrategyLmer$isSingular()`](#method-StrategyLmer-isSingular)

- [`StrategyLmer$contrast_fun()`](#method-StrategyLmer-contrast_fun)

- [`StrategyLmer$df_residual()`](#method-StrategyLmer-df_residual)

- [`StrategyLmer$sigma()`](#method-StrategyLmer-sigma)

- [`StrategyLmer$clone()`](#method-StrategyLmer-clone)

------------------------------------------------------------------------

### Method [`new()`](https://rdrr.io/r/methods/new.html)

Create a new StrategyLmer

#### Usage

    StrategyLmer$new(
      modelstr,
      model_name = "Model",
      report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value",
        "moderated.p.value.adjusted")
    )

#### Arguments

- `modelstr`:

  model formula string

- `model_name`:

  name of model

- `report_columns`:

  columns to report

------------------------------------------------------------------------

### Method `model_fun()`

Fit lmer to one protein's data

#### Usage

    StrategyLmer$model_fun(x, pb, get_formula = FALSE)

#### Arguments

- `x`:

  data.frame for one protein

- `pb`:

  optional progress bar

- `get_formula`:

  if TRUE, return formula instead of fitting

------------------------------------------------------------------------

### Method `isSingular()`

Check if model is singular

#### Usage

    StrategyLmer$isSingular(model)

#### Arguments

- `model`:

  fitted model

------------------------------------------------------------------------

### Method `contrast_fun()`

Compute contrasts from fitted model

#### Usage

    StrategyLmer$contrast_fun(...)

#### Arguments

- `...`:

  passed to
  [`my_contest`](https://wolski.github.io/prolfqua/reference/my_contest.md)

------------------------------------------------------------------------

### Method `df_residual()`

Get residual degrees of freedom

#### Usage

    StrategyLmer$df_residual(model)

#### Arguments

- `model`:

  fitted model

------------------------------------------------------------------------

### Method [`sigma()`](https://rdrr.io/r/stats/sigma.html)

Get residual standard error

#### Usage

    StrategyLmer$sigma(model)

#### Arguments

- `model`:

  fitted model

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    StrategyLmer$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 10, with_missing = FALSE)
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
istar <- prolfqua::LFQData$new(istar$data, istar$config)
istar$data <- istar$data |> dplyr::group_by(protein_Id) |>
  dplyr::mutate(abundanceC = abundance - mean(abundance)) |> dplyr::ungroup()
strat <- StrategyLmer$new("abundanceC ~ group_ + (1|peptide_Id)",
  model_name = "random_example")
mod <- build_model(istar, strat)
#> Warning: There were 4 warnings in `dplyr::mutate()`.
#> The first warning was:
#> ℹ In argument: `linear_model = purrr::map(data, model_strategy$model_fun, pb =
#>   pb)`.
#> ℹ In group 2: `protein_Id = "7cbcrd~5725"`.
#> Caused by warning in `value[[3L]]()`:
#> ! WARN :Error: grouping factors must have > 1 sampled level
#> ℹ Run `dplyr::last_dplyr_warnings()` to see the 3 remaining warnings.
#> Joining with `by = join_by(protein_Id)`
sum(mod$modelDF$exists_lmer)
#> [1] 6
```
