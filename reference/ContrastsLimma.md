# Limma-based contrasts (direct limma pipeline)

Limma-based contrasts (direct limma pipeline)

Limma-based contrasts (direct limma pipeline)

## Details

Uses limma's `contrasts.fit` + `eBayes` pipeline directly, rather than
fitting per-protein lm models and then moderating. Inherits from
[`ContrastsInterface`](https://wolski.github.io/prolfqua/reference/ContrastsInterface.md)
with the same API as
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md).

## See also

Other modelling:
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsFirth`](https://wolski.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsFirthFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsLMImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMImputeFacade.md),
[`ContrastsLMMissingFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMMissingFacade.md),
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

## Super class

[`prolfqua::ContrastsInterface`](https://wolski.github.io/prolfqua/reference/ContrastsInterface.md)
-\> `ContrastsLimma`

## Public fields

- `model`:

  ModelLimma object

- `contrasts`:

  named character vector of contrasts

- `modelName`:

  model name

- `subject_Id`:

  columns with subject_Id (proteinID)

- `p.adjust`:

  function to adjust p-values

- `contrast_result`:

  cached contrast results

## Methods

### Public methods

- [`ContrastsLimma$new()`](#method-ContrastsLimma-new)

- [`ContrastsLimma$get_contrast_sides()`](#method-ContrastsLimma-get_contrast_sides)

- [`ContrastsLimma$get_linfct()`](#method-ContrastsLimma-get_linfct)

- [`ContrastsLimma$get_contrasts()`](#method-ContrastsLimma-get_contrasts)

- [`ContrastsLimma$get_Plotter()`](#method-ContrastsLimma-get_Plotter)

- [`ContrastsLimma$to_wide()`](#method-ContrastsLimma-to_wide)

- [`ContrastsLimma$clone()`](#method-ContrastsLimma-clone)

Inherited methods

- [`prolfqua::ContrastsInterface$column_description()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-column_description)
- [`prolfqua::ContrastsInterface$get_missing()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_missing)

------------------------------------------------------------------------

### Method `new()`

initialize ContrastsLimma

#### Usage

    ContrastsLimma$new(
      model,
      contrasts,
      p.adjust = prolfqua::adjust_p_values,
      modelName = "limma"
    )

#### Arguments

- `model`:

  a
  [`ModelLimma`](https://wolski.github.io/prolfqua/reference/ModelLimma.md)
  object

- `contrasts`:

  named character vector of contrasts

- `p.adjust`:

  function to adjust p-values

- `modelName`:

  name of the contrast method

------------------------------------------------------------------------

### Method `get_contrast_sides()`

get both sides of contrasts

#### Usage

    ContrastsLimma$get_contrast_sides()

------------------------------------------------------------------------

### Method `get_linfct()`

get linear functions from contrasts

#### Usage

    ContrastsLimma$get_linfct(global = TRUE, avg = TRUE)

#### Arguments

- `global`:

  ignored (for API compatibility)

- `avg`:

  logical, also compute avgAbd linfct

------------------------------------------------------------------------

### Method `get_contrasts()`

get table with contrast estimates via limma pipeline

#### Usage

    ContrastsLimma$get_contrasts(all = FALSE)

#### Arguments

- `all`:

  should all columns be returned (default FALSE)

#### Returns

data.frame with contrasts

------------------------------------------------------------------------

### Method `get_Plotter()`

return
[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.md)

#### Usage

    ContrastsLimma$get_Plotter(FCthreshold = 1, FDRthreshold = 0.1)

#### Arguments

- `FCthreshold`:

  fold change threshold to show in plots

- `FDRthreshold`:

  FDR threshold to show in plots

#### Returns

[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.md)

------------------------------------------------------------------------

### Method `to_wide()`

convert to wide format

#### Usage

    ContrastsLimma$to_wide(columns = c("p.value", "FDR", "statistic"))

#### Arguments

- `columns`:

  value column default p.value

#### Returns

data.frame

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsLimma$clone(deep = FALSE)

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
lProt <- LFQData$new(istar$data, istar$config)
lProt$rename_response("transformedIntensity")

strat <- strategy_limma("transformedIntensity ~ group_")
mod_limma <- build_model_limma(lProt, strat)
#> Warning: Partial NA coefficients for 1 probe(s)

Contr <- c("dil.b_vs_a" = "group_A - group_Ctrl")
contr_limma <- ContrastsLimma$new(mod_limma, Contr)
res <- contr_limma$get_contrasts()
head(res)
#> # A tibble: 6 × 13
#>   modelName protein_Id contrast   diff     FDR std.error statistic p.value sigma
#>   <chr>     <chr>      <chr>     <dbl>   <dbl>     <dbl>     <dbl>   <dbl> <dbl>
#> 1 limma     0EfVhX~71… dil.b_v…  3.00  6.99e-4     0.681     4.40  1.43e-5 0.963
#> 2 limma     0m5WN4~35… dil.b_v…  0.222 8.35e-1     0.736     0.301 7.63e-1 0.963
#> 3 limma     76k03k~97… dil.b_v…  0.509 8.35e-1     0.681     0.747 4.55e-1 0.963
#> 4 limma     7QuTub~55… dil.b_v… -1.22  4.39e-1     0.736    -1.66  9.69e-2 0.963
#> 5 limma     7cbcrd~04… dil.b_v…  1.38  5.47e-1     0.963     1.44  1.51e-1 0.963
#> 6 limma     7soopj~34… dil.b_v…  0.822 5.88e-1     0.681     1.21  2.28e-1 0.963
#> # ℹ 4 more variables: df <dbl>, conf.low <dbl>, conf.high <dbl>, avgAbd <dbl>
stopifnot(all(c("diff", "FDR", "p.value", "statistic") %in% colnames(res)))

# Compare with prolfqua's own pipeline
modelFunction <- strategy_lm("transformedIntensity ~ group_")
mod <- build_model(lProt, modelFunction)
#> Joining with `by = join_by(protein_Id)`
contr_prolfqua <- Contrasts$new(mod, Contr)
res_prolfqua <- contr_prolfqua$get_contrasts()
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`

# fold changes should be very similar
merged <- dplyr::inner_join(
  dplyr::select(res, protein_Id, diff_limma = diff),
  dplyr::select(res_prolfqua, protein_Id, diff_prolfqua = diff),
  by = "protein_Id")
#> Adding missing grouping variables: `contrast`
stopifnot(cor(merged$diff_limma, merged$diff_prolfqua, use = "complete.obs") > 0.99)

# Plotter works
pl <- contr_limma$get_Plotter()

# to_wide works
wide <- contr_limma$to_wide()
head(wide)
#> # A tibble: 6 × 5
#>   protein_Id  diff.dil.b_vs_a p.value.dil.b_vs_a FDR.dil.b_vs_a
#>   <chr>                 <dbl>              <dbl>          <dbl>
#> 1 0EfVhX~7161           3.00           0.0000143       0.000699
#> 2 0m5WN4~3543           0.222          0.763           0.835   
#> 3 76k03k~9735           0.509          0.455           0.835   
#> 4 7QuTub~5556          -1.22           0.0969          0.439   
#> 5 7cbcrd~0495           1.38           0.151           0.547   
#> 6 7soopj~3451           0.822          0.228           0.588   
#> # ℹ 1 more variable: statistic.dil.b_vs_a <dbl>

# merge_contrasts_results works
Contr2 <- c("dil.b_vs_a" = "group_A - group_Ctrl")
csi <- ContrastsMissing$new(lProt, contrasts = Contr2)
merged_res <- merge_contrasts_results(contr_limma, csi)
#> completing cases
#> dil.b_vs_a=group_A - group_Ctrl
#> dil.b_vs_a=group_A - group_Ctrl
#> dil.b_vs_a=group_A - group_Ctrl
#> Joining with `by = join_by(protein_Id, contrast)`
#> Joining with `by = join_by(protein_Id, contrast)`
```
