# Compute contrasts with group mean imputation (DEPRECATED)

Compute contrasts with group mean imputation (DEPRECATED)

Compute contrasts with group mean imputation (DEPRECATED)

## Value

An R6 class generator.

## Details

\`r lifecycle::badge("deprecated")\` (placeholder — prolfqua does not
yet depend on the lifecycle package; consider this a soft deprecation
notice. See *Deprecated* below.)

If there are no observations in one of the groups for some of the
proteins, the group mean cannot be estimated. Therefore, assuming that
the observation is missing because the protein abundance is below the
detection limit, we substitute the unobserved group with the median of
protein abundances observed only in one sample of the group. The
variance of a protein is estimated using the pooled variance of all
observations of all groups.

## Deprecated

\`ContrastsMissing\` is a pre-fitting group-mean substitution: it does
not fit a per-protein model, just stamps the group-mean delta with a
pooled-variance t-test. It is superseded by
[`build_model_impute`](https://wolski.github.io/prolfqua/reference/build_model_impute.md)
which refits failed/singular proteins with LOD-imputed values and
borrowed variance per protein. New code should prefer the `lm_impute`
facade
([`ContrastsLMImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMImputeFacade.md))
or any of its `limma_impute` / `limma_voom_impute` cousins. Constructing
`ContrastsMissing` emits a `.Deprecated` warning at `initialize` time.

Still reachable via `model = "lm_missing"` for users who want to
explicitly run the group-mean fallback; the construction emits a
`.Deprecated` warning each time. Will be removed once a merge-style
replacement (e.g. an `lm_impute_missing` facade composing `build_model`
with `build_model_impute`) lands.

## See also

[`build_model_impute`](https://wolski.github.io/prolfqua/reference/build_model_impute.md),
[`ContrastsLMImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMImputeFacade.md)

Other modelling:
[`AnovaExtractor`](https://wolski.github.io/prolfqua/reference/AnovaExtractor.md),
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsDEqMSVoomFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSVoomFacade.md),
[`ContrastsFirth`](https://wolski.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsFirthFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsLMImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMImputeFacade.md),
[`ContrastsLMMissingFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMMissingFacade.md),
[`ContrastsLimma`](https://wolski.github.io/prolfqua/reference/ContrastsLimma.md),
[`ContrastsLimmaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaFacade.md),
[`ContrastsLimmaImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaImputeFacade.md),
[`ContrastsLimmaVoomFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaVoomFacade.md),
[`ContrastsLimmaVoomImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaVoomImputeFacade.md),
[`ContrastsLimpaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimpaFacade.md),
[`ContrastsLmerFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLmerFacade.md),
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
[`StrategyLimpa`](https://wolski.github.io/prolfqua/reference/StrategyLimpa.md),
[`StrategyLmer`](https://wolski.github.io/prolfqua/reference/StrategyLmer.md),
[`StrategyLogistf`](https://wolski.github.io/prolfqua/reference/StrategyLogistf.md),
[`StrategyRLM`](https://wolski.github.io/prolfqua/reference/StrategyRLM.md),
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
[`new_lm_imputed()`](https://wolski.github.io/prolfqua/reference/new_lm_imputed.md),
[`pivot_model_contrasts_to_wide()`](https://wolski.github.io/prolfqua/reference/pivot_model_contrasts_to_wide.md),
[`plot_lmer_peptide_predictions()`](https://wolski.github.io/prolfqua/reference/plot_lmer_peptide_predictions.md),
[`register_facade()`](https://wolski.github.io/prolfqua/reference/register_facade.md),
[`sim_build_models_lm()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lm.md),
[`sim_build_models_lmer()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lmer.md),
[`sim_build_models_logistf()`](https://wolski.github.io/prolfqua/reference/sim_build_models_logistf.md),
[`sim_make_model_lm()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lm.md),
[`sim_make_model_lmer()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lmer.md),
[`strategy_limma()`](https://wolski.github.io/prolfqua/reference/strategy_limma.md),
[`strategy_limpa()`](https://wolski.github.io/prolfqua/reference/strategy_limpa.md),
[`strategy_logistf()`](https://wolski.github.io/prolfqua/reference/strategy.md),
[`summary_ROPECA_median_p.scaled()`](https://wolski.github.io/prolfqua/reference/summary_ROPECA_median_p.scaled.md),
[`unregister_facade()`](https://wolski.github.io/prolfqua/reference/unregister_facade.md)

## Super class

[`prolfqua::ContrastsInterface`](https://wolski.github.io/prolfqua/reference/ContrastsInterface.md)
-\> `ContrastsMissing`

## Public fields

- `subject_id`:

  subject_id e.g. protein_ID column

- `contrasts`:

  array with contrasts (see example)

- `model_name`:

  model name

- `contrast_result`:

  data frame with results of contrast computation

- `lfqdata`:

  data frame

- `confint`:

  confidence interval

- `p.adjust`:

  function to adjust p-values

- `global`:

  Take global or local values for imputation

- `present`:

  default 1, presence in interaction to infer limit of detection.

- `minsd`:

  default 1, if standard deviation can not be estimated, what is the
  prior minimum sd, default = 1s

## Methods

### Public methods

- [`ContrastsMissing$new()`](#method-ContrastsMissing-new)

- [`ContrastsMissing$get_contrast_sides()`](#method-ContrastsMissing-get_contrast_sides)

- [`ContrastsMissing$get_contrasts()`](#method-ContrastsMissing-get_contrasts)

- [`ContrastsMissing$get_Plotter()`](#method-ContrastsMissing-get_Plotter)

- [`ContrastsMissing$to_wide()`](#method-ContrastsMissing-to_wide)

- [`ContrastsMissing$clone()`](#method-ContrastsMissing-clone)

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

initialize

#### Usage

    ContrastsMissing$new(
      lfqdata,
      contrasts,
      confint = 0.95,
      p.adjust = prolfqua::adjust_p_values,
      model_name = "groupAverage"
    )

#### Arguments

- `lfqdata`:

  LFQData

- `contrasts`:

  array of contrasts (see example)

- `confint`:

  confidence interval

- `p.adjust`:

  method for p-value adjustment - default Benjamini Hochberg

- `model_name`:

  default "groupAverage"

------------------------------------------------------------------------

### Method `get_contrast_sides()`

get contrasts sides

#### Usage

    ContrastsMissing$get_contrast_sides()

------------------------------------------------------------------------

### Method `get_contrasts()`

table with results of contrast computation

#### Usage

    ContrastsMissing$get_contrasts(all = FALSE)

#### Arguments

- `all`:

  FALSE, do not show all columns (default)

------------------------------------------------------------------------

### Method `get_Plotter()`

get ContrastsPlotter

#### Usage

    ContrastsMissing$get_Plotter()

#### Returns

Contrast_Plotter

------------------------------------------------------------------------

### Method `to_wide()`

convert contrast results to wide format

#### Usage

    ContrastsMissing$to_wide(columns = c("p.value", "FDR", "statistic"))

#### Arguments

- `columns`:

  value column default p.value

#### Returns

data.frame

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsMissing$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
Nprot <- 120
istar <- prolfqua::sim_lfq_data_protein_config(Nprot = Nprot,weight_missing = .4)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
istar$data$abundance |> is.na() |> sum()
#> [1] 283
protIntensity <- istar$data
config <- istar$config


lProt <- LFQData$new(protIntensity, config)
lProt$rename_response("transformedIntensity")

Contr <- c("dil.b_vs_a" = "group_A - group_Ctrl")
csi <- ContrastsMissing$new(lProt, contrasts = Contr)
#> Warning: ContrastsMissing is deprecated: it substitutes group means rather than fitting a model. Prefer build_model_impute (LOD-imputed per-protein refit with borrowed variance) via the lm_impute / limma_impute facades. See ?ContrastsMissing for details.
csi$get_contrast_sides()
#> # A tibble: 1 × 3
#>   contrast   group_1 group_2   
#>   <chr>      <chr>   <chr>     
#> 1 dil.b_vs_a group_A group_Ctrl

res <- csi$get_contrasts()
#> dil.b_vs_a=group_A - group_Ctrl
#> dil.b_vs_a=group_A - group_Ctrl
#> dil.b_vs_a=group_A - group_Ctrl

stopifnot(nrow(res) ==  (protIntensity$protein_Id |> unique() |> length()))
res$contrast |> table()
#> 
#> dil.b_vs_a 
#>        120 
stopifnot((res$p.value |> is.na() |> sum()) == 0)
plot(res$diff, -log10(res$p.value), pch = ".")

csi$column_description()
#>           column_name
#> modelName   modelName
#> contrast     contrast
#> avgAbd         avgAbd
#> diff             diff
#> FDR               FDR
#> statistic   statistic
#> std.error   std.error
#> df                 df
#> p.value       p.value
#> conf.low     conf.low
#> conf.high   conf.high
#> sigma           sigma
#>                                                                                            description
#> modelName                                                                                type of model
#> contrast                                                      name of difference e.g. group1_vs_group2
#> avgAbd                                                  mean abundance value of protein in all samples
#> diff                                                                       difference among conditions
#> FDR                                                                               false discovery rate
#> statistic                                                                                 t-statistics
#> std.error                                                                               standard error
#> df                                                                                  degrees of freedom
#> p.value                                                                                        p-value
#> conf.low                                                         lower value of 95 confidence interval
#> conf.high                                                         high value of 95 confidence interval
#> sigma     residual standard deviation of linear model (needed for empirical Bayes variance shrinkage).
x<- csi$get_Plotter()
p <- x$volcano()
pdf(file = NULL)
print(p)
#> $p.value
#> 
#> $FDR
#> 
dev.off()
#> agg_record_2c901db433a8 
#>                       2 

dd <- prolfqua::sim_lfq_data_2factor_config(Nprot = 100,weight_missing = 0.1)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done

Contrasts <- c("c1" = "TreatmentA - TreatmentB",
               "C2" = "BackgroundX- BackgroundZ",
               "c3" = "`TreatmentA:BackgroundX` - `TreatmentA:BackgroundZ`",
               "c4" = "`TreatmentB:BackgroundX` - `TreatmentB:BackgroundZ`"
               )
lProt <- LFQData$new(dd$data, dd$config)
lProt$rename_response("transformedIntensity")

csi <- ContrastsMissing$new(lProt, contrasts = Contrasts)
#> Warning: ContrastsMissing is deprecated: it substitutes group means rather than fitting a model. Prefer build_model_impute (LOD-imputed per-protein refit with borrowed variance) via the lm_impute / limma_impute facades. See ?ContrastsMissing for details.
res <- csi$get_contrasts()
#> c1=TreatmentA - TreatmentB
#> C2=BackgroundX- BackgroundZ
#> c3=`TreatmentA:BackgroundX` - `TreatmentA:BackgroundZ`
#> c4=`TreatmentB:BackgroundX` - `TreatmentB:BackgroundZ`
#> c1=TreatmentA - TreatmentB
#> C2=BackgroundX- BackgroundZ
#> c3=`TreatmentA:BackgroundX` - `TreatmentA:BackgroundZ`
#> c4=`TreatmentB:BackgroundX` - `TreatmentB:BackgroundZ`
#> c1=TreatmentA - TreatmentB
#> C2=BackgroundX- BackgroundZ
#> c3=`TreatmentA:BackgroundX` - `TreatmentA:BackgroundZ`
#> c4=`TreatmentB:BackgroundX` - `TreatmentB:BackgroundZ`
pl <- csi$get_Plotter()
pdf(file = NULL)
pl$volcano()
#> $p.value
#> 
#> $FDR
#> 
dev.off()
#> agg_record_2c901db433a8 
#>                       2 
```
