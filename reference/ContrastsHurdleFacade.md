# Hurdle contrasts: intensity and detection combined with Fisher's method

Hurdle contrasts: intensity and detection combined with Fisher's method

Hurdle contrasts: intensity and detection combined with Fisher's method

## Value

An R6 class generator.

## Details

Joins the contrasts of a fitted intensity facade (e.g. `lm_impute`) and
a fitted detection facade (`firth_nested`, `binomial_nested` or
`binomial_firth_nested`) by subject and contrast, and replaces `p.value`
and `FDR` with Fisher's combined test, as in msqrob2's
`hypothesisTestHurdle`. No model is fitted.

For each subject and contrast, \\X^2 = -2 \sum \log p\\ over the
available component p-values, with two degrees of freedom per component,
so a row with a single component keeps that component's p-value. The
detection p-value is left out for subjects whose child features are
detected in every sample, because detection carries no information
there. `FDR` is the Benjamini-Hochberg adjustment within each contrast.

`diff`, `avgAbd`, `std.error`, `statistic`, `df`, the confidence
interval, `sigma` and `estimate_type` come from the intensity component.
They are `NA` for subjects only the detection component estimates; with
`lm` these are the proteins detected in one group only, so
`filter_significant()` and `get_ora()` drop them. Pair with `lm_impute`,
which estimates their `diff`.

Added columns: `p.value.intensity`, `FDR.intensity`, `diff.detection`
(log odds ratio), `p.value.detection`, `FDR.detection` (all as reported
by the components), `chisq.hurdle`, `df.hurdle`, and `hurdle_components`
(`both`, `intensity` or `detection`: the p-values combined).

## See also

Other modelling:
[`AnovaExtractor`](https://wolski.github.io/prolfqua/reference/AnovaExtractor.md),
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsDEqMSVoomFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSVoomFacade.md),
[`ContrastsFacadeBase`](https://wolski.github.io/prolfqua/reference/ContrastsFacadeBase.md),
[`ContrastsFirth`](https://wolski.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsFirthFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsFirthNestedFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthNestedFacade.md),
[`ContrastsHurdleNestedFacade`](https://wolski.github.io/prolfqua/reference/ContrastsHurdleNestedFacade.md),
[`ContrastsLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsLMImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMImputeFacade.md),
[`ContrastsLimma`](https://wolski.github.io/prolfqua/reference/ContrastsLimma.md),
[`ContrastsLimmaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaFacade.md),
[`ContrastsLimmaImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaImputeFacade.md),
[`ContrastsLimmaVoomFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaVoomFacade.md),
[`ContrastsLimmaVoomImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaVoomImputeFacade.md),
[`ContrastsLimpaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimpaFacade.md),
[`ContrastsLimpaNestedFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimpaNestedFacade.md),
[`ContrastsLmerNestedFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLmerNestedFacade.md),
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
[`df.residual.rfit_prolfqua()`](https://wolski.github.io/prolfqua/reference/df.residual.rfit_prolfqua.md),
[`get_anova_df()`](https://wolski.github.io/prolfqua/reference/get_anova_df.md),
[`get_complete_model_fit()`](https://wolski.github.io/prolfqua/reference/get_complete_model_fit.md),
[`get_p_values_pbeta()`](https://wolski.github.io/prolfqua/reference/get_p_values_pbeta.md),
[`group_label()`](https://wolski.github.io/prolfqua/reference/group_label.md),
[`impute_from_model()`](https://wolski.github.io/prolfqua/reference/impute_from_model.md),
[`impute_refit_singular()`](https://wolski.github.io/prolfqua/reference/impute_refit_singular.md),
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
[`new_imputed_model()`](https://wolski.github.io/prolfqua/reference/new_imputed_model.md),
[`pivot_model_contrasts_to_wide()`](https://wolski.github.io/prolfqua/reference/pivot_model_contrasts_to_wide.md),
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
[`summary_ROPECA_median_p_scaled()`](https://wolski.github.io/prolfqua/reference/summary_ROPECA_median_p_scaled.md),
[`vcov.rfit_prolfqua()`](https://wolski.github.io/prolfqua/reference/vcov.rfit_prolfqua.md)

## Super class

[`prolfqua::ContrastsInterface`](https://wolski.github.io/prolfqua/reference/ContrastsInterface.md)
-\> `ContrastsHurdleFacade`

## Public fields

- `intensity`:

  fitted intensity facade

- `detection`:

  fitted detection facade (`firth_nested`, `binomial_nested` or
  `binomial_firth_nested`)

- `model_name`:

  name stamped into the `modelName` column

- `subject_id`:

  columns identifying the subject (e.g. protein_Id)

## Methods

### Public methods

- [`ContrastsHurdleFacade$new()`](#method-ContrastsHurdleFacade-new)

- [`ContrastsHurdleFacade$get_contrasts()`](#method-ContrastsHurdleFacade-get_contrasts)

- [`ContrastsHurdleFacade$clone()`](#method-ContrastsHurdleFacade-clone)

Inherited methods

- [`prolfqua::ContrastsInterface$column_description()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-column_description)
- [`prolfqua::ContrastsInterface$contrast_summary_table()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-contrast_summary_table)
- [`prolfqua::ContrastsInterface$extra_artifacts()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-extra_artifacts)
- [`prolfqua::ContrastsInterface$filter_significant()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-filter_significant)
- [`prolfqua::ContrastsInterface$get_Plotter()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_Plotter)
- [`prolfqua::ContrastsInterface$get_config()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_config)
- [`prolfqua::ContrastsInterface$get_contrast_sides()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_contrast_sides)
- [`prolfqua::ContrastsInterface$get_missing()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_missing)
- [`prolfqua::ContrastsInterface$get_ora()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_ora)
- [`prolfqua::ContrastsInterface$get_rank()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_rank)
- [`prolfqua::ContrastsInterface$to_wide()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-to_wide)

------------------------------------------------------------------------

### Method `new()`

Combine a fitted intensity facade and a fitted detection facade.

#### Usage

    ContrastsHurdleFacade$new(intensity, detection, model_name = "hurdle")

#### Arguments

- `intensity`:

  fitted intensity facade, e.g. from
  `build_contrast_analysis(method = "lm_impute")`

- `detection`:

  fitted `firth_nested`, `binomial_nested` or `binomial_firth_nested`
  facade with the same subject id and contrasts

- `model_name`:

  name stamped into the `modelName` column

------------------------------------------------------------------------

### Method `get_contrasts()`

get contrast results with Fisher-combined `p.value` and `FDR`

#### Usage

    ContrastsHurdleFacade$get_contrasts()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsHurdleFacade$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- sim_lfq_data_peptide_config(Nprot = 20, weight_missing = 0.5, seed = 3)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq_peptide <- LFQData$new(istar$data, istar$config)
lfq_peptide <- lfq_peptide$get_Transformer()$log2()$lfq
#> Column added : log2_abundance
lfq_protein <- lfq_peptide$get_Aggregator()$aggregate()
#> starting aggregation
#> Warning: medpolish() did not converge in 10 iterations
#> completing cases
contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
intensity <- build_contrast_analysis(lfq_protein, "~ group_", contrasts, method = "lm_impute")
detection <- build_contrast_analysis(lfq_peptide, "~ group_", contrasts, method = "binomial_nested")
#> completing cases
hurdle <- ContrastsHurdleFacade$new(intensity, detection)
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`
res <- hurdle$get_contrasts()
#> completing cases
table(res$hurdle_components)
#> 
#>      both intensity 
#>        17         3 
head(res)
#> # A tibble: 6 × 24
#>   modelName estimate_type protein_Id  contrast     diff avgAbd
#>   <chr>     <chr>         <chr>       <chr>       <dbl>  <dbl>
#> 1 hurdle    observed      0GRprF~7339 A_vs_Ctrl -0.231    4.03
#> 2 hurdle    observed      4JK499~3111 A_vs_Ctrl  0.466    3.90
#> 3 hurdle    observed      7IZdVV~6818 A_vs_Ctrl  0.122    4.62
#> 4 hurdle    observed      AZPG26~9461 A_vs_Ctrl -0.0544   4.67
#> 5 hurdle    observed      AoNKbb~3497 A_vs_Ctrl  0.0672   4.52
#> 6 hurdle    observed      CibL2O~2149 A_vs_Ctrl -0.328    4.04
#> # ℹ 18 more variables: std.error.unmoderated <dbl>, df.unmoderated <int>,
#> #   std.error <dbl>, statistic <dbl>, df <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, sigma <dbl>, p.value <dbl>, FDR <dbl>,
#> #   p.value.intensity <dbl>, FDR.intensity <dbl>, diff.detection <dbl>,
#> #   p.value.detection <dbl>, FDR.detection <dbl>, chisq.hurdle <dbl>,
#> #   df.hurdle <dbl>, hurdle_components <chr>
hurdle$to_wide()
#> completing cases
#> # A tibble: 20 × 5
#>    protein_Id diff.A_vs_Ctrl p.value.A_vs_Ctrl FDR.A_vs_Ctrl statistic.A_vs_Ctrl
#>    <chr>               <dbl>             <dbl>         <dbl>               <dbl>
#>  1 0GRprF~73…       -0.231       0.0569         0.142                     -2.74 
#>  2 4JK499~31…        0.466       0.00543        0.0362                     4.45 
#>  3 7IZdVV~68…        0.122       0.0231         0.0769                     3.09 
#>  4 AZPG26~94…       -0.0544      0.263          0.405                     -1.16 
#>  5 AoNKbb~34…        0.0672      0.190          0.316                      1.76 
#>  6 CibL2O~21…       -0.328       0.00959        0.0480                    -2.94 
#>  7 EnbuYk~78…       -0.0358      0.922          0.922                     -0.490
#>  8 JnEDsz~42…        0.0378      0.101          0.202                      1.94 
#>  9 KIBglg~82…       -0.0206      0.657          0.773                     -0.485
#> 10 NNASGK~62…       -0.136       0.0155         0.0620                    -2.74 
#> 11 XYrp6h~48…       -0.0548      0.592          0.740                     -0.927
#> 12 lXM5ll~28…        0.144       0.0390         0.111                      2.76 
#> 13 mkmeN0~24…        0.114       0.135          0.246                      2.40 
#> 14 o3n6lI~39…        0.0536      0.725          0.806                      0.705
#> 15 quTD7H~92…        0.140       0.403          0.575                      1.68 
#> 16 tCZCHm~57…        0.609       0.00000000134  0.0000000268              13.3  
#> 17 tHE075~60…        0.0524      0.0659         0.146                      2.65 
#> 18 uhofDS~09…        0.142       0.000491       0.00491                    5.45 
#> 19 uxVYAA~90…        0.00793     0.853          0.898                      0.175
#> 20 xpDUdv~61…        0.0521      0.498          0.664                      0.824
```
