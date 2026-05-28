# Build a contrast analysis using one of several statistical methods

A builder function that dispatches to the appropriate facade class based
on the chosen method. Each facade encapsulates the full pipeline from
strategy construction through modelling to contrast computation.

## Usage

``` r
build_contrast_analysis(
  lfqdata,
  modelstr,
  contrasts,
  method = c("lm", "lm_impute", "lm_missing", "limma", "limma_impute", "limma_voom",
    "limma_voom_impute", "limpa", "limpa_nested", "rlm", "rfit", "deqms", "deqms_voom",
    "firth", "firth_nested", "lmer_nested", "ropeca_nested"),
  ...
)
```

## Arguments

- lfqdata:

  an [`LFQData`](https://wolski.github.io/prolfqua/reference/LFQData.md)
  object

- modelstr:

  model formula string without the response variable (e.g.
  `"~ group_"`). The response is taken automatically from
  `lfqdata$get_config()$get_response()`.

- contrasts:

  named character vector of contrasts (e.g.
  `c("A_vs_B" = "group_A - group_B")`)

- method:

  one of `"lm"`, `"lm_impute"`, `"lm_missing"`, `"limma"`,
  `"limma_impute"`, `"limma_voom"`, `"limma_voom_impute"`, `"limpa"`,
  `"limpa_nested"`, `"rlm"`, `"deqms"`, `"deqms_voom"`, `"firth"`,
  `"firth_nested"`, `"lmer_nested"`, `"ropeca_nested"`

- ...:

  additional arguments forwarded to the underlying strategy function
  (e.g. `trend`, `robust` for `strategy_limma`)

## Value

one of
[`ContrastsLimmaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimmaFacade.md),
[`ContrastsLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsRLMFacade`](https://wolski.github.io/prolfqua/reference/ContrastsRLMFacade.md),
[`ContrastsLmerNestedFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLmerNestedFacade.md),
[`ContrastsLMMissingFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMMissingFacade.md),
[`ContrastsLMImputeFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLMImputeFacade.md),
[`ContrastsDEqMSFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsROPECANestedFacade`](https://wolski.github.io/prolfqua/reference/ContrastsROPECANestedFacade.md),
[`ContrastsFirthFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsFirthNestedFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthNestedFacade.md),
[`ContrastsLimpaFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimpaFacade.md),
or
[`ContrastsLimpaNestedFacade`](https://wolski.github.io/prolfqua/reference/ContrastsLimpaNestedFacade.md)

## Vectorized mode

Set `options(prolfqua.vectorize = TRUE)` before calling this function to
activate vectorized implementations of
[`compute_contrast`](https://wolski.github.io/prolfqua/reference/compute_contrast.md)
and
[`linfct_matrix_contrasts`](https://wolski.github.io/prolfqua/reference/linfct_matrix_contrasts.md).
This affects all methods that use the Wald test path (lm, rlm, firth,
lmer) and can give a significant speed-up for large datasets. Results
are numerically identical. Example:

    options(prolfqua.vectorize = TRUE)
    fa <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "lm")
    options(prolfqua.vectorize = FALSE)  # restore default

## See also

Other modelling:
[`AnovaExtractor`](https://wolski.github.io/prolfqua/reference/AnovaExtractor.md),
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsDEqMSVoomFacade`](https://wolski.github.io/prolfqua/reference/ContrastsDEqMSVoomFacade.md),
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
[`new_lm_imputed()`](https://wolski.github.io/prolfqua/reference/new_lm_imputed.md),
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

## Examples

``` r
istar <- sim_lfq_data_protein_config(Nprot = 20)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfqdata <- LFQData$new(istar$data, istar$config)
lfqdata$rename_response("transformedIntensity")
contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

fa_lm    <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "lm")
head(fa_lm$get_contrasts())
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`
#> # A tibble: 6 × 14
#>   facade modelName  protein_Id contrast    diff std.error avgAbd statistic    df
#>   <chr>  <chr>      <chr>      <chr>      <dbl>     <dbl>  <dbl>     <dbl> <dbl>
#> 1 lm     WaldTest_… 0EfVhX~59… A_vs_Ct…  2.72       1.14    23.2     2.47  12.1 
#> 2 lm     WaldTest_… 0m5WN4~14… A_vs_Ct…  0.600      0.734   17.4     0.765  8.08
#> 3 lm     WaldTest_… 7cbcrd~83… A_vs_Ct…  2.59       0.572   27.0     3.68  12.1 
#> 4 lm     WaldTest_… 9VUkAq~45… A_vs_Ct…  0.0679     0.760   19.4     0.104 10.1 
#> 5 lm     WaldTest_… At886V~32… A_vs_Ct… -1.01       0.969   19.1    -1.20   9.08
#> 6 lm     WaldTest_… BEJI92~91… A_vs_Ct… -0.873      0.659   20.9    -1.39  11.1 
#> # ℹ 5 more variables: p.value <dbl>, conf.low <dbl>, conf.high <dbl>,
#> #   sigma <dbl>, FDR <dbl>

fa_limma <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "limma")
head(fa_limma$get_contrasts())
#> # A tibble: 6 × 14
#>   facade modelName protein_Id  contrast     diff    FDR std.error statistic
#>   <chr>  <chr>     <chr>       <chr>       <dbl>  <dbl>     <dbl>     <dbl>
#> 1 limma  limma     0EfVhX~5954 A_vs_Ctrl  2.72   0.188      1.09      2.49 
#> 2 limma  limma     0m5WN4~1448 A_vs_Ctrl  0.600  0.623      0.770     0.779
#> 3 limma  limma     7cbcrd~8305 A_vs_Ctrl  2.59   0.0271     0.691     3.75 
#> 4 limma  limma     9VUkAq~4562 A_vs_Ctrl  0.0679 0.967      0.647     0.105
#> 5 limma  limma     At886V~3296 A_vs_Ctrl -1.01   0.623      0.836    -1.21 
#> 6 limma  limma     BEJI92~9143 A_vs_Ctrl -0.873  0.623      0.621    -1.41 
#> # ℹ 6 more variables: p.value <dbl>, sigma <dbl>, df <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, avgAbd <dbl>

fa_miss <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "lm_missing")
#> Warning: ContrastsLMMissingFacade (method = 'lm_missing') is deprecated: its second leg uses ContrastsMissing (group-mean substitution, no model fit). Prefer 'lm_impute' which refits failed/singular proteins with LOD imputation and borrowed variance, tagging rescued rows as 'WaldTest_moderated_imputed'. See ?ContrastsLMMissingFacade for migration.
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`
#> A_vs_Ctrl=group_A - group_Ctrl
#> A_vs_Ctrl=group_A - group_Ctrl
#> A_vs_Ctrl=group_A - group_Ctrl
#> Joining with `by = join_by(protein_Id, contrast)`
#> Joining with `by = join_by(protein_Id, contrast)`
head(fa_miss$get_contrasts())
#> # A tibble: 6 × 14
#>   facade  modelName protein_Id contrast    diff std.error avgAbd statistic    df
#>   <chr>   <fct>     <chr>      <chr>      <dbl>     <dbl>  <dbl>     <dbl> <dbl>
#> 1 lm_mis… WaldTest… 0EfVhX~59… A_vs_Ct…  2.72       1.14    23.2     2.47  12.1 
#> 2 lm_mis… WaldTest… 0m5WN4~14… A_vs_Ct…  0.600      0.734   17.4     0.765  8.08
#> 3 lm_mis… WaldTest… 7cbcrd~83… A_vs_Ct…  2.59       0.572   27.0     3.68  12.1 
#> 4 lm_mis… WaldTest… 9VUkAq~45… A_vs_Ct…  0.0679     0.760   19.4     0.104 10.1 
#> 5 lm_mis… WaldTest… At886V~32… A_vs_Ct… -1.01       0.969   19.1    -1.20   9.08
#> 6 lm_mis… WaldTest… BEJI92~91… A_vs_Ct… -0.873      0.659   20.9    -1.39  11.1 
#> # ℹ 5 more variables: p.value <dbl>, conf.low <dbl>, conf.high <dbl>,
#> #   sigma <dbl>, FDR <dbl>

fa_deqms <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "deqms")
head(fa_deqms$get_contrasts())
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`
#> # A tibble: 6 × 14
#>   facade contrast  modelName protein_Id    diff std.error avgAbd statistic    df
#>   <chr>  <chr>     <chr>     <chr>        <dbl>     <dbl>  <dbl>     <dbl> <int>
#> 1 deqms  A_vs_Ctrl WaldTest… 0EfVhX~59…  2.72       1.14    23.2    4.23       9
#> 2 deqms  A_vs_Ctrl WaldTest… 0m5WN4~14…  0.600      0.734   17.4    0.619      5
#> 3 deqms  A_vs_Ctrl WaldTest… 7cbcrd~83…  2.59       0.572   27.0    4.02       9
#> 4 deqms  A_vs_Ctrl WaldTest… 9VUkAq~45…  0.0679     0.760   19.4    0.0837     7
#> 5 deqms  A_vs_Ctrl WaldTest… At886V~32… -1.01       0.969   19.1   -1.37       6
#> 6 deqms  A_vs_Ctrl WaldTest… BEJI92~91… -0.873      0.659   20.9   -1.31       8
#> # ℹ 5 more variables: p.value <dbl>, conf.low <dbl>, conf.high <dbl>,
#> #   sigma <dbl>, FDR <dbl>

istar_pep <- sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfqdata_pep <- LFQData$new(istar_pep$data, istar_pep$config)
lfqdata_pep <- lfqdata_pep$get_Transformer()$log2()$lfq
#> Column added : log2_abundance

fa_lmer <- build_contrast_analysis(
  lfqdata_pep,
  "~ group_ + (1 | peptide_Id) + (1 | sampleName)",
  contrasts,
  method = "lmer_nested"
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
head(fa_lmer$get_contrasts())
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`
#> # A tibble: 6 × 14
#>   facade modelName protein_Id contrast     diff std.error avgAbd statistic    df
#>   <chr>  <chr>     <chr>      <chr>       <dbl>     <dbl>  <dbl>     <dbl> <dbl>
#> 1 lmer_… WaldTest… 0EfVhX~00… A_vs_Ct… -8.32e-4    0.0730   4.34   -0.0115  28.9
#> 2 lmer_… WaldTest… BEJI92~52… A_vs_Ct…  3.22e-1    0.0832   4.22    2.81    11.6
#> 3 lmer_… WaldTest… Fl4JiV~86… A_vs_Ct… -4.13e-2    0.0850   4.38   -0.503   39.5
#> 4 lmer_… WaldTest… HvIpHG~90… A_vs_Ct… -3.72e-1    0.0616   4.40   -5.65    21.8
#> 5 lmer_… WaldTest… JcKVfU~96… A_vs_Ct… -1.07e-1    0.0577   5.05   -1.88    79.8
#> 6 lmer_… WaldTest… SGIVBl~57… A_vs_Ct…  3.07e-2    0.0695   4.68    0.452   61.0
#> # ℹ 5 more variables: p.value <dbl>, conf.low <dbl>, conf.high <dbl>,
#> #   sigma <dbl>, FDR <dbl>

fa_ropeca <- build_contrast_analysis(lfqdata_pep, "~ group_", contrasts, method = "ropeca_nested")
head(fa_ropeca$get_contrasts())
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, peptide_Id, contrast)`
#> # A tibble: 6 × 14
#> # Groups:   contrast [1]
#>   facade        protein_Id  modelName contrast  avgAbd    diff     FDR statistic
#>   <chr>         <chr>       <chr>     <chr>      <dbl>   <dbl>   <dbl>     <dbl>
#> 1 ropeca_nested 0EfVhX~0087 ROPECA    A_vs_Ctrl   4.27 -0.0742 5.28e-2     -1.75
#> 2 ropeca_nested 7cbcrd~5725 ROPECA    A_vs_Ctrl   4.51  0.741  9.91e-5      8.79
#> 3 ropeca_nested 9VUkAq~4703 ROPECA    A_vs_Ctrl   4.47 -0.598  6.91e-6    -12.7 
#> 4 ropeca_nested BEJI92~5282 ROPECA    A_vs_Ctrl   4.23  0.277  1.87e-3      3.94
#> 5 ropeca_nested CGzoYe~2147 ROPECA    A_vs_Ctrl   4.76 -0.310  3.74e-5     -9.26
#> 6 ropeca_nested DoWup2~5896 ROPECA    A_vs_Ctrl   4.43  0.295  1.38e-6     14.7 
#> # ℹ 6 more variables: std.error <dbl>, df <int>, p.value <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, sigma <dbl>

fa_firth <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "firth")
#> completing cases
#> Joining with `by = join_by(protein_Id)`
head(fa_firth$get_contrasts())
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct_firth
#> Joining with `by = join_by(protein_Id, contrast)`
#> # A tibble: 6 × 14
#> # Groups:   contrast [1]
#>   facade modelName     protein_Id contrast sigma    df      diff   FDR std.error
#>   <chr>  <chr>         <chr>      <chr>    <dbl> <int>     <dbl> <dbl>     <dbl>
#> 1 firth  WaldTestFirth 0EfVhX~59… A_vs_Ct…     1     9  1.07e-15     1      2.11
#> 2 firth  WaldTestFirth 0m5WN4~14… A_vs_Ct…     1     9 -2.20e+ 0     1      1.74
#> 3 firth  WaldTestFirth 7cbcrd~83… A_vs_Ct…     1     9  1.07e-15     1      2.11
#> 4 firth  WaldTestFirth 9VUkAq~45… A_vs_Ct…     1     9 -1.35e+ 0     1      1.78
#> 5 firth  WaldTestFirth At886V~32… A_vs_Ct…     1     9  5.58e-16     1      1.38
#> 6 firth  WaldTestFirth BEJI92~91… A_vs_Ct…     1     9 -1.35e+ 0     1      1.78
#> # ℹ 5 more variables: statistic <dbl>, p.value <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, avgAbd <dbl>
```
