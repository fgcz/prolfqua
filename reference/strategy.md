# Create Firth's logistic regression strategy

Convenience wrapper that creates a
[`StrategyLogistf`](https://fgcz.github.io/prolfqua/reference/StrategyLogistf.md)
object.

Convenience wrapper that creates a
[`StrategyRfit`](https://fgcz.github.io/prolfqua/reference/StrategyRfit.md)
object.

Convenience wrapper that creates a
[`StrategyLmer`](https://fgcz.github.io/prolfqua/reference/StrategyLmer.md)
object.

Convenience wrapper that creates a
[`StrategyLM`](https://fgcz.github.io/prolfqua/reference/StrategyLM.md)
object.

Convenience wrapper that creates a
[`StrategyRLM`](https://fgcz.github.io/prolfqua/reference/StrategyRLM.md)
object.

## Usage

``` r
strategy_logistf(modelstr, model_name = "logistf", control = list())

strategy_rfit(modelstr, model_name = "rfit", control = list())

strategy_lmer(modelstr, model_name = "Model", control = list())

strategy_lm(modelstr, model_name = "Model", weights = NULL)

strategy_rlm(modelstr, model_name = "Model", control = list())
```

## Arguments

- modelstr:

  model formula

- model_name:

  name of model

- control:

  for `strategy_rlm()`: `maxit` (default 20) and `acc` (default 1e-4),
  the iteration limit and the relative change of the IRLS residuals that
  stops [`rlm`](https://rdrr.io/pkg/MASS/man/rlm.html); the criterion is
  relative, so it does not depend on the intensity scale, and 1e-4 is
  far below the precision of log2 intensities.

- weights:

  optional character string naming a column in the data containing
  per-observation weights, passed to
  [`lm`](https://rdrr.io/r/stats/lm.html). Default `NULL` (unweighted).

## Value

a
[`StrategyLogistf`](https://fgcz.github.io/prolfqua/reference/StrategyLogistf.md)
object

a
[`StrategyRfit`](https://fgcz.github.io/prolfqua/reference/StrategyRfit.md)
object

a
[`StrategyLmer`](https://fgcz.github.io/prolfqua/reference/StrategyLmer.md)
object

a
[`StrategyLM`](https://fgcz.github.io/prolfqua/reference/StrategyLM.md)
object

a
[`StrategyRLM`](https://fgcz.github.io/prolfqua/reference/StrategyRLM.md)
object

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
[`vcov.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/vcov.rfit_prolfqua.md)

## Examples

``` r
tmp <- strategy_logistf("bin_resp ~ condition", model_name = "parallel design")
tmp$formula
#> bin_resp ~ condition
#> <environment: 0x615a92fc5928>

istar <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 10, with_missing = TRUE,
  weight_missing = 0.5, seed = 3)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 60 rows with NA response for feature and sample combinations absent from the input.
lfq <- LFQData$new(istar$data, istar$config)
detection <- lfq$data_long()
detection$bin_resp <- as.integer(!is.na(detection[[lfq$response()]]))
config <- lfq$get_config()
config$bin_resp <- "bin_resp"
istar <- LFQData$new(detection, config)
df <- istar$summarize_hierarchy()
df2 <- df[df[[ncol(df)]] > 1, ]
istar2 <- istar$get_subset(df2)
istar2$data_long() |>
  dplyr::group_by(protein_Id) |>
  tidyr::nest() -> nestProtein
modelFunction <- strategy_logistf("bin_resp ~ group_ + peptide_Id",
  model_name = "random_example")
modelFunction$model_fun(nestProtein$data[[1]])
#> logistf::logistf(formula = self$formula, data = DFT, pl = FALSE, 
#>     control = do.call(logistf::logistf.control, self$control), 
#>     weights = Freq)
#> Model fitted by Penalized ML
#> Confidence intervals and p-values by Wald 
#> 
#> Coefficients:
#>        (Intercept)            group_B         group_Ctrl peptide_IdFLq7LKTq 
#>       2.101899e+00       6.773389e-01      -6.663876e-01       9.509484e-16 
#> peptide_IdJYhOpuPH peptide_IdLiw5EMKP peptide_IdVcatZJTa peptide_IdjrLUqOjg 
#>      -1.068335e+00       7.971977e-16      -1.068335e+00       1.197563e+00 
#> peptide_Idq2jTaC1y 
#>      -1.445323e+00 
#> 
#> Likelihood ratio test=10.07111 on 8 df, p=0.2600707, n=84
#> 
modelFunction$model_fun(nestProtein$data[[4]])
#> logistf::logistf(formula = self$formula, data = DFT, pl = FALSE, 
#>     control = do.call(logistf::logistf.control, self$control), 
#>     weights = Freq)
#> Model fitted by Penalized ML
#> Confidence intervals and p-values by Wald 
#> 
#> Coefficients:
#>        (Intercept)            group_B         group_Ctrl peptide_IdWcAw5ozd 
#>       8.375544e-03       8.360364e-11       9.365126e-01      -3.147690e-01 
#> peptide_IdgdnXrza3 peptide_IdxvlVt88v 
#>      -2.035803e-09      -6.303721e-01 
#> 
#> Likelihood ratio test=3.276542 on 5 df, p=0.657435, n=48
#> 
tmp <- strategy_rfit("Intensity ~ condition", model_name = "parallel design")
tmp$formula
#> Intensity ~ condition
#> <environment: 0x615a89949170>
modelFunction <- strategy_lmer("abundanceC ~ group_ + (1|peptide_Id)",
  model_name = "random_example")
modelFunction$formula
#> abundanceC ~ group_ + (1 | peptide_Id)
#> <environment: 0x615a8c51deb8>
tmp <- strategy_lm("Intensity ~ condition", model_name = "parallel design")
tmp$formula
#> Intensity ~ condition
#> <environment: 0x615a879220e0>
tmp$weights
#> NULL
tmp <- strategy_rlm("Intensity ~ condition", model_name = "parallel design")
tmp$formula
#> Intensity ~ condition
#> <environment: 0x615a8a862348>
strategy_rlm("Intensity ~ condition", control = list(maxit = 50))$control
#> $maxit
#> [1] 50
#> 
#> $acc
#> [1] 1e-04
#> 
```
