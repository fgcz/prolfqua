# Build limma model with vooma precision weights (proteomics)

Estimates observation-level precision weights from a mean-variance trend
(vooma) and fits a weighted least squares model via
[`lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html). For proteomics data
that is already log2-transformed.

## Usage

``` r
build_model_limma_voom(
  lfqdata,
  strategy,
  model_name = strategy$model_name,
  span = 0.5,
  plot = FALSE
)
```

## Arguments

- lfqdata:

  an [`LFQData`](https://fgcz.github.io/prolfqua/reference/LFQData.md)
  object (aggregated to protein level)

- strategy:

  output of
  [`strategy_limma`](https://fgcz.github.io/prolfqua/reference/strategy_limma.md)

- model_name:

  name of model (default from strategy)

- span:

  lowess smoother span for mean-variance trend (default 0.5)

- plot:

  logical; if TRUE, plot the mean-variance trend

## Value

a
[`ModelLimma`](https://fgcz.github.io/prolfqua/reference/ModelLimma.md)
object

## Details

The weights follow [`vooma`](https://rdrr.io/pkg/limma/man/vooma.html)
(Law et al. 2014,
[doi:10.1186/gb-2014-15-2-r29](https://doi.org/10.1186/gb-2014-15-2-r29)
):

1.  a preliminary `lmFit()` gives the residual standard deviation
    \\s_g\\ of each protein;

2.  a lowess trend \\f\\ (span `span`) of \\\sqrt{s_g}\\ against the
    mean log2 intensity is evaluated at each fitted value
    \\\hat{y}\_{gj}\\;

3.  since \\f^4\\ estimates the variance, the precision weight is
    \\w\_{gj} = 1 / f(\hat{y}\_{gj})^4\\, and the final fit is weighted
    least squares followed by `eBayes()`.

The data are already on the log2 scale, so there is no library-size
term. With `weights = NULL` and complete data the fit equals
[`limma::vooma()`](https://rdrr.io/pkg/limma/man/vooma.html) followed by
[`limma::lmFit()`](https://rdrr.io/pkg/limma/man/lmFit.html). In designs
with group factors only, all samples of a group share one weight: fold
changes equal those of unweighted limma, while standard errors and
p-values change.

External weights \\n\_{gj}\\ (`strategy$weights`, e.g. the `nr_children`
column, the number of peptides per protein and sample) go beyond
`vooma()`, which takes none. If a protein intensity averages \\n\_{gj}\\
peptide measurements with independent errors, its variance is
\\\sigma^2_g / n\_{gj}\\. The preliminary fit therefore uses \\n\_{gj}\\
as weights, so the trend is estimated from weighted residuals, and the
final fit uses \\n\_{gj} w\_{gj}\\. Only the relative size of
\\n\_{gj}\\ matters: a constant weight gives the same result as none.
`deqms_voom` uses no external weights, because DEqMS already moderates
by the peptide count.

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

## Examples

``` r
istar <- sim_lfq_data_protein_config(n_proteins = 50)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 5 unconfigured columns: average_prot_abundance, sd, mean, N, Replicate.
#> setup_analysis: added 82 rows with NA response for feature and sample combinations absent from the input.
lProt <- LFQData$new(istar$data, istar$config)
lProt$rename_response("transformedIntensity")

strat <- strategy_limma("transformedIntensity ~ group_")
mod <- build_model_limma_voom(lProt, strat)
#> fit (limma): 1 singular fits (e.g. VUnZ6o~1614): coefficients not estimable (Partial NA coefficients)
#> fit (limma): 1 limma warnings (e.g. Partial NA coefficients for 1 probe(s)): raised by limma while fitting all proteins; kept verbatim
#> fit (limma): 82 missing responses dropped (e.g. 0EfVhX~7161, 0m5WN4~3543, 7QuTub~5556, 7cbcrd~0495, 9VUkAq~8655, ...): NA response (missing observation) left out of the fit
#> See $provenance().
mod$get_coefficients()
#> # A tibble: 150 × 6
#>    protein_Id  factor      Estimate Std..Error t.value   Pr...t..
#>    <chr>       <chr>          <dbl>      <dbl>   <dbl>      <dbl>
#>  1 0EfVhX~7161 (Intercept)     20.3      0.465    43.6 2.54 e-147
#>  2 0m5WN4~3543 (Intercept)     20.5      0.451    45.5 3.83 e-153
#>  3 76k03k~9735 (Intercept)     20.2      0.472    42.7 9.11 e-145
#>  4 7QuTub~5556 (Intercept)     22.8      0.516    44.2 2.60 e-149
#>  5 7cbcrd~0495 (Intercept)     17.2      0.670    25.6 8.66 e- 84
#>  6 7soopj~3451 (Intercept)     26.3      0.490    53.8 1.70 e-176
#>  7 9VUkAq~8655 (Intercept)     22.2      0.499    44.6 2.48 e-150
#>  8 At886V~0359 (Intercept)     17.1      0.666    25.7 3.34 e- 84
#>  9 BEJI92~5483 (Intercept)     15.8      0.370    42.7 1.000e-144
#> 10 CGzoYe~1248 (Intercept)     18.2      0.537    34.0 1.84 e-115
#> # ℹ 140 more rows
```
