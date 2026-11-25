# Compare the contrasts of several methods on their subject and contrast pairs

Aligns the contrast results of several methods, e.g. facades of
[`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md)
or
[ContrastsTable](https://fgcz.github.io/prolfqua/reference/ContrastsTable.md)s
holding another tool's results, on the subject (e.g. protein) and the
contrast, and reports which pairs each method estimates and where their
inferences contradict each other. The columns of each result are read
through its
[ContrastConfiguration](https://fgcz.github.io/prolfqua/reference/ContrastConfiguration.md),
so tables with other column names can be compared too.

The status of a pair in a method is `not estimated` (no row or no
effect), `not tested` (an effect but no FDR), `not significant`, or `up`
/ `down` (FDR below `fdr_threshold`, by the sign of the effect).

## Value

An R6 class generator.

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
[`ContrastsROPECAImputeNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAImputeNestedFacade.md),
[`ContrastsROPECANestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECANestedFacade.md),
[`ContrastsROPECAV2`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAV2.md),
[`ContrastsROPECAV2ImputeNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAV2ImputeNestedFacade.md),
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
[`strategy_logistf()`](https://fgcz.github.io/prolfqua/reference/strategy.md),
[`vcov.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/vcov.rfit_prolfqua.md)

## Public fields

- `subject_id`:

  columns identifying the subjects, shared by all methods

- `fdr_threshold`:

  FDR below which a contrast is significant

## Methods

### Public methods

- [`ContrastsComparison$new()`](#method-ContrastsComparison-initialize)

- [`ContrastsComparison$aligned()`](#method-ContrastsComparison-aligned)

- [`ContrastsComparison$pairs()`](#method-ContrastsComparison-pairs)

- [`ContrastsComparison$coverage()`](#method-ContrastsComparison-coverage)

- [`ContrastsComparison$conflicts()`](#method-ContrastsComparison-conflicts)

- [`ContrastsComparison$upset()`](#method-ContrastsComparison-upset)

- [`ContrastsComparison$clone()`](#method-ContrastsComparison-clone)

------------------------------------------------------------------------

### `ContrastsComparison$new()`

align the contrast results of the methods

#### Usage

    ContrastsComparison$new(contrasts, fdr_threshold = 0.05)

#### Arguments

- `contrasts`:

  named list of contrast objects (facades,
  [ContrastsInterface](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.md)
  objects,
  [ContrastsTable](https://fgcz.github.io/prolfqua/reference/ContrastsTable.md)s);
  the names label the methods

- `fdr_threshold`:

  FDR below which a contrast is significant

------------------------------------------------------------------------

### `ContrastsComparison$aligned()`

one row per subject, contrast and method over all pairs any method
reports: `effect`, `statistic`, `p.value`, `FDR`, `estimate_type` (NA
for tables without it) and `status`

#### Usage

    ContrastsComparison$aligned()

#### Returns

data.frame

------------------------------------------------------------------------

### `ContrastsComparison$pairs()`

the pairs of one method with all the given columns defined, e.g.
`c("effect", "p.value")` for the pairs with a fold change and a p-value

#### Usage

    ContrastsComparison$pairs(method, columns = "effect")

#### Arguments

- `method`:

  method label, one of the names of `contrasts`

- `columns`:

  one or more of `"effect"`, `"statistic"`, `"p.value"`

#### Returns

data.frame with the subject columns and `contrast`

------------------------------------------------------------------------

### `ContrastsComparison$coverage()`

number of pairs per set of methods that estimate them (`none`: reported,
but without an effect)

#### Usage

    ContrastsComparison$coverage()

#### Returns

data.frame with columns `estimated_by` and `pairs`

------------------------------------------------------------------------

### `ContrastsComparison$conflicts()`

per pair, the number of methods by status, and where the methods
contradict each other: `sign_conflict` (up in one method, down in
another), `significance_conflict` (significant in one, tested but not
significant in another), `coverage_conflict` (significant in one, not
estimated by another)

#### Usage

    ContrastsComparison$conflicts()

#### Returns

data.frame, one row per pair

------------------------------------------------------------------------

### `ContrastsComparison$upset()`

UpSet plot of the pairs each method estimates

#### Usage

    ContrastsComparison$upset()

#### Returns

the plot of
[`UpSetR::upset()`](https://rdrr.io/pkg/UpSetR/man/upset.html)

------------------------------------------------------------------------

### `ContrastsComparison$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsComparison$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
sim <- sim_lfq_data_protein_config(n_proteins = 20, weight_missing = 0.5, seed = 3)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 5 unconfigured columns: average_prot_abundance, sd, mean, N, Replicate.
#> setup_analysis: added 62 rows with NA response for feature and sample combinations absent from the input.
lfqdata <- LFQData$new(sim$data, sim$config)
contrasts <- c(A_vs_Ctrl = "group_A - group_Ctrl", B_vs_Ctrl = "group_B - group_Ctrl")
methods <- c("lm", "limma", "lm_impute")
fits <- lapply(stats::setNames(methods, methods), function(method) {
  build_contrast_analysis(lfqdata, "~ group_", contrasts, method = method)
})
#> fit (Model): 62 missing responses dropped (e.g. 0GRprF~7339, 4JK499~3111, 7IZdVV~6818, AZPG26~9461, AoNKbb~3497, ...): NA response (missing observation) left out of the fit
#> See $provenance().
#> fit (limma): 3 singular fits (e.g. 4JK499~3111, CibL2O~2149, EnbuYk~7811): coefficients not estimable (Partial NA coefficients)
#> fit (limma): 62 missing responses dropped (e.g. 0GRprF~7339, 4JK499~3111, 7IZdVV~6818, AZPG26~9461, AoNKbb~3497, ...): NA response (missing observation) left out of the fit
#> See $provenance().
#> fit (ModelImputed): 3 refit at the limit of detection (e.g. 4JK499~3111, CibL2O~2149, EnbuYk~7811): fit failed, singular or missing a group; missing responses imputed at the LOD
#> fit (ModelImputed): 39 missing responses dropped (e.g. 0GRprF~7339, 7IZdVV~6818, AZPG26~9461, AoNKbb~3497, JnEDsz~4232, ...): NA response (missing observation) left out of the fit
#> See $provenance().
comparison <- ContrastsComparison$new(fits)
head(comparison$aligned())
#> # A tibble: 6 × 9
#>   protein_Id  contrast  method    effect statistic p.value     FDR estimate_type
#>   <chr>       <chr>     <chr>      <dbl>     <dbl>   <dbl>   <dbl> <chr>        
#> 1 0GRprF~7339 A_vs_Ctrl lm        -0.220    -0.240 8.18e-1 0.842   observed     
#> 2 0GRprF~7339 A_vs_Ctrl limma     -0.220    -0.241 8.17e-1 0.841   observed     
#> 3 0GRprF~7339 A_vs_Ctrl lm_impute -0.220    -0.225 8.28e-1 0.993   observed     
#> 4 7IZdVV~6818 A_vs_Ctrl lm         2.74      5.18  3.69e-4 0.00628 observed     
#> 5 7IZdVV~6818 A_vs_Ctrl limma      2.74      5.20  3.44e-4 0.00585 observed     
#> 6 7IZdVV~6818 A_vs_Ctrl lm_impute  2.74      5.16  3.24e-4 0.00648 observed     
#> # ℹ 1 more variable: status <chr>
comparison$coverage()
#> # A tibble: 2 × 2
#>   estimated_by           pairs
#>   <chr>                  <int>
#> 1 lm & limma & lm_impute    35
#> 2 lm_impute                  5
conflicts <- comparison$conflicts()
conflicts[conflicts$coverage_conflict | conflicts$significance_conflict, ]
#> # A tibble: 1 × 10
#>   protein_Id  contrast     up  down not_significant not_tested not_estimated
#>   <chr>       <chr>     <int> <int>           <int>      <int>         <int>
#> 1 EnbuYk~7811 B_vs_Ctrl     2     0               1          0             0
#> # ℹ 3 more variables: sign_conflict <lgl>, significance_conflict <lgl>,
#> #   coverage_conflict <lgl>
head(comparison$pairs("lm", columns = c("effect", "p.value")))
#> # A tibble: 6 × 2
#>   protein_Id  contrast 
#>   <chr>       <chr>    
#> 1 0GRprF~7339 A_vs_Ctrl
#> 2 7IZdVV~6818 A_vs_Ctrl
#> 3 AZPG26~9461 A_vs_Ctrl
#> 4 AoNKbb~3497 A_vs_Ctrl
#> 5 JnEDsz~4232 A_vs_Ctrl
#> 6 KIBglg~8285 A_vs_Ctrl
if (requireNamespace("UpSetR", quietly = TRUE)) comparison$upset()
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the UpSetR package.
#>   Please report the issue at <https://github.com/hms-dbmi/UpSetR/issues>.
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the UpSetR package.
#>   Please report the issue at <https://github.com/hms-dbmi/UpSetR/issues>.
#> Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
#> ℹ Please use the `linewidth` argument instead.
#> ℹ The deprecated feature was likely used in the UpSetR package.
#>   Please report the issue at <https://github.com/hms-dbmi/UpSetR/issues>.
```
