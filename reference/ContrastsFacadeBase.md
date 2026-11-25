# Base class for contrast analysis facades

The extension point for new modelling backends. Every method of
[`build_contrast_analysis`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md)
is a facade class that inherits from `ContrastsFacadeBase`; a package
adds a method by writing such a class and registering it with
[`register_facade`](https://fgcz.github.io/prolfqua/reference/register_facade.md)
(usually from its `.onLoad()`). The base class provides the uniform API:
`get_contrasts()` (stamped with the method name), `get_missing()`,
`get_Plotter()`, `to_wide()` and `get_config()`.

A subclass implements `initialize(lfqdata, modelstr, contrasts, ...)`:

- call the private `.setup(lfqdata, contrasts, key, modelstr)` first; it
  checks the input shape, the contrasts and the design, stores a copy of
  the input, and returns the full model formula string

- set `self$model` to the fitted model and `self$contrast` to the
  contrast object the facade delegates to (e.g.
  [`ContrastsModerated`](https://fgcz.github.io/prolfqua/reference/ContrastsModerated.md))

- set the private field `.needs` to `"nested"` if the facade takes
  child-level input (e.g. peptides) and reports parent-level contrasts;
  the default is `"aggregated"`

The limma family, `lm`, `deqms` and `deqms_voom` additionally set
`self$.drop_na_diff <- TRUE` to drop rows whose fold change could not be
estimated; `get_missing()` lists them.

## Value

An R6 class generator.

## See also

[`register_facade`](https://fgcz.github.io/prolfqua/reference/register_facade.md),
[`list_facades`](https://fgcz.github.io/prolfqua/reference/list_facades.md)

Other modelling:
[`AnovaExtractor`](https://fgcz.github.io/prolfqua/reference/AnovaExtractor.md),
[`Contrasts`](https://fgcz.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsDEqMSVoomFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsDEqMSVoomFacade.md),
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
[`strategy_logistf()`](https://fgcz.github.io/prolfqua/reference/strategy.md),
[`vcov.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/vcov.rfit_prolfqua.md)

## Super class

[`ContrastsInterface`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.md)
-\> `ContrastsFacadeBase`

## Public fields

- `model`:

  fitted model object

- `contrast`:

  inner contrast object the facade delegates to

- `.lfqdata`:

  snapshot of the input LFQData (used by get_missing)

- `.contrast_names`:

  names of the requested contrasts

- `facade_name`:

  registry key stamped into the modelName column

- `.drop_na_diff`:

  drop rows with NA diff (limma family); default FALSE

## Methods

### Public methods

- [`ContrastsFacadeBase$get_config()`](#method-ContrastsFacadeBase-get_config)

- [`ContrastsFacadeBase$get_contrasts()`](#method-ContrastsFacadeBase-get_contrasts)

- [`ContrastsFacadeBase$get_missing()`](#method-ContrastsFacadeBase-get_missing)

- [`ContrastsFacadeBase$get_Plotter()`](#method-ContrastsFacadeBase-get_Plotter)

- [`ContrastsFacadeBase$to_wide()`](#method-ContrastsFacadeBase-to_wide)

- [`ContrastsFacadeBase$provenance()`](#method-ContrastsFacadeBase-provenance)

- [`ContrastsFacadeBase$.recorded_provenance()`](#method-ContrastsFacadeBase-.recorded_provenance)

- [`ContrastsFacadeBase$clone()`](#method-ContrastsFacadeBase-clone)

Inherited methods

- [`ContrastsInterface$column_description()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-column_description)
- [`ContrastsInterface$contrast_summary_table()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-contrast_summary_table)
- [`ContrastsInterface$extra_artifacts()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-extra_artifacts)
- [`ContrastsInterface$filter_significant()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-filter_significant)
- [`ContrastsInterface$get_contrast_sides()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_contrast_sides)
- [`ContrastsInterface$get_ora()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_ora)
- [`ContrastsInterface$get_rank()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_rank)

------------------------------------------------------------------------

### `ContrastsFacadeBase$get_config()`

get the
[`ContrastConfiguration`](https://fgcz.github.io/prolfqua/reference/ContrastConfiguration.md);
defaults to the inner contrast object's configuration

#### Usage

    ContrastsFacadeBase$get_config()

------------------------------------------------------------------------

### `ContrastsFacadeBase$get_contrasts()`

get contrast results, stamped with the facade key

#### Usage

    ContrastsFacadeBase$get_contrasts(...)

#### Arguments

- `...`:

  passed to the inner contrast object's get_contrasts

------------------------------------------------------------------------

### `ContrastsFacadeBase$get_missing()`

get subject x contrast pairs that could not be estimated

#### Usage

    ContrastsFacadeBase$get_missing()

------------------------------------------------------------------------

### `ContrastsFacadeBase$get_Plotter()`

get
[`ContrastsPlotter`](https://fgcz.github.io/prolfqua/reference/ContrastsPlotter.md)

#### Usage

    ContrastsFacadeBase$get_Plotter(...)

#### Arguments

- `...`:

  passed to the inner contrast object's get_Plotter

------------------------------------------------------------------------

### `ContrastsFacadeBase$to_wide()`

convert results to wide format

#### Usage

    ContrastsFacadeBase$to_wide(...)

#### Arguments

- `...`:

  passed to the inner contrast object's to_wide

------------------------------------------------------------------------

### `ContrastsFacadeBase$provenance()`

problems and dropped items of this analysis, read from the input data,
the aggregated data (`hurdle_nested`, `limpa_nested`) and the model(s),
plus the subjects absent from the contrasts with the reason

#### Usage

    ContrastsFacadeBase$provenance()

#### Returns

tibble with columns source, step, event, n, ids (list), reason, failure,
response

------------------------------------------------------------------------

### `ContrastsFacadeBase$.recorded_provenance()`

internal: the entries recorded in the input data, the aggregated data
and the model(s), without evaluating contrasts (used by
[`print()`](https://rdrr.io/r/base/print.html))

#### Usage

    ContrastsFacadeBase$.recorded_provenance()

#### Returns

tibble as `provenance()`, without the absent subjects

------------------------------------------------------------------------

### `ContrastsFacadeBase$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsFacadeBase$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# a facade: unweighted linear model per protein, moderated contrasts
ContrastsUnweightedLMFacade <- R6::R6Class(
  "ContrastsUnweightedLMFacade",
  inherit = ContrastsFacadeBase,
  public = list(
    initialize = function(lfqdata, modelstr, contrasts) {
      full_formula <- private$.setup(lfqdata, contrasts, "lm_unweighted", modelstr)
      self$model <- build_model(lfqdata, strategy_lm(full_formula))
      self$contrast <- ContrastsModerated$new(Contrasts$new(self$model, contrasts))
    }
  )
)
istar <- sim_lfq_data_protein_config(n_proteins = 10)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 5 unconfigured columns: average_prot_abundance, sd, mean, N, Replicate.
#> setup_analysis: added 16 rows with NA response for feature and sample combinations absent from the input.
lfqdata <- LFQData$new(istar$data, istar$config)
facade <- ContrastsUnweightedLMFacade$new(lfqdata, "~ group_", c(AvsCtrl = "group_A - group_Ctrl"))
#> fit (Model): 16 missing responses dropped (e.g. 0EfVhX~0087, 7cbcrd~5725, 9VUkAq~4703, BEJI92~5282, DoWup2~5896, ...): NA response (missing observation) left out of the fit
#> See $provenance().
head(facade$get_contrasts())
#> # A tibble: 6 × 16
#>   modelName     estimate_type protein_Id  contrast    diff avgAbd
#>   <chr>         <chr>         <chr>       <chr>      <dbl>  <dbl>
#> 1 lm_unweighted observed      0EfVhX~0087 AvsCtrl  -2.62     21.1
#> 2 lm_unweighted observed      7cbcrd~5725 AvsCtrl   2.80     20.7
#> 3 lm_unweighted observed      9VUkAq~4703 AvsCtrl   1.67     20.3
#> 4 lm_unweighted observed      BEJI92~5282 AvsCtrl   0.424    21.0
#> 5 lm_unweighted observed      CGzoYe~2147 AvsCtrl  -0.598    30.8
#> 6 lm_unweighted observed      Fl4JiV~8625 AvsCtrl  -0.0494   21.3
#> # ℹ 10 more variables: std.error.unmoderated <dbl>, df.unmoderated <int>,
#> #   std.error <dbl>, statistic <dbl>, df <dbl>, p.value <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, sigma <dbl>, FDR <dbl>
# registered from a package's .onLoad(), build_contrast_analysis() finds the method by name:
# register_facade("lm_unweighted", class = "ContrastsUnweightedLMFacade", needs = "aggregated",
#                 package = "mypackage")
```
