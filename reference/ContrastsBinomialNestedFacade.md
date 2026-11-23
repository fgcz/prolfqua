# Quasibinomial detection-count facade for nested input

Completes the child-feature rows, encodes detection as 0/1, then
collapses the binary rows into detected and undetected counts per parent
and sample. The resulting quasibinomial model reports parent-level log
odds ratios in `diff`; `avgAbd` is the average linear predictor on the
log-odds scale.

Child features are treated as exchangeable binomial trials. The
symmetric pseudo-count stabilizes complete separation but is not
equivalent to Firth's bias-reducing penalty. Empirical-Bayes dispersion
moderation uses
[`ContrastsModerated`](https://fgcz.github.io/prolfqua/reference/ContrastsModerated.md);
by default posterior dispersion is bounded below by one.

## Value

An R6 class generator.

## Super classes

[`ContrastsInterface`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.md)
-\>
[`ContrastsFacadeBase`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.md)
-\> `ContrastsBinomialNestedFacade`

## Methods

### Public methods

- [`ContrastsBinomialNestedFacade$new()`](#method-ContrastsBinomialNestedFacade-initialize)

- [`ContrastsBinomialNestedFacade$clone()`](#method-ContrastsBinomialNestedFacade-clone)

Inherited methods

- [`ContrastsInterface$column_description()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-column_description)
- [`ContrastsInterface$contrast_summary_table()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-contrast_summary_table)
- [`ContrastsInterface$extra_artifacts()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-extra_artifacts)
- [`ContrastsInterface$filter_significant()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-filter_significant)
- [`ContrastsInterface$get_contrast_sides()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_contrast_sides)
- [`ContrastsInterface$get_ora()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_ora)
- [`ContrastsInterface$get_rank()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_rank)
- [`ContrastsFacadeBase$.recorded_provenance()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-.recorded_provenance)
- [`ContrastsFacadeBase$get_Plotter()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-get_Plotter)
- [`ContrastsFacadeBase$get_config()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-get_config)
- [`ContrastsFacadeBase$get_contrasts()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-get_contrasts)
- [`ContrastsFacadeBase$get_missing()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-get_missing)
- [`ContrastsFacadeBase$provenance()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-provenance)
- [`ContrastsFacadeBase$to_wide()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-to_wide)

------------------------------------------------------------------------

### `ContrastsBinomialNestedFacade$new()`

Fit the nested detection-count model and its contrasts.

#### Usage

    ContrastsBinomialNestedFacade$new(
      lfqdata,
      modelstr,
      contrasts,
      prior_count = 0.1,
      binomial_bound = TRUE,
      control = list()
    )

#### Arguments

- `lfqdata`:

  nested
  [`LFQData`](https://fgcz.github.io/prolfqua/reference/LFQData.md)

- `modelstr`:

  right-hand-side model formula

- `contrasts`:

  named contrast expressions

- `prior_count`:

  non-negative symmetric pseudo-count

- `binomial_bound`:

  bound posterior dispersion below by one

- `control`:

  named list of convergence controls, see
  [`strategy_binomial`](https://fgcz.github.io/prolfqua/reference/strategy_binomial.md)

------------------------------------------------------------------------

### `ContrastsBinomialNestedFacade$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsBinomialNestedFacade$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- sim_lfq_data_peptide_config(n_proteins = 20, weight_missing = 0.5, seed = 3)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 205 rows with NA response for feature and sample combinations absent from the input.
lfqdata <- LFQData$new(istar$data, istar$config)
contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
facade <- ContrastsBinomialNestedFacade$new(lfqdata, "~ group_", contrasts)
head(facade$get_contrasts())
#> # A tibble: 6 × 16
#>   modelName       estimate_type protein_Id  contrast       diff avgAbd
#>   <chr>           <chr>         <chr>       <chr>         <dbl>  <dbl>
#> 1 binomial_nested observed      0GRprF~7339 A_vs_Ctrl -8.87e- 1  0.444
#> 2 binomial_nested observed      4JK499~3111 A_vs_Ctrl  5.18e- 1 -0.722
#> 3 binomial_nested observed      7IZdVV~6818 A_vs_Ctrl -5.85e- 1  1.34 
#> 4 binomial_nested observed      AZPG26~9461 A_vs_Ctrl -7.07e-16  2.40 
#> 5 binomial_nested observed      AoNKbb~3497 A_vs_Ctrl -1.34e+ 0  2.77 
#> 6 binomial_nested observed      CibL2O~2149 A_vs_Ctrl -1.47e+ 0 -0.733
#> # ℹ 10 more variables: std.error.unmoderated <dbl>, df.unmoderated <int>,
#> #   std.error <dbl>, statistic <dbl>, df <dbl>, p.value <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, sigma <dbl>, FDR <dbl>
```
