# Firth-penalized detection-count facade for nested input

Firth-penalized detection-count facade for nested input

Firth-penalized detection-count facade for nested input

## Value

An R6 class generator.

## Details

Same pipeline as
[`ContrastsBinomialNestedFacade`](https://wolski.github.io/prolfqua/reference/ContrastsBinomialNestedFacade.md):
child-feature detection is collapsed into detected and undetected counts
per parent and sample, so the dispersion and the residual degrees of
freedom come from the samples, not from the number of child features.
The fixed pseudo-count is replaced by Firth's bias reduction
([`StrategyBinomialFirth`](https://wolski.github.io/prolfqua/reference/StrategyBinomialFirth.md)),
which gives finite log odds ratios under complete separation, as
[`ContrastsFirthNestedFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthNestedFacade.md)
does. Unlike `firth_nested`, child identity is not modelled.

## Super classes

[`prolfqua::ContrastsInterface`](https://wolski.github.io/prolfqua/reference/ContrastsInterface.md)
-\>
[`prolfqua::ContrastsFacadeBase`](https://wolski.github.io/prolfqua/reference/ContrastsFacadeBase.md)
-\> `ContrastsBinomialFirthNestedFacade`

## Methods

### Public methods

- [`ContrastsBinomialFirthNestedFacade$new()`](#method-ContrastsBinomialFirthNestedFacade-new)

- [`ContrastsBinomialFirthNestedFacade$clone()`](#method-ContrastsBinomialFirthNestedFacade-clone)

Inherited methods

- [`prolfqua::ContrastsInterface$column_description()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-column_description)
- [`prolfqua::ContrastsInterface$contrast_summary_table()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-contrast_summary_table)
- [`prolfqua::ContrastsInterface$extra_artifacts()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-extra_artifacts)
- [`prolfqua::ContrastsInterface$filter_significant()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-filter_significant)
- [`prolfqua::ContrastsInterface$get_contrast_sides()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_contrast_sides)
- [`prolfqua::ContrastsInterface$get_ora()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_ora)
- [`prolfqua::ContrastsInterface$get_rank()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_rank)
- [`prolfqua::ContrastsFacadeBase$get_Plotter()`](https://wolski.github.io/prolfqua/html/ContrastsFacadeBase.html#method-ContrastsFacadeBase-get_Plotter)
- [`prolfqua::ContrastsFacadeBase$get_config()`](https://wolski.github.io/prolfqua/html/ContrastsFacadeBase.html#method-ContrastsFacadeBase-get_config)
- [`prolfqua::ContrastsFacadeBase$get_contrasts()`](https://wolski.github.io/prolfqua/html/ContrastsFacadeBase.html#method-ContrastsFacadeBase-get_contrasts)
- [`prolfqua::ContrastsFacadeBase$get_missing()`](https://wolski.github.io/prolfqua/html/ContrastsFacadeBase.html#method-ContrastsFacadeBase-get_missing)
- [`prolfqua::ContrastsFacadeBase$to_wide()`](https://wolski.github.io/prolfqua/html/ContrastsFacadeBase.html#method-ContrastsFacadeBase-to_wide)

------------------------------------------------------------------------

### Method `new()`

Fit the Firth-penalized nested detection-count model and its contrasts.

#### Usage

    ContrastsBinomialFirthNestedFacade$new(
      lfqdata,
      modelstr,
      contrasts,
      binomial_bound = TRUE,
      ...
    )

#### Arguments

- `lfqdata`:

  nested
  [`LFQData`](https://wolski.github.io/prolfqua/reference/LFQData.md)

- `modelstr`:

  right-hand-side model formula

- `contrasts`:

  named contrast expressions

- `binomial_bound`:

  bound posterior dispersion below by one

- `...`:

  passed to
  [`strategy_binomial_firth`](https://wolski.github.io/prolfqua/reference/strategy_binomial_firth.md)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsBinomialFirthNestedFacade$clone(deep = FALSE)

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
lfqdata <- LFQData$new(istar$data, istar$config)
contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
facade <- ContrastsBinomialFirthNestedFacade$new(lfqdata, "~ group_", contrasts)
#> completing cases
head(facade$get_contrasts())
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`
#> # A tibble: 6 × 16
#>   modelName             estimate_type protein_Id  contrast       diff avgAbd
#>   <chr>                 <chr>         <chr>       <chr>         <dbl>  <dbl>
#> 1 binomial_firth_nested observed      0GRprF~7339 A_vs_Ctrl -8.47e- 1  0.424
#> 2 binomial_firth_nested observed      4JK499~3111 A_vs_Ctrl  5.04e- 1 -0.704
#> 3 binomial_firth_nested observed      7IZdVV~6818 A_vs_Ctrl -5.73e- 1  1.32 
#> 4 binomial_firth_nested observed      AZPG26~9461 A_vs_Ctrl  1.21e-16  2.20 
#> 5 binomial_firth_nested observed      AoNKbb~3497 A_vs_Ctrl -1.18e+ 0  2.63 
#> 6 binomial_firth_nested observed      CibL2O~2149 A_vs_Ctrl -1.44e+ 0 -0.718
#> # ℹ 10 more variables: std.error.unmoderated <dbl>, df.unmoderated <int>,
#> #   std.error <dbl>, statistic <dbl>, df <dbl>, p.value <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, sigma <dbl>, FDR <dbl>
```
