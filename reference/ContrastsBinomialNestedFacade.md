# Quasibinomial detection-count facade for nested input

Quasibinomial detection-count facade for nested input

Quasibinomial detection-count facade for nested input

## Value

An R6 class generator.

## Details

Completes and encodes child-feature detection using the same preparation
as
[`ContrastsFirthNestedFacade`](https://wolski.github.io/prolfqua/reference/ContrastsFirthNestedFacade.md),
then collapses the binary rows into detected and undetected counts per
parent and sample. The resulting quasibinomial model reports
parent-level log odds ratios in `diff`; `avgAbd` is the average linear
predictor on the log-odds scale.

Child features are treated as exchangeable binomial trials. The
symmetric pseudo-count stabilizes complete separation but is not
equivalent to Firth's bias-reducing penalty. Empirical-Bayes dispersion
moderation uses
[`ContrastsModerated`](https://wolski.github.io/prolfqua/reference/ContrastsModerated.md);
by default posterior dispersion is bounded below by one.

## Super classes

[`prolfqua::ContrastsInterface`](https://wolski.github.io/prolfqua/reference/ContrastsInterface.md)
-\>
[`prolfqua::ContrastsFacadeBase`](https://wolski.github.io/prolfqua/reference/ContrastsFacadeBase.md)
-\> `ContrastsBinomialNestedFacade`

## Methods

### Public methods

- [`ContrastsBinomialNestedFacade$new()`](#method-ContrastsBinomialNestedFacade-new)

- [`ContrastsBinomialNestedFacade$clone()`](#method-ContrastsBinomialNestedFacade-clone)

Inherited methods

- [`prolfqua::ContrastsInterface$column_description()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-column_description)
- [`prolfqua::ContrastsInterface$contrast_summary_table()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-contrast_summary_table)
- [`prolfqua::ContrastsInterface$extra_artifacts()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-extra_artifacts)
- [`prolfqua::ContrastsInterface$filter_significant()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-filter_significant)
- [`prolfqua::ContrastsInterface$get_config()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_config)
- [`prolfqua::ContrastsInterface$get_contrast_sides()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_contrast_sides)
- [`prolfqua::ContrastsInterface$get_ora()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_ora)
- [`prolfqua::ContrastsInterface$get_rank()`](https://wolski.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_rank)
- [`prolfqua::ContrastsFacadeBase$get_Plotter()`](https://wolski.github.io/prolfqua/html/ContrastsFacadeBase.html#method-ContrastsFacadeBase-get_Plotter)
- [`prolfqua::ContrastsFacadeBase$get_contrasts()`](https://wolski.github.io/prolfqua/html/ContrastsFacadeBase.html#method-ContrastsFacadeBase-get_contrasts)
- [`prolfqua::ContrastsFacadeBase$get_missing()`](https://wolski.github.io/prolfqua/html/ContrastsFacadeBase.html#method-ContrastsFacadeBase-get_missing)
- [`prolfqua::ContrastsFacadeBase$to_wide()`](https://wolski.github.io/prolfqua/html/ContrastsFacadeBase.html#method-ContrastsFacadeBase-to_wide)

------------------------------------------------------------------------

### Method `new()`

Fit the nested detection-count model and its contrasts.

#### Usage

    ContrastsBinomialNestedFacade$new(
      lfqdata,
      modelstr,
      contrasts,
      prior_count = 0.1,
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

- `prior_count`:

  non-negative symmetric pseudo-count

- `binomial_bound`:

  bound posterior dispersion below by one

- `...`:

  passed to
  [`strategy_binomial`](https://wolski.github.io/prolfqua/reference/strategy_binomial.md)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsBinomialNestedFacade$clone(deep = FALSE)

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
facade <- ContrastsBinomialNestedFacade$new(lfqdata, "~ group_", contrasts)
#> completing cases
head(facade$get_contrasts())
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`
#> # A tibble: 6 × 14
#>   modelName       estimate_type protein_Id  contrast       diff std.error avgAbd
#>   <chr>           <chr>         <chr>       <chr>         <dbl>     <dbl>  <dbl>
#> 1 binomial_nested observed      0GRprF~7339 A_vs_Ctrl -8.87e- 1  1.39e+ 0  0.444
#> 2 binomial_nested observed      4JK499~3111 A_vs_Ctrl  5.18e- 1  1.03e+ 0 -0.722
#> 3 binomial_nested observed      7IZdVV~6818 A_vs_Ctrl -5.85e- 1  7.54e- 1  1.34 
#> 4 binomial_nested observed      AZPG26~9461 A_vs_Ctrl  9.53e-16  1.41e-12  2.40 
#> 5 binomial_nested observed      AoNKbb~3497 A_vs_Ctrl -1.34e+ 0  1.60e+ 0  2.77 
#> 6 binomial_nested observed      CibL2O~2149 A_vs_Ctrl -1.47e+ 0  6.43e- 1 -0.733
#> # ℹ 7 more variables: statistic <dbl>, df <dbl>, p.value <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, sigma <dbl>, FDR <dbl>
```
