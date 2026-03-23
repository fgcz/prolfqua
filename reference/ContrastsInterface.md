# Base class for all Contrasts classes

Base class for all Contrasts classes

Base class for all Contrasts classes

## Methods

### Public methods

- [`ContrastsInterface$get_contrast_sides()`](#method-ContrastsInterface-get_contrast_sides)

- [`ContrastsInterface$get_contrasts()`](#method-ContrastsInterface-get_contrasts)

- [`ContrastsInterface$get_Plotter()`](#method-ContrastsInterface-get_Plotter)

- [`ContrastsInterface$to_wide()`](#method-ContrastsInterface-to_wide)

- [`ContrastsInterface$get_missing()`](#method-ContrastsInterface-get_missing)

- [`ContrastsInterface$column_description()`](#method-ContrastsInterface-column_description)

- [`ContrastsInterface$clone()`](#method-ContrastsInterface-clone)

------------------------------------------------------------------------

### Method `get_contrast_sides()`

get table with sides of the contrast

#### Usage

    ContrastsInterface$get_contrast_sides()

------------------------------------------------------------------------

### Method `get_contrasts()`

get table with contrast results (similar to limma topTable function)

#### Usage

    ContrastsInterface$get_contrasts()

------------------------------------------------------------------------

### Method `get_Plotter()`

initialize plotter

#### Usage

    ContrastsInterface$get_Plotter()

------------------------------------------------------------------------

### Method `to_wide()`

create wide representation of data.

#### Usage

    ContrastsInterface$to_wide()

------------------------------------------------------------------------

### Method `get_missing()`

get protein × contrast pairs that could not be estimated. Returns a
data.frame with hierarchy columns and contrast, or 0 rows.

#### Usage

    ContrastsInterface$get_missing()

------------------------------------------------------------------------

### Method `column_description()`

column description

#### Usage

    ContrastsInterface$column_description()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsInterface$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
int <- ContrastsInterface$new()
testthat::expect_error(int$get_contrast_sides())
testthat::expect_error(int$get_contrasts())
testthat::expect_error(int$get_missing())
testthat::expect_error(int$get_Plotter())
testthat::expect_error(int$to_wide())
int$column_description()
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
```
