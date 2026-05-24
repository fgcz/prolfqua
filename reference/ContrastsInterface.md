# Base class for all Contrasts classes

Base class for all Contrasts classes

Base class for all Contrasts classes

## Value

An R6 class generator.

## Public fields

- `config`:

  [`ContrastConfiguration`](https://wolski.github.io/prolfqua/reference/ContrastConfiguration.md)
  describing the backend's column-role mapping and behaviour flags.
  Subclasses should populate this in `initialize()` so generic consumers
  can resolve column names without hard-coding them.

## Methods

### Public methods

- [`ContrastsInterface$get_config()`](#method-ContrastsInterface-get_config)

- [`ContrastsInterface$get_contrast_sides()`](#method-ContrastsInterface-get_contrast_sides)

- [`ContrastsInterface$get_contrasts()`](#method-ContrastsInterface-get_contrasts)

- [`ContrastsInterface$get_Plotter()`](#method-ContrastsInterface-get_Plotter)

- [`ContrastsInterface$to_wide()`](#method-ContrastsInterface-to_wide)

- [`ContrastsInterface$get_missing()`](#method-ContrastsInterface-get_missing)

- [`ContrastsInterface$get_rank()`](#method-ContrastsInterface-get_rank)

- [`ContrastsInterface$get_ora()`](#method-ContrastsInterface-get_ora)

- [`ContrastsInterface$filter_significant()`](#method-ContrastsInterface-filter_significant)

- [`ContrastsInterface$contrast_summary_table()`](#method-ContrastsInterface-contrast_summary_table)

- [`ContrastsInterface$extra_artifacts()`](#method-ContrastsInterface-extra_artifacts)

- [`ContrastsInterface$column_description()`](#method-ContrastsInterface-column_description)

- [`ContrastsInterface$clone()`](#method-ContrastsInterface-clone)

------------------------------------------------------------------------

### Method `get_config()`

Return the
[`ContrastConfiguration`](https://wolski.github.io/prolfqua/reference/ContrastConfiguration.md),
constructing an LM-flavoured default on demand if a subclass has not set
one.

#### Usage

    ContrastsInterface$get_config()

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

### Method `get_rank()`

get rank-list input table for enrichment tools.

If the backend's
[`ContrastConfiguration`](https://wolski.github.io/prolfqua/reference/ContrastConfiguration.md)
has a p-value column the default rank score is signed
`sign(effect) * -log10(p.value)`. Otherwise the default is the effect
column itself.

#### Usage

    ContrastsInterface$get_rank(score = NULL)

#### Arguments

- `score`:

  column name to use as rank score, or `NULL` to pick from the config

#### Returns

data.frame with subject id, contrast, and score columns.

------------------------------------------------------------------------

### Method `get_ora()`

get ORA input table — features passing the FDR and absolute effect-size
threshold in the requested direction.

#### Usage

    ContrastsInterface$get_ora(up = TRUE, FDR_threshold = 0.05, diff_threshold = 1)

#### Arguments

- `up`:

  if TRUE return positive effects, otherwise negative effects

- `FDR_threshold`:

  false discovery rate threshold

- `diff_threshold`:

  absolute effect-size threshold

#### Returns

filtered contrast data.frame for ORA input.

------------------------------------------------------------------------

### Method `filter_significant()`

Filter contrast rows that pass the FDR and effect-size threshold.

By default this is symmetric on the effect column
(`|effect| > diff_threshold`). For backends whose
[`ContrastConfiguration`](https://wolski.github.io/prolfqua/reference/ContrastConfiguration.md)
has `significance_directional = TRUE` (e.g. SAINTexpress, where only
positive enrichment is biologically meaningful), the filter is one-sided
(`effect > diff_threshold`).

#### Usage

    ContrastsInterface$filter_significant(FDR_threshold = 0.05, diff_threshold = 1)

#### Arguments

- `FDR_threshold`:

  false discovery rate threshold

- `diff_threshold`:

  effect-size threshold

#### Returns

filtered contrast data.frame.

------------------------------------------------------------------------

### Method `contrast_summary_table()`

Per-subject summary table for downstream report grobs. Returns a small
data.frame with canonical column names `contrast`, `effect`, `score`,
`fdr` regardless of backend.

#### Usage

    ContrastsInterface$contrast_summary_table(rounded = TRUE)

#### Arguments

- `rounded`:

  if TRUE round numeric columns to 3 significant digits (matches the
  prolfquapp report convention)

------------------------------------------------------------------------

### Method `extra_artifacts()`

Backend-specific extra artifacts that a report should carry along with
the contrast results (e.g. SAINT input tables). Default returns an empty
named list.

#### Usage

    ContrastsInterface$extra_artifacts()

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
# get_rank / get_ora / filter_significant call get_contrasts() internally,
# which is not implemented on the bare interface, so they surface that error.
testthat::expect_error(int$get_rank())
testthat::expect_error(int$get_ora())
testthat::expect_error(int$filter_significant())
identical(int$extra_artifacts(), list())
#> [1] TRUE
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
