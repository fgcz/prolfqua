# Decorates LFQData with methods to compute statistics of interactions

Decorates LFQData with methods to compute statistics of interactions

Decorates LFQData with methods to compute statistics of interactions

## Details

compute stdv, mean and CV per peptide or protein and condition.

## See also

Other LFQData:
[`LFQData`](https://wolski.github.io/prolfqua/reference/LFQData.md),
[`LFQDataAggregator`](https://wolski.github.io/prolfqua/reference/LFQDataAggregator.md),
[`LFQDataImp`](https://wolski.github.io/prolfqua/reference/LFQDataImp.md),
[`LFQDataPlotter`](https://wolski.github.io/prolfqua/reference/LFQDataPlotter.md),
[`LFQDataSummariser`](https://wolski.github.io/prolfqua/reference/LFQDataSummariser.md),
[`LFQDataToSummarizedExperiment()`](https://wolski.github.io/prolfqua/reference/LFQDataToSummarizedExperiment.md)

## Public fields

- `lfq`:

  LFQData

- `stat`:

  either CV or sd (if is_transformed)

- `statsdf`:

  frame with statistics.

## Methods

### Public methods

- [`LFQDataStats$new()`](#method-LFQDataStats-new)

- [`LFQDataStats$stats()`](#method-LFQDataStats-stats)

- [`LFQDataStats$stats_wide()`](#method-LFQDataStats-stats_wide)

- [`LFQDataStats$stats_quantiles()`](#method-LFQDataStats-stats_quantiles)

- [`LFQDataStats$density()`](#method-LFQDataStats-density)

- [`LFQDataStats$density_median()`](#method-LFQDataStats-density_median)

- [`LFQDataStats$violin()`](#method-LFQDataStats-violin)

- [`LFQDataStats$violin_median()`](#method-LFQDataStats-violin_median)

- [`LFQDataStats$stdv_vs_mean()`](#method-LFQDataStats-stdv_vs_mean)

- [`LFQDataStats$power_t_test_quantiles()`](#method-LFQDataStats-power_t_test_quantiles)

- [`LFQDataStats$power_t_test()`](#method-LFQDataStats-power_t_test)

- [`LFQDataStats$clone()`](#method-LFQDataStats-clone)

------------------------------------------------------------------------

### Method `new()`

create analyse variances and CV

#### Usage

    LFQDataStats$new(lfqdata, stats = c("everything", "interaction", "all"))

#### Arguments

- `lfqdata`:

  LFQData object

- `stats`:

  if \`interaction\`, compute within-group stats; if \`all\`, compute
  overall CV; if \`pooled\`, use pooled variance with grouping
  information

------------------------------------------------------------------------

### Method `stats()`

access data.frame with statistics

#### Usage

    LFQDataStats$stats()

#### Returns

data.frame with computed statistics

------------------------------------------------------------------------

### Method `stats_wide()`

access data.frame with statistics in wide format

#### Usage

    LFQDataStats$stats_wide()

#### Returns

data.frame with computed statistics in wide format

------------------------------------------------------------------------

### Method `stats_quantiles()`

Determine CV or sd for the quantiles

#### Usage

    LFQDataStats$stats_quantiles(probs = c(0.1, 0.25, 0.5, 0.75, 0.9))

#### Arguments

- `probs`:

  for which quantile to determine CV or sd

------------------------------------------------------------------------

### Method [`density()`](https://rdrr.io/r/stats/density.html)

plots density or ecdf

#### Usage

    LFQDataStats$density(ggstat = c("density", "ecdf"))

#### Arguments

- `ggstat`:

  either density or ecdf

#### Returns

ggplot

------------------------------------------------------------------------

### Method `density_median()`

plot density or ecdf of CV or sd for the 50

#### Usage

    LFQDataStats$density_median(ggstat = c("density", "ecdf"))

#### Arguments

- `ggstat`:

  either density of ecdf

#### Returns

ggplot

------------------------------------------------------------------------

### Method `violin()`

plot violinplot of CV or sd

#### Usage

    LFQDataStats$violin()

#### Arguments

- `ggstat`:

  either density of ecdf

#### Returns

ggplot

------------------------------------------------------------------------

### Method `violin_median()`

plot violinplot of CV or sd for the 50

#### Usage

    LFQDataStats$violin_median()

#### Returns

ggplot

------------------------------------------------------------------------

### Method `stdv_vs_mean()`

plot sd vs mean

#### Usage

    LFQDataStats$stdv_vs_mean(size = 200)

#### Arguments

- `size`:

  number of points to sample (default 200)

#### Returns

ggplot

------------------------------------------------------------------------

### Method `power_t_test_quantiles()`

compute sample size for entire dataset

#### Usage

    LFQDataStats$power_t_test_quantiles(
      probs = c(0.1, 0.25, 0.5, 0.75, 0.9),
      delta = c(0.59, 1, 2),
      power = 0.8,
      sig.level = 0.05
    )

#### Arguments

- `probs`:

  quantiles of sd for which sample size should be computed

- `delta`:

  effect size

- `power`:

  power of test

- `sig.level`:

  significance level.

------------------------------------------------------------------------

### Method `power_t_test()`

compute sample for each protein

#### Usage

    LFQDataStats$power_t_test(delta = c(0.59, 1, 2), power = 0.8, sig.level = 0.05)

#### Arguments

- `delta`:

  effect size

- `power`:

  power of test

- `sig.level`:

  significance level.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    LFQDataStats$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# study variance of not normalized data
#source("c:/Users/wewol/prog/prolfqua/R/LFQData.R")
runallfuncs <- function(x){

  stopifnot("data.frame" %in% class(x$stats()))
  stopifnot("data.frame" %in% class(x$stats_wide()))
  stopifnot(c("long", "wide") %in% names(x$stats_quantiles()))
  stopifnot("ggplot" %in% class(x$density()))
  stopifnot("ggplot" %in% class(x$density_median()))
  stopifnot("ggplot" %in% class(x$density("ecdf")))
  stopifnot("ggplot" %in% class(x$density_median("ecdf")))
  stopifnot("ggplot" %in% class(x$violin()))
  stopifnot("ggplot" %in% class(x$violin_median()))
  stopifnot("ggplot" %in% class(x$stdv_vs_mean(size = 400)))
  if(!x$lfq$is_transformed()){
    stopifnot(is.null(x$power_t_test()))
    stopifnot(is.null(x$power_t_test_quantiles()))
  }
}
bb <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done



lfqdata <- LFQData$new(bb$data, bb$config)
lfqstats <- lfqdata$get_Stats()
#> completing cases
stopifnot(ncol(lfqstats$stats_wide()) == 30)
lfqstats$violin()
#> Warning: Removed 2 rows containing non-finite outside the scale range
#> (`stat_ydensity()`).
#> Warning: Removed 2 rows containing non-finite outside the scale range
#> (`stat_summary()`).

runallfuncs(lfqstats)
#> Warning: data is not transformed - aborting
#> Warning: data is not transformed - aborting
x <- lfqstats

#study variance of normalized data


lfqdata <- LFQData$new(bb$data, bb$config)
lfqdata$is_transformed(TRUE)
lfqstats <- lfqdata$get_Stats()
#> completing cases
stopifnot(ncol(lfqstats$stats_wide()) == 26)
runallfuncs(lfqstats)

#Slightly different dataset


# estimates statistics for all samples
lfqstats <- lfqdata$get_Stats(stats = "all")
#> completing cases
stopifnot(ncol(lfqstats$stats_wide()) == 8)
runallfuncs(lfqstats)
lfqstats <- lfqdata$get_Stats(stats = "interaction")
#> completing cases
stopifnot(ncol(lfqstats$stats_wide()) == 20)
runallfuncs(lfqstats)

# Group size 1
bb <- prolfqua::sim_lfq_data_peptide_config(N=1)
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
table_factors_size(bb$data,bb$config )
#> # A tibble: 3 × 2
#>   group_     n
#>   <chr>  <int>
#> 1 A          1
#> 2 B          1
#> 3 Ctrl       1
lfqdata <- LFQData$new(bb$data, bb$config)
lfqstats <- lfqdata$get_Stats()
#> completing cases

# stopifnot(ncol(lfqstats$stats_wide()) == 30)
lfqstats$violin()
#> Warning: Removed 2 rows containing non-finite outside the scale range
#> (`stat_ydensity()`).
#> Warning: Removed 2 rows containing non-finite outside the scale range
#> (`stat_summary()`).

runallfuncs(lfqstats)
#> Warning: data is not transformed - aborting
#> Warning: data is not transformed - aborting
```
