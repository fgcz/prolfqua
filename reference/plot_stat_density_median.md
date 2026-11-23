# Plot the density or ECDF of sd, mean or CV for intensities below and above the median

Plot the density or ECDF of sd, mean or CV for intensities below and
above the median

## Usage

``` r
plot_stat_density_median(
  pdata,
  factor_key,
  stat = c("CV", "meanAbundance", "sd"),
  ggstat = c("density", "ecdf")
)
```

## Arguments

- pdata:

  data.frame with statistics

- factor_key:

  character - factor column name for faceting

- stat:

  sd, mean or CV

- ggstat:

  either density or ecdf

## Value

a ggplot object with density or ECDF lines of `stat` for features above
and below the median mean intensity, faceted by `factor_key`

## Note

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Stats()$density_median()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
