# Violin plot of sd, CV or mean for intensities below and above the median

Violin plot of sd, CV or mean for intensities below and above the median

## Usage

``` r
plot_stat_violin_median(
  pdata,
  factor_key,
  stat = c("CV", "meanAbundance", "sd")
)
```

## Arguments

- pdata:

  data.frame with statistics

- factor_key:

  character - factor column name for x-axis

- stat:

  either CV, mean or sd

## Value

a ggplot object with violins of `stat` per level of `factor_key` and the
median and interquartile range marked, faceted into features below and
above the median mean intensity

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$get_Stats()$violin_median()`.
See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
