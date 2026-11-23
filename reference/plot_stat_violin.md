# Violin plot of sd, CV or mean

Violin plot of sd, CV or mean

## Usage

``` r
plot_stat_violin(
  pdata,
  factor_keys_depth,
  stat = c("CV", "meanAbundance", "sd")
)
```

## Arguments

- pdata:

  data.frame with statistics

- factor_keys_depth:

  character vector - factor columns for grouping

- stat:

  either CV, mean or sd

## Value

a ggplot object with one violin of `stat` per combination of factor
levels, with the median marked

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$get_Stats()$violin()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
