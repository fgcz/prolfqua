# Plot the density or ECDF of sd, mean or CV

Plot the density or ECDF of sd, mean or CV

## Usage

``` r
plot_stat_density(
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

  character - factor column name for colouring

- stat:

  sd, mean or CV

- ggstat:

  either density or ecdf

## Value

a ggplot object with one density or ECDF line of `stat` per level of
`factor_key`

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$get_Stats()$density()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
