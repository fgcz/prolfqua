# Plot standard deviation against mean to assess variance stability

Plot standard deviation against mean to assess variance stability

## Usage

``` r
plot_stdv_vs_mean(pdata, factor_keys_depth, size = 2000)
```

## Arguments

- pdata:

  data.frame with statistics

- factor_keys_depth:

  character vector - factor columns for faceting

- size:

  how many points to sample (since scatter plot to slow for all)

## Value

a ggplot scatter plot of `sd` against `meanAbundance` with a loess fit,
one facet per combination of factor levels

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$get_Stats()$stdv_vs_mean()`.
See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
