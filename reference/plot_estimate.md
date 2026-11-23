# Plot feature data and result of aggregation

Plot feature data and result of aggregation

## Usage

``` r
plot_estimate(lfqdata, lfqdata_agg, show.legend = FALSE)
```

## Arguments

- lfqdata:

  LFQData object (original data)

- lfqdata_agg:

  LFQData object (aggregated data)

- show.legend:

  logical, show legend in plot

## Value

a tibble with one row per aggregated element and the columns
`hierarchy_id`, `data` (the feature intensities), `other` (the
aggregated intensities) and `plots` (ggplot objects of the feature
intensities with the aggregate drawn in black)

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$get_Aggregator()$plot()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
