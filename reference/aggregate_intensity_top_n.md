# Aggregates top N intensities

run
[rank_peptide_by_intensity](https://fgcz.github.io/prolfqua/reference/rank_peptide_by_intensity.md)
first

## Usage

``` r
aggregate_intensity_top_n(ranked_data, lfqdata, .func, N = 3)
```

## Arguments

- ranked_data:

  data.frame with ranked peptides

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

- .func:

  function to use for aggregation

- N:

  default 3 top intensities.

## Value

list with data and new reduced configuration (config)

## Note

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Aggregator("topN")$aggregate()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
