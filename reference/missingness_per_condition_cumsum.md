# Cumulative counts of missing values per condition

Cumulative counts of missing values per condition

## Usage

``` r
missingness_per_condition_cumsum(
  lfqdata,
  factors = lfqdata$relevant_factor_keys()
)
```

## Arguments

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

- factors:

  factor columns to use

## Value

a list with `data`, a wide tibble with one row per condition and one
column per number of missing values holding the cumulative feature
count, and `figure`, a ggplot bar chart of these counts

## Note

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Summariser()$missingness_per_group_cumsum()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
