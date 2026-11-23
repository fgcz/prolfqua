# Complete cases

The tidy table does not need to contain missing data. This function
re-establishes the missing observations in a sample.

## Usage

``` r
complete_cases(lfqdata)
```

## Arguments

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

## Value

a tibble of the long-format data with one row for every combination of
sample and hierarchy element; added rows have `NA` intensities

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$complete_cases()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
