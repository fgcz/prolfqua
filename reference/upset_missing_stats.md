# Prepare data for an UpSetR plot of all samples

Prepare data for an UpSetR plot of all samples

## Usage

``` r
upset_missing_stats(lfqdata)
```

## Arguments

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

## Value

a list with `data`, a data.frame with one row per hierarchy element and
one 0/1 presence column per sample, and `nsets`, the number of sample
columns

## Note

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Plotter()$upset_missing()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
