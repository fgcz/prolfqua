# Sample subset of proteins/peptides/precursors

Sample subset of proteins/peptides/precursors

## Usage

``` r
sample_subset(size, pdata, hierarchy_keys_depth)
```

## Arguments

- size:

  size of sample

- pdata:

  tidy table

- hierarchy_keys_depth:

  character vector - hierarchy columns at current depth

## Value

a tibble with all rows of `pdata` for `size` randomly drawn hierarchy
elements

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$get_sample()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
