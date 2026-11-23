# Count distinct elements for each level of hierarchy

E.g. number of proteins, peptides, precursors in the dataset. Counts
over the whole dataset, whether observed or not.

## Usage

``` r
hierarchy_counts(pdata, hierarchy_keys)
```

## Arguments

- pdata:

  data.frame

- hierarchy_keys:

  character vector - all hierarchy column names

## Value

a one-row tibble with one column per hierarchy key holding the number of
distinct elements

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$hierarchy_counts()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
