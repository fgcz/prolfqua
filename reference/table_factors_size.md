# Table of distinct factors with group sizes

Table of distinct factors with group sizes

## Usage

``` r
table_factors_size(
  pdata,
  file_name,
  sample_name,
  factor_keys,
  factor_keys_depth
)
```

## Arguments

- pdata:

  data.frame

- file_name:

  character - file name column

- sample_name:

  character - sample name column

- factor_keys:

  character vector - all factor column names

- factor_keys_depth:

  character vector - factor columns at current depth

## Value

a tibble with one row per combination of the `factor_keys_depth` levels
and a column `n` with the number of samples

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$factors()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
