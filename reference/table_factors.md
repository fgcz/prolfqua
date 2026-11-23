# Table of distinct factors (sample annotation)

Table of distinct factors (sample annotation)

## Usage

``` r
table_factors(pdata, file_name, sample_name, factor_keys)
```

## Arguments

- pdata:

  data.frame

- file_name:

  character - file name column

- sample_name:

  character - sample name column

- factor_keys:

  character vector - factor column names

## Value

a tibble with the distinct combinations of the file name, sample name
and factor columns, sorted by sample name

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$factors()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
