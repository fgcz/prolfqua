# Separates hierarchy columns into starting columns

Separates hierarchy columns into starting columns

## Usage

``` r
separate_hierarchy(data, config)
```

## Arguments

- data:

  data.frame

- config:

  AnalysisConfiguration

## Value

`data` with each hierarchy key at the current depth split at
`config$sep` into its component columns, which are added after the key
column

## Note

Deprecated in prolfqua 1.9.5; use `prolfquapp::separate_hierarchy()`.
See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
