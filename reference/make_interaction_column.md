# Create an interaction column from factors

Create an interaction column from factors

## Usage

``` r
make_interaction_column(data, columns, sep = ".")
```

## Arguments

- data:

  data frame

- columns:

  factor columns to combine

- sep:

  separator between the levels in the combined labels

## Value

`data` with an added factor column `interaction` combining the values of
`columns`, each prefixed with its column name and joined by `sep`

## Note

Deprecated in prolfqua 1.9.5; use
[`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md).
See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
