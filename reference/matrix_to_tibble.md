# Convert a matrix or data.frame to a tibble (taken from tidyquant)

Convert a matrix or data.frame to a tibble (taken from tidyquant)

## Usage

``` r
matrix_to_tibble(x, preserve_row_names = "row.names", ...)
```

## Arguments

- x:

  a matrix

- preserve_row_names:

  give name to rownames column, if NULL discard rownames

- ...:

  further parameters passed to as_tibble

## Value

a tibble with the columns of `x`, preceded by a column
`preserve_row_names` holding the row names unless `preserve_row_names`
is `NULL`.

## Note

Deprecated in prolfqua 1.9.5; use `tibble::as_tibble(x, rownames = )`.
See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
