# Takes matrix of responses and converts into tibble

Takes matrix of responses and converts into tibble

## Usage

``` r
response_matrix_as_tibble(pdata, value, config, data = NULL, sep = "~lfq~")
```

## Arguments

- pdata:

  matrix with rownames encoding hierarchy keys

- value:

  name of column to store values in

- config:

  AnalysisConfiguration (needed for column name mapping)

- data:

  optional data.frame to join back to

- sep:

  separator used to unite the hierarchy keys

## Value

a long-format tibble with the hierarchy key, sample name and `value`
columns; if `data` is given, `data` joined with these values

## Note

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Transformer()$intensity_matrix()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
