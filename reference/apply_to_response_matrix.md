# Apply function requiring a matrix to tidy table

Apply function requiring a matrix to tidy table

## Usage

``` r
apply_to_response_matrix(lfqdata, .func, .funcname = NULL)
```

## Arguments

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

- .func:

  function taking and returning a matrix

- .funcname:

  name of function (used for creating new column)

## Value

list with `data` (data.frame) and `colname` (new column name)

## Note

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Transformer()$intensity_matrix()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
