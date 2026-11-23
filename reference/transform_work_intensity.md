# Transform intensity column by applying a function

Transform intensity column by applying a function

## Usage

``` r
transform_work_intensity(pdata, response, .func, .funcname = NULL)
```

## Arguments

- pdata:

  data.frame

- response:

  character, name of the response column to transform

- .func:

  function to transform intensities e.g. log2

- .funcname:

  name of function (used for creating new column name)

## Value

list with `data` (data.frame), `colname` (new column name) and `zeros`
(logical, rows whose zero intensity became `NA`)

## Details

For `log`, `log2` and `log10`, zero intensities mean "not detected" and
become `NA`; LFQData records them in its provenance. Negative
intensities give `NaN` with a warning, which fitting and aggregation
then reject.

## Note

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Transformer()$intensity_array()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
