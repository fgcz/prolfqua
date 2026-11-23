# Robustly scale a matrix by medians and MADs

Robustly scale a matrix by medians and MADs

## Usage

``` r
robust_scale(data, dim = 2, preserve_mean = FALSE)
```

## Arguments

- data:

  numeric matrix

- dim:

  margin to scale: 2 scales each column (sample), 1 each row

- preserve_mean:

  if TRUE, add the mean of the medians back

## Value

a matrix like `data`, centred by the median of each column (`dim = 2`)
or row (`dim = 1`) and divided by its MAD relative to the mean MAD; only
centred, with a warning, if a MAD is zero; with `preserve_mean = TRUE`
the mean of the medians is added back

## Note

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Transformer()$robscale()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
