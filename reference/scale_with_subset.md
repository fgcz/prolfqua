# Scale data using a subset of the data

this should reduce the overall variance.

## Usage

``` r
scale_with_subset(lfqdata, lfqsubset, preserve_mean = FALSE, colname = NULL)
```

## Arguments

- lfqdata:

  LFQData object with full dataset

- lfqsubset:

  LFQData object with subset for computing scales

- preserve_mean:

  default FALSE - sets mean to zero

- colname:

  name of the scaled intensity column; default NULL uses
  `<response>_subset_scaled`

## Value

list with data, scales, and colname

## Note

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Transformer()$robscale_subset()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
