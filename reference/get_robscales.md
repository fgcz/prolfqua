# Compute the median and MAD of each sample

Compute the median and MAD of each sample

## Usage

``` r
get_robscales(lfqdata)
```

## Arguments

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

## Value

a list with the numeric vectors `medians` and `mads`, with one element
per sample

## Note

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Transformer()$get_scales()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
