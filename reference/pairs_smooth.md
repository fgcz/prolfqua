# Pairs plot with smoothed scatter panels

Pairs plot with smoothed scatter panels

## Usage

``` r
pairs_smooth(data, legend = FALSE, ...)
```

## Arguments

- data:

  data matrix or data.frame as normally passed to pairs

- legend:

  add legend to plots

- ...:

  params usually passed to pairs

## Value

`NULL`, invisibly; called for its side effect of drawing the pairs plot

## Note

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Plotter()$pairs_smooth()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
