# Visualize correlation among samples

Visualize correlation among samples

## Usage

``` r
plot_sample_correlation(matrix, method = "pearson", ...)
```

## Arguments

- matrix:

  numeric matrix - wide-format intensity data (samples as columns)

- method:

  correlation method, as for
  [`plot_heatmap_cor()`](https://fgcz.github.io/prolfqua/reference/plot_heatmap_cor.md)

- ...:

  passed to
  [`corrplot::corrplot.mixed()`](https://rdrr.io/pkg/corrplot/man/corrplot.mixed.html),
  e.g. `tl.cex` and `number.cex` for text size (default 1)

## Value

invisibly, the list returned by
[`corrplot::corrplot.mixed()`](https://rdrr.io/pkg/corrplot/man/corrplot.mixed.html)
with the correlation matrix `corr` and the cell positions `corrPos`;
called for its side effect of drawing the plot

## Note

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Plotter()$sample_correlation()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
