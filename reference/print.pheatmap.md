# Print method for pheatmap objects

Fixes pheatmap's print method which omits `grid.newpage()`, causing the
heatmap to draw on top of the previous plot when multiple plots are
produced in the same knitr chunk or graphics device.

## Usage

``` r
# S3 method for class 'pheatmap'
print(x, ...)
```

## Arguments

- x:

  a pheatmap object

- ...:

  ignored

## Value

The computed result.
