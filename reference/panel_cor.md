# Correlation panel for pairs plots (default lower panel of `pairs_smooth()`)

Correlation panel for pairs plots (default lower panel of
[`pairs_smooth()`](https://fgcz.github.io/prolfqua/reference/pairs_smooth.md))

## Usage

``` r
panel_cor(x, y, digits = 2, ...)
```

## Arguments

- x:

  numeric data

- y:

  numeric data

- digits:

  number of digits to display

- ...:

  not used

## Value

`NULL`, invisibly; called for its side effect of writing the
correlation, R^2 and p-value into the current panel

## Note

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Plotter()$pairs_smooth()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
