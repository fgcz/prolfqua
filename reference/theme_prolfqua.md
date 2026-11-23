# Prism-style ggplot theme used by all prolfqua plots

White panel, black axis lines and outward ticks, bold axis titles, bold
facet labels on a light grey strip, legend without box. Add it to your
own ggplots to match the look of prolfqua figures.

## Usage

``` r
theme_prolfqua(grid = FALSE, base_size = 12)
```

## Arguments

- grid:

  add faint major gridlines, for plots read against reference values
  (scatter, density)

- base_size:

  base font size

## Value

a ggplot2 theme

## See also

Other plotting:
[`ContrastsPlotter`](https://fgcz.github.io/prolfqua/reference/ContrastsPlotter.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`upset_interaction_missing_stats()`](https://fgcz.github.io/prolfqua/reference/upset_interaction_missing_stats.md)

## Examples

``` r
ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
  ggplot2::geom_point() +
  theme_prolfqua(grid = TRUE)
```
