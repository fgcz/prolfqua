# Plot a heatmap without clustering (use to show NAs)

Plot a heatmap without clustering (use to show NAs)

## Usage

``` r
plot_raster(
  matrix,
  annotation,
  factor_keys,
  sample_name,
  arrange = c("mean", "var"),
  not_na = FALSE,
  show_rownames = FALSE,
  max_rownames_chars = 60,
  max_sample_label_chars = 20,
  ...
)
```

## Arguments

- matrix:

  numeric matrix - wide-format intensity data

- annotation:

  data.frame - sample annotation

- factor_keys:

  character vector - factor column names for annotation

- sample_name:

  character - sample name column

- arrange:

  either mean or var

- not_na:

  if true than arrange by nr of NA's first and then by arrange

- show_rownames:

  logical, show row names in heatmap

- max_rownames_chars:

  maximum displayed row label length

- max_sample_label_chars:

  maximum displayed sample label length. Labels keep their suffix
  because sample prefixes are often shared.

- ...:

  additional arguments passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)

## Value

invisibly, a ComplexHeatmap `Heatmap` object of the abundances with rows
ordered by `arrange` and unclustered columns; `NULL` with a warning if
`matrix` has fewer than two rows

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$get_Plotter()$raster()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
