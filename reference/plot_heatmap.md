# Plot a row z-score heatmap with annotations

Plot a row z-score heatmap with annotations

## Usage

``` r
plot_heatmap(
  matrix,
  annotation,
  factor_keys,
  sample_name,
  na_fraction = 0.4,
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

- na_fraction:

  maximum fraction of missing values per row. Rows that cannot be
  meaningfully z-scored are always removed.

- show_rownames:

  if TRUE shows row names, default FALSE

- max_rownames_chars:

  maximum displayed row label length

- max_sample_label_chars:

  maximum displayed sample label length. Labels keep their suffix
  because sample prefixes are often shared.

- ...:

  passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)

## Value

invisibly, a ComplexHeatmap `Heatmap` object of the row z-scored data;
`NULL` with a warning if `matrix` has no rows or no row can be z-scored

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$get_Plotter()$heatmap()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
