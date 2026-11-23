# Plot a sample correlation heatmap with annotations

Plot a sample correlation heatmap with annotations

## Usage

``` r
plot_heatmap_cor(
  matrix,
  annotation,
  factor_keys,
  sample_name,
  R2 = FALSE,
  color = c("white", "red"),
  max_sample_label_chars = 20,
  method = "pearson",
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

- R2:

  logical - plot R-squared instead of correlation

- color:

  colours, spread evenly from the lowest to the highest correlation
  (default white to red)

- max_sample_label_chars:

  maximum displayed sample label length. Labels keep their suffix
  because sample prefixes are often shared.

- method:

  correlation method: `"pearson"` (default), `"spearman"` or `"kendall"`
  for [`stats::cor()`](https://rdrr.io/r/stats/cor.html) on pairwise
  complete observations, or a function of the matrix returning a
  sample-by-sample correlation matrix

- ...:

  passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)

## Value

invisibly, a ComplexHeatmap `Heatmap` object of the clustered sample
correlations (or R^2 values)

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$get_Plotter()$heatmap_cor()`.
See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
