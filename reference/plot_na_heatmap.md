# Plot a heatmap of missing values

Plot a heatmap of missing values

## Usage

``` r
plot_na_heatmap(
  matrix,
  annotation,
  factor_keys,
  sample_name,
  limitrows = 10000,
  distance = "binary"
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

- limitrows:

  max rows to display; above it, `limitrows` evenly spaced rows

- distance:

  distance method for clustering

## Value

a ComplexHeatmap `Heatmap` object marking the missing values of the rows
that have any; `NULL` if fewer than two such rows

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$get_Plotter()$na_heatmap()`.
See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
