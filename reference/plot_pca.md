# Plot a principal component analysis of the samples

Plot a principal component analysis of the samples

## Usage

``` r
plot_pca(
  matrix,
  annotation,
  sample_name,
  factor_keys,
  PC = c(1, 2),
  add_txt = FALSE,
  nudge = 0.1
)
```

## Arguments

- matrix:

  numeric matrix - wide-format intensity data

- annotation:

  data.frame - sample annotation

- sample_name:

  character - sample name column

- factor_keys:

  character vector - factor column names (first for color, second for
  shape)

- PC:

  which principal components to plot

- add_txt:

  show sample labels

- nudge:

  label nudge distance

## Value

a ggplot object of the sample scores on the principal components `PC`,
coloured by the first and shaped by the second of `factor_keys`

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$get_Plotter()$pca()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
