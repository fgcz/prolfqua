# Visualize intensity distributions

Visualize intensity distributions

visualize intensity distributions

## Usage

``` r
plot_intensity_distribution_violin(
  pdata,
  sample_name,
  response,
  is_transformed = FALSE
)

plot_intensity_distribution_density(
  pdata,
  sample_name,
  response,
  is_transformed = FALSE,
  legend = NA,
  max_legend_samples = 16
)
```

## Arguments

- pdata:

  data.frame

- sample_name:

  character - sample column name

- response:

  character - intensity column name

- is_transformed:

  logical - is intensity log-transformed?

- legend:

  show legend. If `NA`, hide automatically when the number of samples is
  larger than `max_legend_samples`.

- max_legend_samples:

  maximum number of samples for automatic legend display.

## Value

`plot_intensity_distribution_violin()`: a ggplot object with one violin
of the intensities per sample and the median marked.

`plot_intensity_distribution_density()`: a ggplot object with one
density line of the intensities per sample.

## Note

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Plotter()$intensity_distribution_violin()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Plotter()$intensity_distribution_density()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
