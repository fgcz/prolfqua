# Boxplot of feature intensities by factor levels

Boxplot of feature intensities by factor levels

## Usage

``` r
plot_hierarchies_boxplot(
  pdata,
  title,
  lfqdata,
  facet_grid_on = NULL,
  beeswarm = TRUE,
  show_mean = TRUE,
  pb
)
```

## Arguments

- pdata:

  data.frame

- title:

  name to show

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

- facet_grid_on:

  on which variable to run facet_grid

- beeswarm:

  use beeswarm default FALSE

- show_mean:

  show mean values

- pb:

  progress bar

## Value

a ggplot object with boxplots of the intensities in `pdata` per
combination of factor levels, optionally with the individual points and
the means

## Note

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Plotter()$write_boxplots()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
