# Generate feature-level boxplots for all proteins

Generate feature-level boxplots for all proteins

## Usage

``` r
plot_hierarchies_boxplot_df(
  pdata,
  lfqdata,
  hierarchy = lfqdata$relevant_hierarchy_keys(),
  facet_grid_on = NULL
)
```

## Arguments

- pdata:

  data.frame

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

- hierarchy:

  e.g. protein_Id default relevant_hierarchy_keys

- facet_grid_on:

  default NULL

## Value

a tibble with one row per element of `hierarchy`, a column (named after
the `hierarchy` columns joined by `+`) identifying it, and a list column
`boxplot` of ggplot objects made by
[`plot_hierarchies_boxplot()`](https://fgcz.github.io/prolfqua/reference/plot_hierarchies_boxplot.md)

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$get_Plotter()$boxplots()`.
See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
