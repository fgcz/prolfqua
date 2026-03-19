# plot peptides by factors and it's levels.

plot peptides by factors and it's levels.

## Usage

``` r
plot_hierarchies_boxplot(
  pdata,
  title,
  config,
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

- config:

  AnalysisConfiguration

- facet_grid_on:

  on which variable to run facet_grid

- beeswarm:

  use beeswarm default FALSE

## Examples

``` r


bb <- prolfqua_data('data_skylineSRM_HL_A')
config <- bb$config_f()
analysis <- bb$analysis(bb$data, config)
#> creating sampleName from fileName column
#> Warning: no nr_children column specified in the data, adding column nr_children and setting to 1.
#> completing cases
#> completing cases done
#> setup done
data <- prolfqua::transform_work_intensity(analysis, config, log2)
#> Column added : log2_Area
res <- plot_hierarchies_boxplot_df(data, config)
res$boxplot[[1]]
#> Warning: Removed 143 rows containing non-finite outside the scale
#> range (`stat_boxplot()`).
#> Warning: Removed 143 rows containing non-finite outside the scale
#> range (`stat_summary()`).
#> Warning: Removed 143 rows containing non-finite outside the scale
#> range (`stat_summary()`).
#> Warning: Removed 143 rows containing missing values or values outside
#> the scale range (`position_quasirandom()`).


hierarchy = config$hierarchy_keys_depth()
xnested <- data |> dplyr::group_by(across(all_of(hierarchy))) |> tidyr::nest()
p <- plot_hierarchies_boxplot(xnested$data[[1]], xnested$protein_Id[[1]],
  config, beeswarm = FALSE, show_mean = TRUE)
p <- plot_hierarchies_boxplot(xnested$data[[1]], xnested$protein_Id[[1]],
  config, beeswarm = TRUE)
p <- plot_hierarchies_boxplot(xnested$data[[1]], xnested$protein_Id[[1]],
  config, beeswarm = TRUE, facet_grid_on = "precursor_Id")
p
#> Warning: Removed 143 rows containing non-finite outside the scale
#> range (`stat_boxplot()`).
#> Warning: Removed 143 rows containing non-finite outside the scale
#> range (`stat_summary()`).
#> Warning: Removed 143 rows containing non-finite outside the scale
#> range (`stat_summary()`).
#> Warning: Removed 143 rows containing missing values or values outside
#> the scale range (`position_quasirandom()`).
```
