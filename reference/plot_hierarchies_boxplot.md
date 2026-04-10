# plot peptides by factors and it's levels.

plot peptides by factors and it's levels.

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

  LFQData object

- facet_grid_on:

  on which variable to run facet_grid

- beeswarm:

  use beeswarm default FALSE

- show_mean:

  show mean values

- pb:

  progress bar

## Examples

``` r
istar <- sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq <- LFQData$new(istar$data, istar$config)
res <- plot_hierarchies_boxplot_df(lfq$get_data(), lfq)
res$boxplot[[1]]
#> Warning: Removed 7 rows containing non-finite outside the scale range
#> (`stat_boxplot()`).
#> Warning: Removed 7 rows containing non-finite outside the scale range
#> (`stat_summary()`).
#> Warning: Removed 7 rows containing non-finite outside the scale range
#> (`stat_summary()`).
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`position_quasirandom()`).


xnested <- lfq$get_data() |>
  dplyr::group_by(across(all_of(lfq$relevant_hierarchy_keys()))) |> tidyr::nest()
p <- plot_hierarchies_boxplot(xnested$data[[1]], xnested$protein_Id[[1]],
  lfq, beeswarm = FALSE, show_mean = TRUE)
p <- plot_hierarchies_boxplot(xnested$data[[1]], xnested$protein_Id[[1]],
  lfq, beeswarm = TRUE)
p
#> Warning: Removed 7 rows containing non-finite outside the scale range
#> (`stat_boxplot()`).
#> Warning: Removed 7 rows containing non-finite outside the scale range
#> (`stat_summary()`).
#> Warning: Removed 7 rows containing non-finite outside the scale range
#> (`stat_summary()`).
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`position_quasirandom()`).
```
