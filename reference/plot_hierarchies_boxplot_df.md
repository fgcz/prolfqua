# generates peptide level plots for all Proteins

generates peptide level plots for all Proteins

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

  LFQData object

- hierarchy:

  e.g. protein_Id default relevant_hierarchy_keys

- facet_grid_on:

  default NULL

## See also

Other plotting:
[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`medpolish_estimate_df()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate_df.md),
[`missigness_histogram()`](https://wolski.github.io/prolfqua/reference/missigness_histogram.md),
[`missingness_per_condition()`](https://wolski.github.io/prolfqua/reference/missingness_per_condition.md),
[`missingness_per_condition_cumsum()`](https://wolski.github.io/prolfqua/reference/missingness_per_condition_cumsum.md),
[`plot_estimate()`](https://wolski.github.io/prolfqua/reference/plot_estimate.md),
[`plot_heatmap()`](https://wolski.github.io/prolfqua/reference/plot_heatmap.md),
[`plot_heatmap_cor()`](https://wolski.github.io/prolfqua/reference/plot_heatmap_cor.md),
[`plot_hierarchies_add_quantline()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_add_quantline.md),
[`plot_hierarchies_line()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line.md),
[`plot_hierarchies_line_df()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line_df.md),
[`plot_intensity_distribution_violin()`](https://wolski.github.io/prolfqua/reference/plot_intensity_distribution_violin.md),
[`plot_na_heatmap()`](https://wolski.github.io/prolfqua/reference/plot_NA_heatmap.md),
[`plot_pca()`](https://wolski.github.io/prolfqua/reference/plot_pca.md),
[`plot_raster()`](https://wolski.github.io/prolfqua/reference/plot_raster.md),
[`plot_sample_correlation()`](https://wolski.github.io/prolfqua/reference/plot_sample_correlation.md),
[`upset_interaction_missing_stats()`](https://wolski.github.io/prolfqua/reference/UpSet_interaction_missing_stats.md),
[`upset_missing_stats()`](https://wolski.github.io/prolfqua/reference/UpSet_missing_stats.md)

## Examples

``` r
istar <- sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq <- LFQData$new(istar$data, istar$config)
res <- plot_hierarchies_boxplot_df(lfq$data_long(), lfq)
res$boxplot[[1]]
#> Warning: Removed 7 rows containing non-finite outside the scale range
#> (`stat_boxplot()`).
#> Warning: Removed 7 rows containing non-finite outside the scale range
#> (`stat_summary()`).
#> Warning: Removed 7 rows containing non-finite outside the scale range
#> (`stat_summary()`).
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`position_quasirandom()`).


lfq2 <- LFQData$new(
  istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 2)),
  istar$config)
res <- plot_hierarchies_boxplot_df(lfq2$data_long(), lfq2)
res$boxplot[[1]]
#> Warning: Removed 1 row containing non-finite outside the scale range (`stat_boxplot()`).
#> Warning: Removed 1 row containing non-finite outside the scale range (`stat_summary()`).
#> Warning: Removed 1 row containing non-finite outside the scale range (`stat_summary()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`position_quasirandom()`).
```
