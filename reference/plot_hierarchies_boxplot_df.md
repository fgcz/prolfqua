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

## See also

Other plotting:
[`ContrastsPlotter`](https://fgcz.github.io/prolfqua/reference/ContrastsPlotter.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`missingness_histogram()`](https://fgcz.github.io/prolfqua/reference/missingness_histogram.md),
[`missingness_per_condition()`](https://fgcz.github.io/prolfqua/reference/missingness_per_condition.md),
[`missingness_per_condition_cumsum()`](https://fgcz.github.io/prolfqua/reference/missingness_per_condition_cumsum.md),
[`plot_estimate()`](https://fgcz.github.io/prolfqua/reference/plot_estimate.md),
[`plot_heatmap()`](https://fgcz.github.io/prolfqua/reference/plot_heatmap.md),
[`plot_heatmap_cor()`](https://fgcz.github.io/prolfqua/reference/plot_heatmap_cor.md),
[`plot_hierarchies_line()`](https://fgcz.github.io/prolfqua/reference/plot_hierarchies_line.md),
[`plot_intensity_distribution_violin()`](https://fgcz.github.io/prolfqua/reference/plot_intensity_distribution_violin.md),
[`plot_na_heatmap()`](https://fgcz.github.io/prolfqua/reference/plot_na_heatmap.md),
[`plot_pca()`](https://fgcz.github.io/prolfqua/reference/plot_pca.md),
[`plot_raster()`](https://fgcz.github.io/prolfqua/reference/plot_raster.md),
[`plot_sample_correlation()`](https://fgcz.github.io/prolfqua/reference/plot_sample_correlation.md),
[`theme_prolfqua()`](https://fgcz.github.io/prolfqua/reference/theme_prolfqua.md),
[`upset_interaction_missing_stats()`](https://fgcz.github.io/prolfqua/reference/upset_interaction_missing_stats.md),
[`upset_missing_stats()`](https://fgcz.github.io/prolfqua/reference/upset_missing_stats.md)

## Examples

``` r
istar <- sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
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
```
