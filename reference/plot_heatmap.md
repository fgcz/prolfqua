# plot heatmap with annotations

plot heatmap with annotations

## Usage

``` r
plot_heatmap(data, config, na_fraction = 0.4, show_rownames = FALSE, ...)
```

## Arguments

- na_fraction:

  fraction of NA values per row

- show_rownames:

  if TRUE shows row names, default FALSE

## See also

Other plotting:
[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`UpSet_interaction_missing_stats()`](https://wolski.github.io/prolfqua/reference/UpSet_interaction_missing_stats.md),
[`UpSet_missing_stats()`](https://wolski.github.io/prolfqua/reference/UpSet_missing_stats.md),
[`medpolish_estimate_df()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate_df.md),
[`missigness_histogram()`](https://wolski.github.io/prolfqua/reference/missigness_histogram.md),
[`missingness_per_condition()`](https://wolski.github.io/prolfqua/reference/missingness_per_condition.md),
[`missingness_per_condition_cumsum()`](https://wolski.github.io/prolfqua/reference/missingness_per_condition_cumsum.md),
[`plot_NA_heatmap()`](https://wolski.github.io/prolfqua/reference/plot_NA_heatmap.md),
[`plot_estimate()`](https://wolski.github.io/prolfqua/reference/plot_estimate.md),
[`plot_heatmap_cor()`](https://wolski.github.io/prolfqua/reference/plot_heatmap_cor.md),
[`plot_heatmap_cor_iheatmap()`](https://wolski.github.io/prolfqua/reference/plot_heatmap_cor_iheatmap.md),
[`plot_hierarchies_add_quantline()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_add_quantline.md),
[`plot_hierarchies_boxplot_df()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_boxplot_df.md),
[`plot_hierarchies_line()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line.md),
[`plot_hierarchies_line_df()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line_df.md),
[`plot_intensity_distribution_violin()`](https://wolski.github.io/prolfqua/reference/plot_intensity_distribution_violin.md),
[`plot_pca()`](https://wolski.github.io/prolfqua/reference/plot_pca.md),
[`plot_raster()`](https://wolski.github.io/prolfqua/reference/plot_raster.md),
[`plot_sample_correlation()`](https://wolski.github.io/prolfqua/reference/plot_sample_correlation.md),
[`plot_screeplot()`](https://wolski.github.io/prolfqua/reference/plot_screeplot.md)

## Examples

``` r

istar <- sim_lfq_data_protein_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
config <- istar$config
analysis <- istar$data


p  <- plot_heatmap(analysis, config)
stopifnot(class(p) == "pheatmap")
p2 <- plot_heatmap(analysis, config, show_rownames = TRUE)
stopifnot(class(p) == "pheatmap")
```
