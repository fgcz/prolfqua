# plot PCA

plot PCA

## Usage

``` r
plot_pca(
  data,
  config,
  PC = c(1, 2),
  add_txt = FALSE,
  plotly = FALSE,
  nudge = 0.1
)
```

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
[`plot_heatmap()`](https://wolski.github.io/prolfqua/reference/plot_heatmap.md),
[`plot_heatmap_cor()`](https://wolski.github.io/prolfqua/reference/plot_heatmap_cor.md),
[`plot_hierarchies_add_quantline()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_add_quantline.md),
[`plot_hierarchies_boxplot_df()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_boxplot_df.md),
[`plot_hierarchies_line()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line.md),
[`plot_hierarchies_line_df()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line_df.md),
[`plot_intensity_distribution_violin()`](https://wolski.github.io/prolfqua/reference/plot_intensity_distribution_violin.md),
[`plot_raster()`](https://wolski.github.io/prolfqua/reference/plot_raster.md),
[`plot_sample_correlation()`](https://wolski.github.io/prolfqua/reference/plot_sample_correlation.md)

## Examples

``` r

istar <- sim_lfq_data_protein_config(with_missing = TRUE, weight_missing = .8, Nprot = 3000)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
config <- istar$config
analysis <- istar$data
tmp <- plot_pca(analysis, config, add_txt= TRUE, nudge = 0.01)
#> PCA: removed 2459 of 2986 features with missing values. For PCA with all features, impute first using impute_with_zcomp().
#> Joining with `by = join_by(sampleName)`
print(tmp)


stopifnot("ggplot" %in% class(tmp) )
tmp <- plot_pca(analysis, config, add_txt= FALSE)
#> PCA: removed 2459 of 2986 features with missing values. For PCA with all features, impute first using impute_with_zcomp().
#> Joining with `by = join_by(sampleName)`
stopifnot("ggplot" %in% class(tmp) )
tmp <- plot_pca(analysis, config, PC = c(1,2))
#> PCA: removed 2459 of 2986 features with missing values. For PCA with all features, impute first using impute_with_zcomp().
#> Joining with `by = join_by(sampleName)`
stopifnot("ggplot" %in% class(tmp) )
tmp <- plot_pca(analysis, config, PC = c(2,40))
#> PCA: removed 2459 of 2986 features with missing values. For PCA with all features, impute first using impute_with_zcomp().
#> Warning: nr of PCs: 13
print(tmp)
#> NULL
```
