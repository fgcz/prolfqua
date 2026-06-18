# plot PCA

plot PCA

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

  numeric matrix — wide-format intensity data

- annotation:

  data.frame — sample annotation

- sample_name:

  character — sample name column

- factor_keys:

  character vector — factor column names (first for color, second for
  shape)

- PC:

  which principal components to plot

- add_txt:

  show sample labels

- nudge:

  label nudge distance

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
[`plot_hierarchies_boxplot_df()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_boxplot_df.md),
[`plot_hierarchies_line()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line.md),
[`plot_hierarchies_line_df()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line_df.md),
[`plot_intensity_distribution_violin()`](https://wolski.github.io/prolfqua/reference/plot_intensity_distribution_violin.md),
[`plot_na_heatmap()`](https://wolski.github.io/prolfqua/reference/plot_NA_heatmap.md),
[`plot_raster()`](https://wolski.github.io/prolfqua/reference/plot_raster.md),
[`plot_sample_correlation()`](https://wolski.github.io/prolfqua/reference/plot_sample_correlation.md),
[`upset_interaction_missing_stats()`](https://wolski.github.io/prolfqua/reference/UpSet_interaction_missing_stats.md),
[`upset_missing_stats()`](https://wolski.github.io/prolfqua/reference/UpSet_missing_stats.md)

## Examples

``` r
istar <- sim_lfq_data_protein_config(with_missing = TRUE, weight_missing = .8, Nprot = 3000)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq <- LFQData$new(istar$data, istar$config)
wide <- lfq$data_wide(as.matrix = TRUE)
tmp <- plot_pca(wide$data, wide$annotation, lfq$sample_name(), lfq$factor_keys(),
  add_txt = TRUE, nudge = 0.01)
#> PCA: removed 2459 of 2986 features with missing values. To keep all features, impute missing values first, e.g. AggregateLimpa$new(lfqdata, impute_only = TRUE)$aggregate().
stopifnot("ggplot" %in% class(tmp))
tmp <- plot_pca(wide$data, wide$annotation, lfq$sample_name(), lfq$factor_keys())
#> PCA: removed 2459 of 2986 features with missing values. To keep all features, impute missing values first, e.g. AggregateLimpa$new(lfqdata, impute_only = TRUE)$aggregate().
stopifnot("ggplot" %in% class(tmp))
```
