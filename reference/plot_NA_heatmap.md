# plot heatmap of NA values

plot heatmap of NA values

## Usage

``` r
plot_na_heatmap(
  matrix,
  annotation,
  factor_keys,
  sample_name,
  limitrows = 10000,
  distance = "binary"
)
```

## Arguments

- matrix:

  numeric matrix — wide-format intensity data

- annotation:

  data.frame — sample annotation

- factor_keys:

  character vector — factor column names for annotation

- sample_name:

  character — sample name column

- limitrows:

  max rows to display

- distance:

  distance method for clustering

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
wide <- lfq$data_wide(as.matrix = TRUE)
tmp <- plot_na_heatmap(wide$data, wide$annotation, lfq$factor_keys(), lfq$sample_name())
#> rows with NA's: 16; all rows :28
stopifnot(methods::is(tmp, "Heatmap"))
```
