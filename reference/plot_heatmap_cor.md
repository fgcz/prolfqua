# plot correlation heatmap with annotations

plot correlation heatmap with annotations

## Usage

``` r
plot_heatmap_cor(
  matrix,
  annotation,
  factor_keys,
  sample_name,
  R2 = FALSE,
  color = colorRampPalette(c("white", "red"))(1024),
  max_sample_label_chars = 20,
  ...
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

- R2:

  logical — plot R-squared instead of correlation

- color:

  color palette

- max_sample_label_chars:

  maximum displayed sample label length. Labels keep their suffix
  because sample prefixes are often shared.

- ...:

  passed to \[ComplexHeatmap::Heatmap()\]

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
[`plot_hierarchies_add_quantline()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_add_quantline.md),
[`plot_hierarchies_boxplot_df()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_boxplot_df.md),
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
istar <- sim_lfq_data_protein_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq <- LFQData$new(istar$data, istar$config)
wide <- lfq$data_wide(as.matrix = TRUE)
pheat_map <- plot_heatmap_cor(wide$data, wide$annotation,
  lfq$factor_keys(), lfq$sample_name())
stopifnot(methods::is(pheat_map, "Heatmap"))
pheat_map <- plot_heatmap_cor(wide$data, wide$annotation,
  lfq$factor_keys(), lfq$sample_name(), R2 = TRUE)
stopifnot(methods::is(pheat_map, "Heatmap"))
```
