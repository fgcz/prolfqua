# Plot a sample correlation heatmap with annotations

Plot a sample correlation heatmap with annotations

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
  method = "pearson",
  ...
)
```

## Arguments

- matrix:

  numeric matrix - wide-format intensity data

- annotation:

  data.frame - sample annotation

- factor_keys:

  character vector - factor column names for annotation

- sample_name:

  character - sample name column

- R2:

  logical - plot R-squared instead of correlation

- color:

  color palette

- max_sample_label_chars:

  maximum displayed sample label length. Labels keep their suffix
  because sample prefixes are often shared.

- method:

  correlation method: `"pearson"` (default), `"spearman"` or `"kendall"`
  for [`stats::cor()`](https://rdrr.io/r/stats/cor.html) on pairwise
  complete observations, or a function of the matrix returning a
  sample-by-sample correlation matrix

- ...:

  passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)

## Value

invisibly, a ComplexHeatmap `Heatmap` object of the clustered sample
correlations (or R^2 values)

## See also

Other plotting:
[`ContrastsPlotter`](https://fgcz.github.io/prolfqua/reference/ContrastsPlotter.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`missingness_histogram()`](https://fgcz.github.io/prolfqua/reference/missingness_histogram.md),
[`missingness_per_condition()`](https://fgcz.github.io/prolfqua/reference/missingness_per_condition.md),
[`missingness_per_condition_cumsum()`](https://fgcz.github.io/prolfqua/reference/missingness_per_condition_cumsum.md),
[`plot_estimate()`](https://fgcz.github.io/prolfqua/reference/plot_estimate.md),
[`plot_heatmap()`](https://fgcz.github.io/prolfqua/reference/plot_heatmap.md),
[`plot_hierarchies_boxplot_df()`](https://fgcz.github.io/prolfqua/reference/plot_hierarchies_boxplot_df.md),
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
istar <- sim_lfq_data_protein_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 5 unconfigured columns: average_prot_abundance, sd, mean, N, Replicate.
#> setup_analysis: added 16 rows with NA response for feature and sample combinations absent from the input.
lfq <- LFQData$new(istar$data, istar$config)
wide <- lfq$data_wide(as.matrix = TRUE)
pheat_map <- plot_heatmap_cor(wide$data, wide$annotation,
  lfq$factor_keys(), lfq$sample_name())
stopifnot(methods::is(pheat_map, "Heatmap"))
pheat_map <- plot_heatmap_cor(wide$data, wide$annotation,
  lfq$factor_keys(), lfq$sample_name(), R2 = TRUE)
stopifnot(methods::is(pheat_map, "Heatmap"))
# rank correlation, or any function returning a correlation matrix
spearman <- plot_heatmap_cor(wide$data, wide$annotation,
  lfq$factor_keys(), lfq$sample_name(), method = "spearman")
complete <- plot_heatmap_cor(wide$data, wide$annotation,
  lfq$factor_keys(), lfq$sample_name(), method = function(x) stats::cor(x, use = "complete.obs"))
```
