# visualize correlation among samples

visualize correlation among samples

## Usage

``` r
plot_sample_correlation(matrix, method = "pearson", ...)
```

## Arguments

- matrix:

  numeric matrix — wide-format intensity data (samples as columns)

- method:

  correlation method, as for \[plot_heatmap_cor()\]

- ...:

  passed to \[corrplot::corrplot.mixed()\], e.g. \`tl.cex\` and
  \`number.cex\` for text size (default 1)

## See also

Other plotting:
[`ContrastsPlotter`](https://fgcz.github.io/prolfqua/reference/ContrastsPlotter.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`missigness_histogram()`](https://fgcz.github.io/prolfqua/reference/missigness_histogram.md),
[`missingness_per_condition()`](https://fgcz.github.io/prolfqua/reference/missingness_per_condition.md),
[`missingness_per_condition_cumsum()`](https://fgcz.github.io/prolfqua/reference/missingness_per_condition_cumsum.md),
[`plot_estimate()`](https://fgcz.github.io/prolfqua/reference/plot_estimate.md),
[`plot_heatmap()`](https://fgcz.github.io/prolfqua/reference/plot_heatmap.md),
[`plot_heatmap_cor()`](https://fgcz.github.io/prolfqua/reference/plot_heatmap_cor.md),
[`plot_hierarchies_boxplot_df()`](https://fgcz.github.io/prolfqua/reference/plot_hierarchies_boxplot_df.md),
[`plot_hierarchies_line()`](https://fgcz.github.io/prolfqua/reference/plot_hierarchies_line.md),
[`plot_intensity_distribution_violin()`](https://fgcz.github.io/prolfqua/reference/plot_intensity_distribution_violin.md),
[`plot_na_heatmap()`](https://fgcz.github.io/prolfqua/reference/plot_na_heatmap.md),
[`plot_pca()`](https://fgcz.github.io/prolfqua/reference/plot_pca.md),
[`plot_raster()`](https://fgcz.github.io/prolfqua/reference/plot_raster.md),
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
lfq <- lfq$get_Transformer()$log2()$lfq
#> Column added : log2_abundance
plot_sample_correlation(lfq$data_wide(as.matrix = TRUE)$data)
```
