# Histogram summarizing missigness

Histogram summarizing missigness

## Usage

``` r
missigness_histogram(
  x,
  config,
  showempty = FALSE,
  factors = config$factor_keys_depth(),
  alpha = 0.1
)
```

## See also

Other plotting:
[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`UpSet_interaction_missing_stats()`](https://wolski.github.io/prolfqua/reference/UpSet_interaction_missing_stats.md),
[`UpSet_missing_stats()`](https://wolski.github.io/prolfqua/reference/UpSet_missing_stats.md),
[`medpolish_estimate_df()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate_df.md),
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
[`plot_pca()`](https://wolski.github.io/prolfqua/reference/plot_pca.md),
[`plot_raster()`](https://wolski.github.io/prolfqua/reference/plot_raster.md),
[`plot_sample_correlation()`](https://wolski.github.io/prolfqua/reference/plot_sample_correlation.md)

Other imputation:
[`UpSet_interaction_missing_stats()`](https://wolski.github.io/prolfqua/reference/UpSet_interaction_missing_stats.md),
[`UpSet_missing_stats()`](https://wolski.github.io/prolfqua/reference/UpSet_missing_stats.md),
[`get_contrast()`](https://wolski.github.io/prolfqua/reference/get_contrast.md),
[`missingness_per_condition()`](https://wolski.github.io/prolfqua/reference/missingness_per_condition.md),
[`missingness_per_condition_cumsum()`](https://wolski.github.io/prolfqua/reference/missingness_per_condition_cumsum.md)

## Examples

``` r
istar <- sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
config <- istar$config
analysis <- istar$data
xx <- complete_cases(analysis, config)
#> completing cases
pl <- missigness_histogram(xx, config)
#> Warning: >>>> deprecated! <<<< 
#> 
#>           use summarize_stats_factors instead.
#> completing cases
#> isotopeLabel ~ group_

pl <- missigness_histogram(analysis, config, showempty=FALSE)
#> Warning: >>>> deprecated! <<<< 
#> 
#>           use summarize_stats_factors instead.
#> completing cases
#> isotopeLabel ~ group_
stopifnot("ggplot" %in% class(pl))
pl <- missigness_histogram(analysis, config, showempty=TRUE)
#> Warning: >>>> deprecated! <<<< 
#> 
#>           use summarize_stats_factors instead.
#> completing cases
#> isotopeLabel ~ group_
stopifnot("ggplot" %in% class(pl))
```
