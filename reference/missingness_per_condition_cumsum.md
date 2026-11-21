# Cumulative counts of missing values per condition

Cumulative counts of missing values per condition

## Usage

``` r
missingness_per_condition_cumsum(
  lfqdata,
  factors = lfqdata$relevant_factor_keys()
)
```

## Arguments

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

- factors:

  factor columns to use

## Value

a list with `data`, a wide tibble with one row per condition and one
column per number of missing values holding the cumulative feature
count, and `figure`, a ggplot bar chart of these counts

## See also

Other plotting:
[`ContrastsPlotter`](https://fgcz.github.io/prolfqua/reference/ContrastsPlotter.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`missingness_histogram()`](https://fgcz.github.io/prolfqua/reference/missingness_histogram.md),
[`missingness_per_condition()`](https://fgcz.github.io/prolfqua/reference/missingness_per_condition.md),
[`plot_estimate()`](https://fgcz.github.io/prolfqua/reference/plot_estimate.md),
[`plot_heatmap()`](https://fgcz.github.io/prolfqua/reference/plot_heatmap.md),
[`plot_heatmap_cor()`](https://fgcz.github.io/prolfqua/reference/plot_heatmap_cor.md),
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

Other imputation:
[`get_contrast()`](https://fgcz.github.io/prolfqua/reference/get_contrast.md),
[`missingness_histogram()`](https://fgcz.github.io/prolfqua/reference/missingness_histogram.md),
[`missingness_per_condition()`](https://fgcz.github.io/prolfqua/reference/missingness_per_condition.md),
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
res <- missingness_per_condition_cumsum(lfq)
stopifnot("ggplot" %in% class(res$figure))
stopifnot(ncol(res$data) >= 6)
```
