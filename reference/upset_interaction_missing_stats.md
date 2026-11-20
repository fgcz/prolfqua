# UpSetR plot from interaction_missing_stats

UpSetR plot from interaction_missing_stats

## Usage

``` r
upset_interaction_missing_stats(lfqdata, tr = 2)
```

## Arguments

- lfqdata:

  LFQData object

- tr:

  if less than tr observations in condition then missing

## Value

The computed result.

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
[`plot_sample_correlation()`](https://fgcz.github.io/prolfqua/reference/plot_sample_correlation.md),
[`upset_missing_stats()`](https://fgcz.github.io/prolfqua/reference/upset_missing_stats.md)

Other imputation:
[`get_contrast()`](https://fgcz.github.io/prolfqua/reference/get_contrast.md),
[`missigness_histogram()`](https://fgcz.github.io/prolfqua/reference/missigness_histogram.md),
[`missingness_per_condition()`](https://fgcz.github.io/prolfqua/reference/missingness_per_condition.md),
[`missingness_per_condition_cumsum()`](https://fgcz.github.io/prolfqua/reference/missingness_per_condition_cumsum.md),
[`upset_missing_stats()`](https://fgcz.github.io/prolfqua/reference/upset_missing_stats.md)

## Examples

``` r
istar <- sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
lfq <- LFQData$new(istar$data, istar$config)
pups <- upset_interaction_missing_stats(lfq)
stopifnot(ncol(pups$data) == 5)
UpSetR::upset(pups$data, order.by = "freq", nsets = pups$nsets)
```
