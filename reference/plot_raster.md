# plot heatmap without any clustering (use to show NA's)

plot heatmap without any clustering (use to show NA's)

## Usage

``` r
plot_raster(
  data,
  config,
  arrange = c("mean", "var"),
  not_na = FALSE,
  show_rownames = FALSE,
  ...
)
```

## Arguments

- data:

  dataframe

- config:

  dataframe configuration

- arrange:

  either mean or var

- not_na:

  if true than arrange by nr of NA's first and then by arrange

- show_rownames:

  logical, show row names in heatmap

- ...:

  additional arguments passed to pheatmap

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
[`plot_pca()`](https://wolski.github.io/prolfqua/reference/plot_pca.md),
[`plot_sample_correlation()`](https://wolski.github.io/prolfqua/reference/plot_sample_correlation.md)

## Examples

``` r

istar <- sim_lfq_data_protein_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
config <- istar$config
analysis <- istar$data
rs <- plot_raster(analysis, config, show_rownames=FALSE)
stopifnot(class(rs) == "pheatmap")
rs <- plot_raster(analysis[1,], config)
#> Warning: The dataset has :1
stopifnot(is.null(rs))
rs <- plot_raster(analysis, config, "var")
stopifnot(class(rs) == "pheatmap")
rs <- plot_raster(analysis, config, show_rownames = TRUE)
stopifnot(class(rs) == "pheatmap")
```
