# Extract response column of a protein into matrix Median polish estimates of e.g. protein abundances for entire data.frame

Extract response column of a protein into matrix Median polish estimates
of e.g. protein abundances for entire data.frame

## Usage

``` r
medpolish_estimate_df(pdata, response, feature, sample_name)
```

## Arguments

- pdata:

  data.frame

- response:

  column name with intensities

- feature:

  column name e.g. peptide ids

- sample_name:

  column name e.g. sample_name

## Value

data.frame

## See also

[`medpolish_estimate`](https://wolski.github.io/prolfqua/reference/medpolish_estimate.md)

Other aggregation:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`aggregate_intensity_top_n()`](https://wolski.github.io/prolfqua/reference/aggregate_intensity_top_n.md),
[`estimate_intensity()`](https://wolski.github.io/prolfqua/reference/estimate_intensity.md),
[`medpolish_estimate()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate.md),
[`medpolish_estimate_dfconfig()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate_dfconfig.md),
[`plot_estimate()`](https://wolski.github.io/prolfqua/reference/plot_estimate.md),
[`plot_hierarchies_add_quantline()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_add_quantline.md),
[`plot_hierarchies_line()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line.md),
[`plot_hierarchies_line_df()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line_df.md),
[`rlm_estimate()`](https://wolski.github.io/prolfqua/reference/rlm_estimate.md),
[`rlm_estimate_dfconfig()`](https://wolski.github.io/prolfqua/reference/rlm_estimate_dfconfig.md)

Other plotting:
[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
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
[`plot_pca()`](https://wolski.github.io/prolfqua/reference/plot_pca.md),
[`plot_raster()`](https://wolski.github.io/prolfqua/reference/plot_raster.md),
[`plot_sample_correlation()`](https://wolski.github.io/prolfqua/reference/plot_sample_correlation.md),
[`upset_interaction_missing_stats()`](https://wolski.github.io/prolfqua/reference/UpSet_interaction_missing_stats.md),
[`upset_missing_stats()`](https://wolski.github.io/prolfqua/reference/UpSet_missing_stats.md)

## Examples

``` r
bb <- sim_lfq_data_peptide_config(Nprot = 20)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
conf <- bb$config
data <- bb$data

conf$hierarchy_depth <- 1
xnested <- data |>
  dplyr::group_by(across(all_of(conf$hierarchy_keys_depth()))) |>
  tidyr::nest()

feature <- base::setdiff(
  conf$hierarchy_keys(),
  conf$hierarchy_keys_depth()
)
x <- xnested$data[[1]]
bb <- medpolish_estimate_df(x,
  response = conf$get_response(),
  feature = feature,
  sample_name = conf$sample_name
)
prolfqua:::.reestablish_condition(x, bb, conf$sample_name,
  conf$factor_keys(), conf$file_name, conf$isotope_label)
#> # A tibble: 12 × 5
#>    sampleName group_ sample  isotopeLabel medpolish
#>    <chr>      <chr>  <chr>   <chr>            <dbl>
#>  1 A_V1       A      A_V1    light             30.0
#>  2 A_V2       A      A_V2    light             29.8
#>  3 A_V3       A      A_V3    light             29.1
#>  4 A_V4       A      A_V4    light             29.9
#>  5 B_V1       B      B_V1    light             32.1
#>  6 B_V2       B      B_V2    light             32.4
#>  7 B_V3       B      B_V3    light             31.7
#>  8 B_V4       B      B_V4    light             33.5
#>  9 Ctrl_V1    Ctrl   Ctrl_V1 light             22.1
#> 10 Ctrl_V2    Ctrl   Ctrl_V2 light             18.5
#> 11 Ctrl_V3    Ctrl   Ctrl_V3 light             19.8
#> 12 Ctrl_V4    Ctrl   Ctrl_V4 light             21.4
```
