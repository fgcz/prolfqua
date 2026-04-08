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
[`aggregate_intensity_topN()`](https://wolski.github.io/prolfqua/reference/aggregate_intensity_topN.md),
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
bb <- prolfqua_data("data_ionstar")$filtered()
#> Column added : nr_peptide_Id_IN_protein_Id
stopifnot(nrow(bb$data) == 25780)
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
prolfqua:::.reestablish_condition(x, bb, conf)
#> # A tibble: 20 × 6
#>    sampleName dilution. run_Id raw.file                        isotope medpolish
#>    <chr>      <chr>     <chr>  <chr>                           <chr>       <dbl>
#>  1 b~02       b         02     b03_02_150304_human_ecoli_b_3u… light   53281919.
#>  2 c~03       c         03     b03_03_150304_human_ecoli_c_3u… light   53108622.
#>  3 d~04       d         04     b03_04_150304_human_ecoli_d_3u… light   51768575 
#>  4 e~05       e         05     b03_05_150304_human_ecoli_e_3u… light   48923212.
#>  5 e~06       e         06     b03_06_150304_human_ecoli_e_3u… light   46757962.
#>  6 d~07       d         07     b03_07_150304_human_ecoli_d_3u… light   43139562.
#>  7 c~08       c         08     b03_08_150304_human_ecoli_c_3u… light   49301991.
#>  8 b~09       b         09     b03_09_150304_human_ecoli_b_3u… light   46173388.
#>  9 a~10       a         10     b03_10_150304_human_ecoli_a_3u… light   47789438.
#> 10 a~11       a         11     b03_11_150304_human_ecoli_a_3u… light   61505569.
#> 11 b~12       b         12     b03_12_150304_human_ecoli_b_3u… light   57534412.
#> 12 c~13       c         13     b03_13_150304_human_ecoli_c_3u… light   50100838.
#> 13 d~14       d         14     b03_14_150304_human_ecoli_d_3u… light   58970725 
#> 14 e~15       e         15     b03_15_150304_human_ecoli_e_3u… light   61280338.
#> 15 e~16       e         16     b03_16_150304_human_ecoli_e_3u… light   56780612.
#> 16 d~17       d         17     b03_17_150304_human_ecoli_d_3u… light   59887225 
#> 17 c~18       c         18     b03_18_150304_human_ecoli_c_3u… light   57788409.
#> 18 b~19       b         19     b03_19_150304_human_ecoli_b_3u… light   60040838.
#> 19 a~20       a         20     b03_20_150304_human_ecoli_a_3u… light   57032150 
#> 20 a~21       a         21     b03_21_150304_human_ecoli_a_3u… light   55993412.
```
