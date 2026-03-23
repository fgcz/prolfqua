# Estimate protein abundance from peptide abundances using MASS::rlm

Estimate protein abundance from peptide abundances using MASS::rlm

## Usage

``` r
rlm_estimate_dfconfig(pdata, config, name = FALSE)
```

## Arguments

- pdata:

  data.frame

- config:

  [`AnalysisConfiguration`](https://wolski.github.io/prolfqua/reference/AnalysisConfiguration.md)

## See also

[`rlm_estimate`](https://wolski.github.io/prolfqua/reference/rlm_estimate.md)

Other aggregation:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`aggregate_intensity_topN()`](https://wolski.github.io/prolfqua/reference/aggregate_intensity_topN.md),
[`estimate_intensity()`](https://wolski.github.io/prolfqua/reference/estimate_intensity.md),
[`medpolish_estimate()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate.md),
[`medpolish_estimate_df()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate_df.md),
[`medpolish_estimate_dfconfig()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate_dfconfig.md),
[`plot_estimate()`](https://wolski.github.io/prolfqua/reference/plot_estimate.md),
[`plot_hierarchies_add_quantline()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_add_quantline.md),
[`plot_hierarchies_line()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line.md),
[`plot_hierarchies_line_df()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line_df.md),
[`rlm_estimate()`](https://wolski.github.io/prolfqua/reference/rlm_estimate.md)

## Examples

``` r

bb <- prolfqua_data("data_ionstar")$filtered()
#> Column added : nr_peptide_Id_IN_protein_Id
conf <- old2new(bb$config)
data <- bb$data
conf$hierarchyDepth <- 1
xnested <- data |>
  dplyr::group_by(across(all_of(conf$hierarchy_keys_depth()))) |>
  tidyr::nest()

feature <- base::setdiff(conf$hierarchy_keys(), conf$hierarchy_keys_depth())
x <- xnested$data[[1]]
bb <- rlm_estimate_dfconfig(x, conf)

prolfqua:::.reestablish_condition(x, bb, conf)
#> # A tibble: 20 × 8
#>    sampleName dilution. run_Id raw.file  isotope mean.peptide.intensity  weights
#>    <chr>      <chr>     <chr>  <chr>     <chr>                    <dbl>    <dbl>
#>  1 b~02       b         02     b03_02_1… light                70845993. 3.98e-14
#>  2 c~03       c         03     b03_03_1… light                69363043. 1.21e-14
#>  3 d~04       d         04     b03_04_1… light                67994529. 1.49e-14
#>  4 e~05       e         05     b03_05_1… light                62812643. 1.16e-14
#>  5 e~06       e         06     b03_06_1… light                67705308. 1.80e-14
#>  6 d~07       d         07     b03_07_1… light                54308813. 1.00e-14
#>  7 c~08       c         08     b03_08_1… light                61674100  8.96e-15
#>  8 b~09       b         09     b03_09_1… light                70944417. 3.09e-14
#>  9 a~10       a         10     b03_10_1… light                66899623. 8.27e-15
#> 10 a~11       a         11     b03_11_1… light                75611447. 1.13e-14
#> 11 b~12       b         12     b03_12_1… light                79307067. 6.70e-15
#> 12 c~13       c         13     b03_13_1… light                67809215. 3.96e-14
#> 13 d~14       d         14     b03_14_1… light                78312000  2.70e-14
#> 14 e~15       e         15     b03_15_1… light                77760429. 3.99e-14
#> 15 e~16       e         16     b03_16_1… light                70743500  7.85e-15
#> 16 d~17       d         17     b03_17_1… light                81735462. 9.39e-15
#> 17 c~18       c         18     b03_18_1… light                74486688. 3.50e-15
#> 18 b~19       b         19     b03_19_1… light                74329840  1.19e-14
#> 19 a~20       a         20     b03_20_1… light                73915533. 5.14e-15
#> 20 a~21       a         21     b03_21_1… light                70448067. 7.99e-15
#> # ℹ 1 more variable: lmrob <dbl>
```
