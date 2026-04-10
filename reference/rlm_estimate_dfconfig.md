# Estimate protein abundance from peptide abundances using MASS::rlm

Estimate protein abundance from peptide abundances using MASS::rlm

## Usage

``` r
rlm_estimate_dfconfig(
  pdata,
  response,
  hierarchy_keys,
  hierarchy_keys_depth,
  sample_name,
  name = FALSE
)
```

## Arguments

- pdata:

  data.frame

- response:

  character — intensity column name

- hierarchy_keys:

  character vector — all hierarchy column names

- hierarchy_keys_depth:

  character vector — hierarchy columns at current depth

- sample_name:

  character — sample name column

- name:

  if TRUE return function name only

## See also

[`rlm_estimate`](https://wolski.github.io/prolfqua/reference/rlm_estimate.md)

Other aggregation:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`aggregate_intensity_top_n()`](https://wolski.github.io/prolfqua/reference/aggregate_intensity_top_n.md),
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

feature <- base::setdiff(conf$hierarchy_keys(), conf$hierarchy_keys_depth())
x <- xnested$data[[1]]
bb <- rlm_estimate_dfconfig(x, conf$get_response(),
  conf$hierarchy_keys(), conf$hierarchy_keys_depth(), conf$sample_name)
prolfqua:::.reestablish_condition(x, bb, conf$sample_name,
  conf$factor_keys(), conf$file_name, conf$isotope_label)
#> # A tibble: 12 × 7
#>    sampleName group_ sample  isotopeLabel mean.abundance lmrob weights
#>    <chr>      <chr>  <chr>   <chr>                 <dbl> <dbl>   <dbl>
#>  1 A_V1       A      A_V1    light                  30.0  30.0       1
#>  2 A_V2       A      A_V2    light                  29.8  29.8       1
#>  3 A_V3       A      A_V3    light                  29.1  29.1       1
#>  4 A_V4       A      A_V4    light                  29.9  29.9       1
#>  5 B_V1       B      B_V1    light                  32.1  32.1       1
#>  6 B_V2       B      B_V2    light                  32.4  32.4       1
#>  7 B_V3       B      B_V3    light                  31.7  31.7       1
#>  8 B_V4       B      B_V4    light                  33.5  33.5       1
#>  9 Ctrl_V1    Ctrl   Ctrl_V1 light                  22.1  22.1       1
#> 10 Ctrl_V2    Ctrl   Ctrl_V2 light                  18.5  18.5       1
#> 11 Ctrl_V3    Ctrl   Ctrl_V3 light                  19.8  19.8       1
#> 12 Ctrl_V4    Ctrl   Ctrl_V4 light                  21.4  21.4       1
```
