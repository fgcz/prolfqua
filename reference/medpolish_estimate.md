# Median polish estimate of e.g. protein from peptide intensities

Compute Tukey's median polish estimate of a protein from peptide or
precursor intensities

## Usage

``` r
medpolish_estimate(x, name = FALSE, sample_name = "sample_name")
```

## Arguments

- x:

  a matrix

- name:

  if TRUE returns the name of the summary column

## Value

data.frame with number of rows equal to number of columns of input
matrix.

## See also

Other aggregation:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`aggregate_intensity_top_n()`](https://wolski.github.io/prolfqua/reference/aggregate_intensity_top_n.md),
[`estimate_intensity()`](https://wolski.github.io/prolfqua/reference/estimate_intensity.md),
[`medpolish_estimate_df()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate_df.md),
[`medpolish_estimate_dfconfig()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate_dfconfig.md),
[`plot_estimate()`](https://wolski.github.io/prolfqua/reference/plot_estimate.md),
[`plot_hierarchies_add_quantline()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_add_quantline.md),
[`plot_hierarchies_line()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line.md),
[`plot_hierarchies_line_df()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line_df.md),
[`rlm_estimate()`](https://wolski.github.io/prolfqua/reference/rlm_estimate.md),
[`rlm_estimate_dfconfig()`](https://wolski.github.io/prolfqua/reference/rlm_estimate_dfconfig.md)

## Examples

``` r
medpolish_estimate(name = TRUE)
#> [1] "medpolish"
gg <- matrix(runif(20), 4, 5)
rownames(gg) <- paste0("A", 1:4)
colnames(gg) <- make.names(1:5)
gg
#>           X1        X2        X3        X4        X5
#> A1 0.3170962 0.3568517 0.2082721 0.3275074 0.8280330
#> A2 0.4747718 0.5299904 0.2738367 0.5142372 0.9245008
#> A3 0.5459187 0.5604496 0.7831665 0.4032164 0.4215357
#> A4 0.9429698 0.3253429 0.5733794 0.9108078 0.2704623
mx <- medpolish_estimate(gg)
```
