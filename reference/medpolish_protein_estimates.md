# median polish from normalized peptide intensities

median polish from normalized peptide intensities

## Usage

``` r
medpolish_protein_estimates(data, config)
```

## See also

Other aggregation:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`aggregate_intensity_topN()`](https://wolski.github.io/prolfqua/reference/aggregate_intensity_topN.md),
[`estimate_intensity()`](https://wolski.github.io/prolfqua/reference/estimate_intensity.md),
[`intensity_summary_by_hkeys()`](https://wolski.github.io/prolfqua/reference/intensity_summary_by_hkeys.md),
[`medpolish_estimate()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate.md),
[`medpolish_estimate_df()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate_df.md),
[`medpolish_estimate_dfconfig()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate_dfconfig.md),
[`plot_estimate()`](https://wolski.github.io/prolfqua/reference/plot_estimate.md),
[`plot_hierarchies_add_quantline()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_add_quantline.md),
[`plot_hierarchies_line()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line.md),
[`plot_hierarchies_line_df()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line_df.md),
[`rlm_estimate()`](https://wolski.github.io/prolfqua/reference/rlm_estimate.md),
[`rlm_estimate_dfconfig()`](https://wolski.github.io/prolfqua/reference/rlm_estimate_dfconfig.md)

Other deprecated:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`intensity_summary_by_hkeys()`](https://wolski.github.io/prolfqua/reference/intensity_summary_by_hkeys.md)

## Examples

``` r
istar <- prolfqua_data("data_ionstar")$normalized()
istar$config <- old2new(istar$config)
istar_data <- istar$data
res <- medpolish_protein_estimates(
  istar_data,
  istar$config
)
#> starting aggregation
#> Warning: medpolish() did not converge in 10 iterations
#> Warning: medpolish() did not converge in 10 iterations

dr <- res("unnest")$data
stopifnot(nrow(dr) ==
  length(unique(istar$data$protein_Id)) * length(unique(istar$data$raw.file)))
```
