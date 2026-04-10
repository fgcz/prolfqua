# Compute theoretical sample sizes from factor level standard deviations

Compute theoretical sample sizes from factor level standard deviations

## Usage

``` r
lfq_power_t_test_quantiles(
  lfqdata,
  delta = 1,
  power = 0.8,
  sig.level = 0.05,
  probs = seq(0.5, 0.9, by = 0.1)
)
```

## Arguments

- lfqdata:

  LFQData object

- delta:

  effect size you are interested in

- power:

  of test

- sig.level:

  P-Value

- probs:

  numeric vector of quantile probabilities

## See also

Other stats:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`lfq_power_t_test_proteins()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_proteins.md),
[`lfq_power_t_test_quantiles_V2()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_quantiles_V2.md),
[`plot_stat_density()`](https://wolski.github.io/prolfqua/reference/plot_stat_density.md),
[`plot_stat_density_median()`](https://wolski.github.io/prolfqua/reference/plot_stat_density_median.md),
[`plot_stat_violin()`](https://wolski.github.io/prolfqua/reference/plot_stat_violin.md),
[`plot_stat_violin_median()`](https://wolski.github.io/prolfqua/reference/plot_stat_violin_median.md),
[`plot_stdv_vs_mean()`](https://wolski.github.io/prolfqua/reference/plot_stdv_vs_mean.md),
[`pooled_V2()`](https://wolski.github.io/prolfqua/reference/pooled_var.md),
[`summarize_stats()`](https://wolski.github.io/prolfqua/reference/summarize_stats.md)

## Examples

``` r
bb1 <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq <- LFQData$new(bb1$data, bb1$config)
res <- lfq_power_t_test_quantiles(lfq)
#> Warning: Intensities are not transformed yet.
res$summary
#> # A tibble: 5 × 4
#>   quantile probs    sd `FC=2`
#>   <chr>    <dbl> <dbl>  <dbl>
#> 1 50%        0.5 0.861     13
#> 2 60%        0.6 1.00      17
#> 3 70%        0.7 1.08      20
#> 4 80%        0.8 1.28      27
#> 5 90%        0.9 1.48      36
res <- lfq_power_t_test_quantiles(lfq, delta = 2)
#> Warning: Intensities are not transformed yet.
res <- lfq_power_t_test_quantiles(lfq, delta = c(0.5, 1, 2))
#> Warning: Intensities are not transformed yet.
```
