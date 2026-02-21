# estimate sample sizes

estimate sample sizes

## Usage

``` r
lfq_power_t_test_quantiles_V2(
  quantile_sd,
  delta = c(0.59, 1, 2),
  power = 0.8,
  sig.level = 0.05,
  min.n = 1.5
)
```

## Arguments

- quantile_sd:

  output of \`summarize_stats_quantiles\`

- delta:

  effect size you are interested in

- power:

  of test

- sig.level:

  P-Value

- min.n:

  smallest n to determine

## See also

Other stats:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`lfq_power_t_test_proteins()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_proteins.md),
[`lfq_power_t_test_quantiles()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_quantiles.md),
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
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
config <- bb1$config
data2 <- bb1$data
stats_res <- summarize_stats(data2, config)
#> completing cases
xx <- summarize_stats_quantiles(stats_res, config, probs = c(0.5,0.8))
bbb <- lfq_power_t_test_quantiles_V2(xx$long)
bbb <- dplyr::bind_rows(bbb)
summary <- bbb |>
 dplyr::select( -N_exact, -quantiles, -sdtrimmed ) |>
 tidyr::pivot_wider(names_from = delta, values_from = N)
```
