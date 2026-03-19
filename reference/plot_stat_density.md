# plot density distribution or ecdf of sd, mean or CV

plot density distribution or ecdf of sd, mean or CV

## Usage

``` r
plot_stat_density(
  pdata,
  config,
  stat = c("CV", "meanAbundance", "sd"),
  ggstat = c("density", "ecdf")
)
```

## Arguments

- pdata:

  data.frame

- config:

  AnalysisConfiguration

- stat:

  sd, mean or CV

- ggstat:

  either density of ecdf

## See also

Other stats:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`lfq_power_t_test_proteins()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_proteins.md),
[`lfq_power_t_test_quantiles()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_quantiles.md),
[`lfq_power_t_test_quantiles_V2()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_quantiles_V2.md),
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
config <-  bb1$config
data <- bb1$data
res <- summarize_stats(data, config)
#> completing cases
plot_stat_density(res, config, stat = "meanAbundance")

plot_stat_density(res, config, stat = "sd")
#> Warning: Removed 2 rows containing non-finite outside the scale range
#> (`stat_density()`).

plot_stat_density(res, config, stat = "CV")
#> Warning: Removed 2 rows containing non-finite outside the scale range
#> (`stat_density()`).
```
