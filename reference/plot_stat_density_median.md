# plot density distribution or ecdf of sd, mean or cv given intensity below and above median

plot density distribution or ecdf of sd, mean or cv given intensity
below and above median

## Usage

``` r
plot_stat_density_median(
  pdata,
  factor_key,
  stat = c("CV", "meanAbundance", "sd"),
  ggstat = c("density", "ecdf")
)
```

## Arguments

- pdata:

  data.frame with statistics

- factor_key:

  character — factor column name for faceting

- stat:

  sd, mean or CV

- ggstat:

  either density or ecdf

## See also

Other stats:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`lfq_power_t_test_proteins()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_proteins.md),
[`lfq_power_t_test_quantiles()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_quantiles.md),
[`lfq_power_t_test_quantiles_V2()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_quantiles_V2.md),
[`plot_stat_density()`](https://wolski.github.io/prolfqua/reference/plot_stat_density.md),
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
res <- lfq$get_Stats()$stats()
plot_stat_density_median(res, lfq$factor_keys()[1], "CV")

plot_stat_density_median(res, lfq$factor_keys()[1], "meanAbundance")

plot_stat_density_median(res, lfq$factor_keys()[1], "sd")
```
