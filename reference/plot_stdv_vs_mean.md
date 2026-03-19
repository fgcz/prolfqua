# plot stddev vs mean to asses stability of variance

plot stddev vs mean to asses stability of variance

## Usage

``` r
plot_stdv_vs_mean(pdata, config, size = 2000)
```

## Arguments

- pdata:

  data.frame

- config:

  AnalysisConfiguration

- size:

  how many points to sample (since scatter plot to slow for all)

## See also

Other stats:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`lfq_power_t_test_proteins()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_proteins.md),
[`lfq_power_t_test_quantiles()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_quantiles.md),
[`lfq_power_t_test_quantiles_V2()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_quantiles_V2.md),
[`plot_stat_density()`](https://wolski.github.io/prolfqua/reference/plot_stat_density.md),
[`plot_stat_density_median()`](https://wolski.github.io/prolfqua/reference/plot_stat_density_median.md),
[`plot_stat_violin()`](https://wolski.github.io/prolfqua/reference/plot_stat_violin.md),
[`plot_stat_violin_median()`](https://wolski.github.io/prolfqua/reference/plot_stat_violin_median.md),
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
data <- bb1$data
res <- summarize_stats(data, config)
#> completing cases

plot_stdv_vs_mean(res, config)
#> `geom_smooth()` using formula = 'y ~ x'
#> Warning: Removed 2 rows containing non-finite outside the scale range
#> (`stat_smooth()`).
#> Warning: Removed 2 rows containing missing values or values outside
#> the scale range (`geom_point()`).

datalog2 <- transform_work_intensity(data, config, log2)
#> Column added : log2_abundance
statlog2 <- summarize_stats(datalog2, config)
#> completing cases
plot_stdv_vs_mean(statlog2, config)
#> `geom_smooth()` using formula = 'y ~ x'
#> Warning: Removed 2 rows containing non-finite outside the scale range
#> (`stat_smooth()`).
#> Warning: Removed 2 rows containing missing values or values outside
#> the scale range (`geom_point()`).

config$get_response()
#> [1] "log2_abundance"
config$pop_response()
#> [1] "log2_abundance"
datasqrt <- transform_work_intensity(data, config, sqrt)
#> Column added : sqrt_abundance
ressqrt <- summarize_stats(datasqrt, config)
#> completing cases
plot_stdv_vs_mean(ressqrt, config)
#> `geom_smooth()` using formula = 'y ~ x'
#> Warning: Removed 2 rows containing non-finite outside the scale range
#> (`stat_smooth()`).
#> Warning: Removed 2 rows containing missing values or values outside
#> the scale range (`geom_point()`).

```
