# plot stddev vs mean to asses stability of variance

plot stddev vs mean to asses stability of variance

## Usage

``` r
plot_stdv_vs_mean(pdata, factor_keys_depth, size = 2000)
```

## Arguments

- pdata:

  data.frame with statistics

- factor_keys_depth:

  character vector — factor columns for faceting

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
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq <- LFQData$new(bb1$data, bb1$config)
res <- lfq$get_Stats()$stats()
plot_stdv_vs_mean(res, lfq$relevant_factor_keys())
#> `geom_smooth()` using formula = 'y ~ x'
#> Warning: Removed 2 rows containing non-finite outside the scale range (`stat_smooth()`).
#> Warning: Removed 2 rows containing missing values or values outside the scale range
#> (`geom_point()`).

```
