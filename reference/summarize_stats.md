# Compute mean, sd, and CV for all Peptides, or proteins, for all interactions and all samples.

Compute mean, sd, and CV for all Peptides, or proteins, for all
interactions and all samples.

Compute mean, sd, and CV for e.g. Peptides, or proteins, for all
samples.

summarize stats output (compute quantiles)

## Usage

``` r
summarize_stats(
  pdata,
  config,
  factor_key = config$factor_keys_depth(),
  .completed = FALSE
)

summarize_stats_all(pdata, config, .completed = FALSE)

summarize_stats_quantiles(
  stats_res,
  config,
  stats = c("sd", "CV"),
  probs = c(0.1, 0.25, 0.5, 0.75, 0.9)
)
```

## Arguments

- pdata:

  data.frame

- config:

  AnalysisConfiguration

- stats_res:

  result of running \`summarize_stats\`

- stats:

  summarize either sd or CV

- probs:

  for which quantiles 10, 20 etc.

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
[`plot_stdv_vs_mean()`](https://wolski.github.io/prolfqua/reference/plot_stdv_vs_mean.md),
[`pooled_V2()`](https://wolski.github.io/prolfqua/reference/pooled_var.md)

Other stats:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`lfq_power_t_test_proteins()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_proteins.md),
[`lfq_power_t_test_quantiles()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_quantiles.md),
[`lfq_power_t_test_quantiles_V2()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_quantiles_V2.md),
[`plot_stat_density()`](https://wolski.github.io/prolfqua/reference/plot_stat_density.md),
[`plot_stat_density_median()`](https://wolski.github.io/prolfqua/reference/plot_stat_density_median.md),
[`plot_stat_violin()`](https://wolski.github.io/prolfqua/reference/plot_stat_violin.md),
[`plot_stat_violin_median()`](https://wolski.github.io/prolfqua/reference/plot_stat_violin_median.md),
[`plot_stdv_vs_mean()`](https://wolski.github.io/prolfqua/reference/plot_stdv_vs_mean.md),
[`pooled_V2()`](https://wolski.github.io/prolfqua/reference/pooled_var.md)

Other stats:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`lfq_power_t_test_proteins()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_proteins.md),
[`lfq_power_t_test_quantiles()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_quantiles.md),
[`lfq_power_t_test_quantiles_V2()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_quantiles_V2.md),
[`plot_stat_density()`](https://wolski.github.io/prolfqua/reference/plot_stat_density.md),
[`plot_stat_density_median()`](https://wolski.github.io/prolfqua/reference/plot_stat_density_median.md),
[`plot_stat_violin()`](https://wolski.github.io/prolfqua/reference/plot_stat_violin.md),
[`plot_stat_violin_median()`](https://wolski.github.io/prolfqua/reference/plot_stat_violin_median.md),
[`plot_stdv_vs_mean()`](https://wolski.github.io/prolfqua/reference/plot_stdv_vs_mean.md),
[`pooled_V2()`](https://wolski.github.io/prolfqua/reference/pooled_var.md)

## Examples

``` r

bb <- prolfqua::sim_lfq_data_protein_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
config <- bb$config
data <- bb$data

res1 <- summarize_stats(data, config)
#> completing cases

res2 <- prolfqua::sim_lfq_data_2Factor_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
res2$config$factorDepth <- 2
stats <- summarize_stats(res2$data, res2$config)
#> completing cases
stopifnot(nrow(stats) == 40)

stats <- summarize_stats(res2$data, res2$config, factor_key = res2$config$factor_keys()[1])
#> completing cases
stopifnot(nrow(stats) == 20)
stats <- summarize_stats(res2$data, res2$config, factor_key = res2$config$factor_keys()[2])
#> completing cases
stopifnot(nrow(stats) == 20)
stats <- summarize_stats(res2$data, res2$config, factor_key = NULL)
#> completing cases
stopifnot(nrow(stats) == 10)
# TODO (WEW) add test when there is one level per group.


bb <- prolfqua::sim_lfq_data_protein_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done

res1 <- summarize_stats_all(bb$data, bb$config)
#> completing cases

stopifnot((res1 |> dplyr::filter(group_ == "All") |> nrow()) == (res1 |> nrow()))
res2 <- prolfqua::sim_lfq_data_2Factor_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
resSt <- summarize_stats_all(res2$data, res2$config)
#> completing cases
library(ggplot2)
bb1 <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
config <- bb1$config
data <- bb1$data
stats_res <- summarize_stats(data, config)
#> completing cases
sq <- summarize_stats_quantiles(stats_res, config)
sq <- summarize_stats_quantiles(stats_res, config, stats = "CV")
bb <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
config <- bb$config
data <- bb$data
config$get_response()
#> [1] "abundance"
stats_res <- summarize_stats(data, config)
#> completing cases
sq <- summarize_stats_quantiles(stats_res, config)
sq <- summarize_stats_quantiles(stats_res, config, stats = "sd")
stats_res <- summarize_stats(data, config)
#> completing cases
xx <- summarize_stats_quantiles(stats_res, config, probs = seq(0,1,by = 0.1))
ggplot2::ggplot(xx$long, aes(x = probs, y = quantiles, color = group_)) + geom_line() + geom_point()


```
