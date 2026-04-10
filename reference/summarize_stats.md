# Compute mean, sd, and CV for all Peptides, or proteins, for all interactions and all samples.

Compute mean, sd, and CV for all Peptides, or proteins, for all
interactions and all samples.

Compute mean, sd, and CV for e.g. Peptides, or proteins, for all
samples.

summarize stats output (compute quantiles)

## Usage

``` r
summarize_stats(lfqdata, factor_key = lfqdata$relevant_factor_keys())

summarize_stats_all(lfqdata)

summarize_stats_quantiles(
  stats_res,
  factor_keys_depth,
  stats = c("sd", "CV"),
  probs = c(0.1, 0.25, 0.5, 0.75, 0.9)
)
```

## Arguments

- lfqdata:

  LFQData object

- factor_key:

  character vector — factor columns to group by (default:
  relevant_factor_keys)

- stats_res:

  result of running \`summarize_stats\`

- factor_keys_depth:

  character vector — factor columns at current depth

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
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq <- LFQData$new(bb$data, bb$config)
res1 <- summarize_stats(lfq)

res2 <- prolfqua::sim_lfq_data_2factor_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
res2$config$factor_depth <- 2
lfq2 <- LFQData$new(res2$data, res2$config)
stats <- summarize_stats(lfq2)
stopifnot(nrow(stats) == 40)

stats <- summarize_stats(lfq2, factor_key = lfq2$factor_keys()[1])
stopifnot(nrow(stats) == 20)
stats <- summarize_stats(lfq2, factor_key = lfq2$factor_keys()[2])
stopifnot(nrow(stats) == 20)
stats <- summarize_stats(lfq2, factor_key = NULL)
stopifnot(nrow(stats) == 10)

bb <- prolfqua::sim_lfq_data_protein_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq <- LFQData$new(bb$data, bb$config)
res1 <- summarize_stats_all(lfq)
stopifnot((res1 |> dplyr::filter(group_ == "All") |> nrow()) == (res1 |> nrow()))
library(ggplot2)
bb1 <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq <- LFQData$new(bb1$data, bb1$config)
stats_res <- summarize_stats(lfq)
sq <- summarize_stats_quantiles(stats_res, lfq$relevant_factor_keys())
sq <- summarize_stats_quantiles(stats_res, lfq$relevant_factor_keys(), stats = "CV")
sq <- summarize_stats_quantiles(stats_res, lfq$relevant_factor_keys(), stats = "sd")
xx <- summarize_stats_quantiles(stats_res, lfq$relevant_factor_keys(), probs = seq(0, 1, by = 0.1))
ggplot2::ggplot(xx$long, aes(x = probs, y = quantiles, color = group_)) + geom_line() + geom_point()

```
