# Violin plot of sd, CV or mean

Violin plot of sd, CV or mean

## Usage

``` r
plot_stat_violin(
  pdata,
  factor_keys_depth,
  stat = c("CV", "meanAbundance", "sd")
)
```

## Arguments

- pdata:

  data.frame with statistics

- factor_keys_depth:

  character vector - factor columns for grouping

- stat:

  either CV, mean or sd

## Value

a ggplot object with one violin of `stat` per combination of factor
levels, with the median marked

## See also

Other stats:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`compute_pooled()`](https://fgcz.github.io/prolfqua/reference/pooled_var.md),
[`fdr_sig_level()`](https://fgcz.github.io/prolfqua/reference/fdr_sig_level.md),
[`lfq_power_t_test_proteins()`](https://fgcz.github.io/prolfqua/reference/lfq_power_t_test_proteins.md),
[`lfq_power_t_test_quantiles()`](https://fgcz.github.io/prolfqua/reference/lfq_power_t_test_quantiles.md),
[`null_split_check()`](https://fgcz.github.io/prolfqua/reference/null_split_check.md),
[`plot_stat_density()`](https://fgcz.github.io/prolfqua/reference/plot_stat_density.md),
[`plot_stat_density_median()`](https://fgcz.github.io/prolfqua/reference/plot_stat_density_median.md),
[`plot_stat_violin_median()`](https://fgcz.github.io/prolfqua/reference/plot_stat_violin_median.md),
[`plot_stdv_vs_mean()`](https://fgcz.github.io/prolfqua/reference/plot_stdv_vs_mean.md),
[`pooled_variance_explained()`](https://fgcz.github.io/prolfqua/reference/pooled_variance_explained.md),
[`power_fdr_prior_n()`](https://fgcz.github.io/prolfqua/reference/power_fdr_prior_n.md),
[`summarize_stats()`](https://fgcz.github.io/prolfqua/reference/summarize_stats.md),
[`summarize_variance_explained()`](https://fgcz.github.io/prolfqua/reference/summarize_variance_explained.md)

## Examples

``` r
bb1 <- prolfqua::sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
lfq <- LFQData$new(bb1$data, bb1$config)
res <- lfq$get_Stats()$stats()
plot_stat_violin(res, lfq$relevant_factor_keys(), stat = "meanAbundance")

plot_stat_violin(res, lfq$relevant_factor_keys(), stat = "sd")
#> Warning: Removed 2 rows containing non-finite outside the scale range
#> (`stat_ydensity()`).
#> Warning: Removed 2 rows containing non-finite outside the scale range
#> (`stat_summary()`).

plot_stat_violin(res, lfq$relevant_factor_keys(), stat = "CV")
#> Warning: Removed 2 rows containing non-finite outside the scale range
#> (`stat_ydensity()`).
#> Warning: Removed 2 rows containing non-finite outside the scale range
#> (`stat_summary()`).

```
