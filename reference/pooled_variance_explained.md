# Variance explained by each term, pooled over features

Sums the sums of squares of
[`summarize_variance_explained()`](https://fgcz.github.io/prolfqua/reference/summarize_variance_explained.md)
over features. For each term the result is
`sum(Sum.Sq) / sum(SS_total)`; one minus the share of the residuals is
the pooled R^2. On features without missing values the pooled R^2 equals
the PERMANOVA R^2 with Euclidean distances between samples. Features
with a large total variance weigh more than in the median of the
per-feature `R2`.

## Usage

``` r
pooled_variance_explained(x)
```

## Arguments

- x:

  the tibble returned by
  [`summarize_variance_explained()`](https://fgcz.github.io/prolfqua/reference/summarize_variance_explained.md).

## Value

A data.frame with columns `factor` and `eta2`, one row per term in the
order of the ANOVA table, residuals included.

## See also

[`summarize_variance_explained()`](https://fgcz.github.io/prolfqua/reference/summarize_variance_explained.md)

Other stats:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`compute_pooled()`](https://fgcz.github.io/prolfqua/reference/pooled_var.md),
[`fdr_sig_level()`](https://fgcz.github.io/prolfqua/reference/fdr_sig_level.md),
[`lfq_power_t_test_proteins()`](https://fgcz.github.io/prolfqua/reference/lfq_power_t_test_proteins.md),
[`lfq_power_t_test_quantiles()`](https://fgcz.github.io/prolfqua/reference/lfq_power_t_test_quantiles.md),
[`null_split_check()`](https://fgcz.github.io/prolfqua/reference/null_split_check.md),
[`plot_stat_density()`](https://fgcz.github.io/prolfqua/reference/plot_stat_density.md),
[`plot_stat_density_median()`](https://fgcz.github.io/prolfqua/reference/plot_stat_density_median.md),
[`plot_stat_violin()`](https://fgcz.github.io/prolfqua/reference/plot_stat_violin.md),
[`plot_stat_violin_median()`](https://fgcz.github.io/prolfqua/reference/plot_stat_violin_median.md),
[`plot_stdv_vs_mean()`](https://fgcz.github.io/prolfqua/reference/plot_stdv_vs_mean.md),
[`power_fdr_prior_n()`](https://fgcz.github.io/prolfqua/reference/power_fdr_prior_n.md),
[`summarize_stats()`](https://fgcz.github.io/prolfqua/reference/summarize_stats.md),
[`summarize_variance_explained()`](https://fgcz.github.io/prolfqua/reference/summarize_variance_explained.md)

## Examples

``` r
bb <- prolfqua::sim_lfq_data_2factor_config(n_proteins = 20)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, group, mean, N, Replicate.
#> setup_analysis: added 44 rows with NA response for feature and sample combinations absent from the input.
bb$config$factor_depth <- 2
lfqdata <- LFQData$new(bb$data, bb$config)
lfqdata$is_transformed(TRUE)
pooled <- pooled_variance_explained(summarize_variance_explained(lfqdata))
#> fit (Model): 44 missing responses dropped (e.g. 0EfVhX~5954, 0m5WN4~1448, 7cbcrd~8305, 9VUkAq~4562, At886V~3296, ...): NA response (missing observation) left out of the fit
#> See $provenance().
pooled
#>                 factor      eta2
#> 1            Treatment 0.4180986
#> 2           Background 0.1680539
#> 3 Treatment:Background 0.1393417
#> 4            Residuals 0.2745057
# pooled R^2
1 - pooled$eta2[pooled$factor == "Residuals"]
#> [1] 0.7254943
```
