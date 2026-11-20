# Null-split calibration check of normalization and model

Splits the samples within every stratum at random into two
pseudo-groups, A and B, and tests B against A with
[`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md).
No feature differs between the pseudo-groups, so a calibrated pipeline
gives uniform p-values: about 5% below 0.05, 1% below 0.01, and no
feature with FDR \< 0.05. The check uses no labels beyond the strata, so
it compares normalizations without favouring any of them by
construction.

## Usage

``` r
null_split_check(
  lfqdata,
  method = "lm",
  n_splits = 20,
  strata = lfqdata$relevant_factor_keys(),
  n_features = NULL,
  seed = NULL,
  ...
)
```

## Arguments

- lfqdata:

  the model input: a normalized
  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md)
  object, aggregated or not as `method` requires.

- method:

  facade name passed to
  [`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md).

- n_splits:

  number of distinct random splits.

- strata:

  factor keys within whose levels the samples are split; the model is
  `~ <strata> + pseudo_group_`. `character(0)` splits all samples as one
  group; subset `lfqdata` to one condition first.

- n_features:

  number of randomly chosen top-level features (`hierarchy_keys()[1]`)
  to test; `NULL` (default) tests all.

- seed:

  `NULL` (default) uses the current random number generator state;
  otherwise the draws run under
  [`withr::with_seed()`](https://withr.r-lib.org/reference/with_seed.html).

- ...:

  further arguments to
  [`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md).

## Value

A tibble with one row per split: `split`, `n` (tests with a p-value),
`p_below_0.05`, `p_below_0.01` (shares of p-values), `n_FDR_below_0.05`,
`ks_D` (Kolmogorov-Smirnov distance of the p-values from the uniform
distribution), and the list-columns `assignment` (sample name and
`pseudo_group_` per sample) and `contrasts` (subject keys, `diff`,
`p.value` and `FDR`).

## Details

Within each stratum the samples are split into halves; an odd number of
samples gives halves that differ by one. Splits that differ only by
swapping A and B count as one. When the design has fewer distinct splits
than `n_splits`, all of them are returned with a message.

Read the median and range over splits. Splits share samples, so they are
not independent replicates. A split that follows a batch or the run
order puts a real technical effect into the contrast: compare
`assignment` with the batch of each sample, and persistent hits point to
a batch effect the normalization left in. With few samples per
pseudo-group, variance moderation and imputation behave differently than
in the real comparison; calibration carries over, power does not.

## See also

[`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md)

Other stats:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`compute_pooled()`](https://fgcz.github.io/prolfqua/reference/pooled_var.md),
[`fdr_sig_level()`](https://fgcz.github.io/prolfqua/reference/fdr_sig_level.md),
[`lfq_power_t_test_proteins()`](https://fgcz.github.io/prolfqua/reference/lfq_power_t_test_proteins.md),
[`lfq_power_t_test_quantiles_V2()`](https://fgcz.github.io/prolfqua/reference/lfq_power_t_test_quantiles_V2.md),
[`plot_stat_density()`](https://fgcz.github.io/prolfqua/reference/plot_stat_density.md),
[`plot_stat_density_median()`](https://fgcz.github.io/prolfqua/reference/plot_stat_density_median.md),
[`plot_stat_violin()`](https://fgcz.github.io/prolfqua/reference/plot_stat_violin.md),
[`plot_stat_violin_median()`](https://fgcz.github.io/prolfqua/reference/plot_stat_violin_median.md),
[`plot_stdv_vs_mean()`](https://fgcz.github.io/prolfqua/reference/plot_stdv_vs_mean.md),
[`pooled_variance_explained()`](https://fgcz.github.io/prolfqua/reference/pooled_variance_explained.md),
[`power_fdr_prior_n()`](https://fgcz.github.io/prolfqua/reference/power_fdr_prior_n.md),
[`summarize_stats()`](https://fgcz.github.io/prolfqua/reference/summarize_stats.md),
[`summarize_variance_explained()`](https://fgcz.github.io/prolfqua/reference/summarize_variance_explained.md)

## Examples

``` r
bb <- sim_lfq_data_protein_config(Nprot = 50, with_missing = FALSE)
#> creating sampleName from file_name column
lfqdata <- LFQData$new(bb$data, bb$config)
lfqdata$is_transformed(TRUE)
res <- null_split_check(lfqdata, method = "limma", n_splits = 5, seed = 1)
res[, c("split", "n", "p_below_0.05", "n_FDR_below_0.05", "ks_D")]
#> # A tibble: 5 × 5
#>   split     n p_below_0.05 n_FDR_below_0.05   ks_D
#>   <int> <int>        <dbl>            <int>  <dbl>
#> 1     1    50         0.08                0 0.154 
#> 2     2    50         0.04                0 0.132 
#> 3     3    50         0.06                0 0.107 
#> 4     4    50         0.04                0 0.0632
#> 5     5    50         0.02                0 0.0614
```
