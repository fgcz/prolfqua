# Per-test significance level for FDR-controlled sample-size planning

Convert a target false discovery rate and average power into the
per-test significance level used by an ordinary two-sided t-test. Pass
the result to
[`stats::power.t.test()`](https://rdrr.io/r/stats/power.t.test.html),
[`lfq_power_t_test_proteins()`](https://fgcz.github.io/prolfqua/reference/lfq_power_t_test_proteins.md),
or
[`lfq_power_t_test_quantiles()`](https://fgcz.github.io/prolfqua/reference/lfq_power_t_test_quantiles.md).

## Usage

``` r
fdr_sig_level(fdr = 0.05, pi0 = 0.95, power = 0.8)
```

## Arguments

- fdr:

  one target false discovery rate, strictly between zero and one.

- pi0:

  one assumed proportion of truly null features, strictly between zero
  and one; it is not estimated from the pilot.

- power:

  one target average power among non-null features, strictly between
  zero and one.

## Value

The scalar per-test significance level \\\alpha\_\* = power \times fdr
\times (1-pi0) / (pi0 \times (1-fdr))\\. Combinations giving a
significance level outside (0, 1) are rejected.

## Details

The Jung/Liu–Hwang planning relationship assumes independent tests or an
adequate approximation under weak dependence, and many tested features.
It controls the expected false-discovery proportion under the assumed
null fraction and effect distribution, not the realized FDR of every
experiment. All sample sizes are per group in a balanced two-group
design with equal within-group variances.

## References

Jung SH (2005). Sample size for FDR-control in microarray data analysis.
Bioinformatics 21:3097–3104.
[doi:10.1093/bioinformatics/bti456](https://doi.org/10.1093/bioinformatics/bti456)
. Liu P, Hwang JTG (2007). Quick calculation for sample size while
controlling false discovery rate with application to microarray
analysis. Bioinformatics 23:739–746.
[doi:10.1093/bioinformatics/btl664](https://doi.org/10.1093/bioinformatics/btl664)
.

## See also

Other stats:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`compute_pooled()`](https://fgcz.github.io/prolfqua/reference/pooled_var.md),
[`lfq_power_t_test_proteins()`](https://fgcz.github.io/prolfqua/reference/lfq_power_t_test_proteins.md),
[`lfq_power_t_test_quantiles()`](https://fgcz.github.io/prolfqua/reference/lfq_power_t_test_quantiles.md),
[`null_split_check()`](https://fgcz.github.io/prolfqua/reference/null_split_check.md),
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
alpha <- fdr_sig_level(fdr = 0.05, pi0 = 0.9, power = 0.8)
alpha
#> [1] 0.004678363
ceiling(stats::power.t.test(delta = 1, sd = 1, sig.level = alpha, power = 0.8)$n)
#> [1] 29
```
