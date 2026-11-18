# Plan sample size from a pilot variance prior and a planned effect distribution

Fit a variance prior to per-feature pilot variances and residual degrees
of freedom, then integrate ordinary two-sided t-test power over that
prior and the user's planned distribution of non-null effects.

## Usage

``` r
power_fdr_prior_n(
  var,
  df,
  delta_mean,
  delta_se = 0,
  fdr = 0.05,
  power = 0.8,
  pi0 = 0.95,
  max_n = 100
)
```

## Arguments

- var:

  real numeric vector of unweighted, unmoderated within-group variance
  estimates on the response scale (usually log2 intensity). Non-finite
  values are excluded; finite negative variances are errors.

- df:

  real numeric vector of residual degrees of freedom, either length one
  (shared by all features) or length(var). Non-finite and non-positive
  values are excluded together with their variance estimate.

- delta_mean:

  one finite mean difference among non-null features, on the same scale
  as the pilot response; negative or zero means are allowed.

- delta_se:

  one finite, non-negative standard deviation of the true differences
  among non-null features. This is not the standard error of an
  estimated mean effect. Zero specifies one fixed planned effect.

- fdr:

  one target false discovery rate, strictly between zero and one.

- power:

  one target average power among non-null features, strictly between
  zero and one.

- pi0:

  one assumed proportion of truly null features, strictly between zero
  and one; it is not estimated from the pilot.

- max_n:

  one whole number of at least two: the largest sample count per group
  to consider.

## Value

A list containing `N` (first integer sample count per group reaching the
target, or `NA_integer_`), `power_achieved` (at `N`, or at `max_n` when
the target is not reached), `power_curve` (a tibble of `n` and `power`
for 2 through `max_n`), `sig_level`, `var_prior`, `df_prior`,
`n_features_used`, and the original inputs `var`, `df`, `delta_mean`,
`delta_se`, `fdr`, `power`, `pi0`, `max_n`.

## Details

At least three usable features, including two positive variances, are
needed to estimate the prior with
`limma::squeezeVar(var, df, robust = FALSE)`. Zero variances are
retained; any limma warnings about unreliable or offset variances remain
visible. If the estimated prior parameters are \\s_0^2\\ and \\d_0\\,
true variances follow an inverse-gamma distribution with shape \\d_0/2\\
and scale \\d_0 s_0^2/2\\. True non-null differences independently
follow `Normal(delta_mean, delta_se)` (mean and standard deviation).

Deterministic midpoint quadrature uses gamma-precision and normal-effect
quantiles: 256 variance nodes (1024 when prior df is below four) and 512
effect nodes. Infinite prior df gives the single variance `var_prior`;
zero `delta_se` gives the single effect `delta_mean`.

The power curve uses the fixed threshold from
[`fdr_sig_level()`](https://fgcz.github.io/prolfqua/reference/fdr_sig_level.md)
at the requested target power. The chosen sample size achieves the FDR
planning bound under those assumptions; points below the target power
need not. Prospective tests use `2 * n - 2` degrees of freedom, not
moderated df. A moderated test gains prior degrees of freedom, so this
ordinary t-test sample size is conservative relative to that test under
the same model. This function estimates the variance distribution, not
effects or `pi0`.

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
[`fdr_sig_level()`](https://fgcz.github.io/prolfqua/reference/fdr_sig_level.md),
[`lfq_power_t_test_proteins()`](https://fgcz.github.io/prolfqua/reference/lfq_power_t_test_proteins.md),
[`lfq_power_t_test_quantiles_V2()`](https://fgcz.github.io/prolfqua/reference/lfq_power_t_test_quantiles_V2.md),
[`plot_stat_density()`](https://fgcz.github.io/prolfqua/reference/plot_stat_density.md),
[`plot_stat_density_median()`](https://fgcz.github.io/prolfqua/reference/plot_stat_density_median.md),
[`plot_stat_violin()`](https://fgcz.github.io/prolfqua/reference/plot_stat_violin.md),
[`plot_stat_violin_median()`](https://fgcz.github.io/prolfqua/reference/plot_stat_violin_median.md),
[`plot_stdv_vs_mean()`](https://fgcz.github.io/prolfqua/reference/plot_stdv_vs_mean.md),
[`summarize_stats()`](https://fgcz.github.io/prolfqua/reference/summarize_stats.md)

## Examples

``` r
variances <- 0.01 * stats::qf((seq_len(200) - 0.5) / 200, df1 = 8, df2 = 10)
planned <- power_fdr_prior_n(variances, df = 8, delta_mean = 0.15, delta_se = 0.02)
planned[c("N", "power_achieved", "var_prior", "df_prior", "n_features_used")]
#> $N
#> [1] 21
#> 
#> $power_achieved
#> [1] 0.8177378
#> 
#> $var_prior
#> [1] 0.01000914
#> 
#> $df_prior
#> [1] 10.07975
#> 
#> $n_features_used
#> [1] 200
#> 
```
