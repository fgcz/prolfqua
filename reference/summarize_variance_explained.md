# Variance explained by the design, per feature

Fits a least-squares model to every feature of a log-transformed
[LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) and
splits each feature's variance with its ANOVA table. Use it to check a
normalization: a good normalization lowers the within-group variance
more than the overall variance, so the share explained by the design,
the model R^2, rises.

## Usage

``` r
summarize_variance_explained(lfqdata, modelstr = NULL)
```

## Arguments

- lfqdata:

  a log-transformed
  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md)
  object.

- modelstr:

  model formula string, as in
  [`strategy_lm()`](https://fgcz.github.io/prolfqua/reference/strategy.md).
  `NULL` (default) uses the full factorial of
  `lfqdata$relevant_factor_keys()`, for example
  `"log2_abundance ~ Treatment * Background"`.

## Value

A tibble with one row per feature and ANOVA term, residuals included:
the feature keys, `factor`, `Df`, `Sum.Sq`, `SS_total` (the feature's
total sum of squares), `eta2` (`Sum.Sq / SS_total`) and `R2` (the
feature's model R^2, `1 - SS_residual / SS_total`, repeated on each of
its rows).

## Details

One model is fitted per lowest-level feature (all hierarchy keys), the
grouping of
[`summarize_stats()`](https://fgcz.github.io/prolfqua/reference/summarize_stats.md).
Missing responses are omitted within each fit. Only full-rank fits with
more than one residual degree of freedom are kept, the set returned by
[`get_complete_model_fit()`](https://fgcz.github.io/prolfqua/reference/get_complete_model_fit.md),
so all features share the same terms; a message reports how many
features were left out.

Read the model R^2 first. The ANOVA uses sequential (type I) sums of
squares: when a feature is unbalanced, for example through missing
values, the shares of the individual factors depend on the order of the
terms in `modelstr`. R^2 does not.

Without any effect R^2 is positive on average: E\[R^2\] = p / (n - 1)
for p model degrees of freedom and n observed values. Compare
normalizations of the same samples with the same missing values;
imputation changes n.

Use R^2 as a check of a normalization, not as a criterion to pick the
normalization with the highest value: that choice favours the group
differences tested afterwards.

## See also

[`pooled_variance_explained()`](https://fgcz.github.io/prolfqua/reference/pooled_variance_explained.md)

Other stats:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`fdr_sig_level()`](https://fgcz.github.io/prolfqua/reference/fdr_sig_level.md),
[`null_split_check()`](https://fgcz.github.io/prolfqua/reference/null_split_check.md),
[`pooled_variance_explained()`](https://fgcz.github.io/prolfqua/reference/pooled_variance_explained.md),
[`power_fdr_prior_n()`](https://fgcz.github.io/prolfqua/reference/power_fdr_prior_n.md)

## Examples

``` r
bb <- prolfqua::sim_lfq_data_protein_config(n_proteins = 20)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 5 unconfigured columns: average_prot_abundance, sd, mean, N, Replicate.
#> setup_analysis: added 34 rows with NA response for feature and sample combinations absent from the input.
lfqdata <- LFQData$new(bb$data, bb$config)
# the simulated intensities are already on the log2 scale
lfqdata$is_transformed(TRUE)
ve <- summarize_variance_explained(lfqdata)
#> fit (Model): 34 missing responses dropped (e.g. 0m5WN4~1448, 9VUkAq~4562, At886V~3296, BEJI92~9143, CtOJ9t~2837, ...): NA response (missing observation) left out of the fit
#> See $provenance().
head(ve)
#> # A tibble: 6 × 7
#>   protein_Id  factor       Df Sum.Sq SS_total  eta2    R2
#>   <chr>       <chr>     <int>  <dbl>    <dbl> <dbl> <dbl>
#> 1 0EfVhX~5954 group_        2 53.0      76.3  0.694 0.694
#> 2 0EfVhX~5954 Residuals     9 23.3      76.3  0.306 0.694
#> 3 0m5WN4~1448 group_        2  0.892     4.48 0.199 0.199
#> 4 0m5WN4~1448 Residuals     5  3.59      4.48 0.801 0.199
#> 5 7cbcrd~8305 group_        2 19.1      25.0  0.765 0.765
#> 6 7cbcrd~8305 Residuals     9  5.89     25.0  0.235 0.765
pooled_variance_explained(ve)
#>      factor     eta2
#> 1    group_ 0.681619
#> 2 Residuals 0.318381
```
