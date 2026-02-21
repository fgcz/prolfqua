# apply multcomp::glht method to linfct

apply multcomp::glht method to linfct

## Usage

``` r
my_glht(model, linfct, sep = TRUE)
```

## See also

Other modelling:
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsFirth`](https://wolski.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsMissing`](https://wolski.github.io/prolfqua/reference/ContrastsMissing.md),
[`ContrastsModerated`](https://wolski.github.io/prolfqua/reference/ContrastsModerated.md),
[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.md),
[`ContrastsProDA`](https://wolski.github.io/prolfqua/reference/ContrastsProDA.md),
[`ContrastsROPECA`](https://wolski.github.io/prolfqua/reference/ContrastsROPECA.md),
[`ContrastsTable`](https://wolski.github.io/prolfqua/reference/ContrastsTable.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`LR_test()`](https://wolski.github.io/prolfqua/reference/LR_test.md),
[`Model`](https://wolski.github.io/prolfqua/reference/Model.md),
[`ModelFirth`](https://wolski.github.io/prolfqua/reference/ModelFirth.md),
[`build_model()`](https://wolski.github.io/prolfqua/reference/build_model.md),
[`build_model_logistf()`](https://wolski.github.io/prolfqua/reference/build_model_logistf.md),
[`contrasts_fisher_exact()`](https://wolski.github.io/prolfqua/reference/contrasts_fisher_exact.md),
[`get_anova_df()`](https://wolski.github.io/prolfqua/reference/get_anova_df.md),
[`get_complete_model_fit()`](https://wolski.github.io/prolfqua/reference/get_complete_model_fit.md),
[`get_p_values_pbeta()`](https://wolski.github.io/prolfqua/reference/get_p_values_pbeta.md),
[`group_label()`](https://wolski.github.io/prolfqua/reference/group_label.md),
[`isSingular_lm()`](https://wolski.github.io/prolfqua/reference/isSingular_lm.md),
[`linfct_all_possible_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_all_possible_contrasts.md),
[`linfct_factors_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_factors_contrasts.md),
[`linfct_from_model()`](https://wolski.github.io/prolfqua/reference/linfct_from_model.md),
[`linfct_matrix_contrasts()`](https://wolski.github.io/prolfqua/reference/linfct_matrix_contrasts.md),
[`merge_contrasts_results()`](https://wolski.github.io/prolfqua/reference/merge_contrasts_results.md),
[`model_analyse()`](https://wolski.github.io/prolfqua/reference/model_analyse.md),
[`model_summary()`](https://wolski.github.io/prolfqua/reference/model_summary.md),
[`moderated_p_limma()`](https://wolski.github.io/prolfqua/reference/moderated_p_limma.md),
[`moderated_p_limma_long()`](https://wolski.github.io/prolfqua/reference/moderated_p_limma_long.md),
[`my_contest()`](https://wolski.github.io/prolfqua/reference/my_contest.md),
[`my_contrast()`](https://wolski.github.io/prolfqua/reference/my_contrast.md),
[`my_contrast_V1()`](https://wolski.github.io/prolfqua/reference/my_contrast_V1.md),
[`my_contrast_V2()`](https://wolski.github.io/prolfqua/reference/my_contrast_V2.md),
[`pivot_model_contrasts_2_Wide()`](https://wolski.github.io/prolfqua/reference/pivot_model_contrasts_2_Wide.md),
[`plot_lmer_peptide_predictions()`](https://wolski.github.io/prolfqua/reference/plot_lmer_peptide_predictions.md),
[`sim_build_models_lm()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lm.md),
[`sim_build_models_lmer()`](https://wolski.github.io/prolfqua/reference/sim_build_models_lmer.md),
[`sim_build_models_logistf()`](https://wolski.github.io/prolfqua/reference/sim_build_models_logistf.md),
[`sim_make_model_lm()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lm.md),
[`sim_make_model_lmer()`](https://wolski.github.io/prolfqua/reference/sim_make_model_lmer.md),
[`strategy_logistf()`](https://wolski.github.io/prolfqua/reference/strategy.md),
[`summary_ROPECA_median_p.scaled()`](https://wolski.github.io/prolfqua/reference/summary_ROPECA_median_p.scaled.md)

## Examples

``` r
mb <- sim_make_model_lm( "interaction")
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
#> Joining with `by = join_by(protein_Id)`
linfct <- linfct_from_model(mb)
names(linfct)
#> [1] "linfct_factors"      "linfct_interactions"
my_glht(mb, linfct$linfct_factors)
#> # A tibble: 4 × 10
#>   contrast    null.value estimate std.error statistic adj.p.value conf.low
#>   <chr>            <dbl>    <dbl>     <dbl>     <dbl>       <dbl>    <dbl>
#> 1 BackgroundX          0     18.6     0.322      57.8    4.44e-16     17.9
#> 2 BackgroundZ          0     18.2     0.322      56.4    6.66e-16     17.5
#> 3 TreatmentA           0     18.8     0.322      58.2    4.44e-16     18.1
#> 4 TreatmentB           0     18.1     0.322      56.1    6.66e-16     17.4
#> # ℹ 3 more variables: conf.high <dbl>, df <int>, sigma <dbl>

m <-  sim_make_model_lm( "factors")
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
#> Joining with `by = join_by(protein_Id)`
linfct <- linfct_from_model(m)$linfct_factors
my_glht(m, linfct)
#> # A tibble: 4 × 10
#>   contrast    null.value estimate std.error statistic adj.p.value conf.low
#>   <chr>            <dbl>    <dbl>     <dbl>     <dbl>       <dbl>    <dbl>
#> 1 BackgroundX          0     18.6     0.551      33.8    4.62e-14     17.4
#> 2 BackgroundZ          0     18.2     0.551      33.0    6.29e-14     17.0
#> 3 TreatmentA           0     18.8     0.551      34.1    4.26e-14     17.6
#> 4 TreatmentB           0     18.1     0.551      32.8    6.83e-14     16.9
#> # ℹ 3 more variables: conf.high <dbl>, df <int>, sigma <dbl>
```
