# get linfct from model

get linfct from model

## Usage

``` r
linfct_from_model(m, as_list = TRUE)
```

## Arguments

- m:

  linear model

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
[`my_glht()`](https://wolski.github.io/prolfqua/reference/my_glht.md),
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
m <- sim_make_model_lm("factors")
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
#> Joining with `by = join_by(protein_Id)`
linfct <- linfct_from_model(m, as_list = TRUE)

linfct$linfct_factors
#>             (Intercept) TreatmentB BackgroundZ
#> BackgroundX           1        0.5         0.0
#> BackgroundZ           1        0.5         1.0
#> TreatmentA            1        0.0         0.5
#> TreatmentB            1        1.0         0.5
linfct$linfct_interactions
#>                        (Intercept) TreatmentB BackgroundZ
#> TreatmentA:BackgroundX           1          0           0
#> TreatmentA:BackgroundZ           1          0           1
#> TreatmentB:BackgroundX           1          1           0
#> TreatmentB:BackgroundZ           1          1           1
lf <- matrix(
c(1, 1, 1, 1, 0.5, 0.5, 0, 1, 0, 1, 0.5, 0.5),
nrow = 4,
byrow = FALSE,
dimnames = list(c("BackgroundX", "BackgroundZ", "TreatmentA", "TreatmentB"),
                c("(Intercept)", "TreatmentB", "BackgroundZ"))
)
stopifnot(lf == linfct$linfct_factors)
m <- sim_make_model_lm("interaction")
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
#> Joining with `by = join_by(protein_Id)`
linfct <- linfct_from_model(m)

m <- lm(Petal.Width ~ Species, data = iris)
linfct_from_model(m)
#> $linfct_factors
#>                   (Intercept) Speciesversicolor Speciesvirginica
#> Speciessetosa               1                 0                0
#> Speciesversicolor           1                 1                0
#> Speciesvirginica            1                 0                1
#> 
#> $linfct_interactions
#>                   (Intercept) Speciesversicolor Speciesvirginica
#> Speciessetosa               1                 0                0
#> Speciesversicolor           1                 1                0
#> Speciesvirginica            1                 0                1
#> 
xx <- data.frame( Y = 1:10 , Condition = c(rep("a",5), rep("b",5)) )
m <- lm(Y ~ Condition, data = xx)
linfct_from_model(m)
#> $linfct_factors
#>            (Intercept) Conditionb
#> Conditiona           1          0
#> Conditionb           1          1
#> 
#> $linfct_interactions
#>            (Intercept) Conditionb
#> Conditiona           1          0
#> Conditionb           1          1
#> 
xx <- data.frame( Y = 1:10 , Condition = c(rep("a",5), rep("b.b",5)) )
m <- lm(Y ~ Condition, data = xx)
linfct_from_model(m)
#> $linfct_factors
#>              (Intercept) Conditionb.b
#> Conditiona             1            0
#> Conditionb.b           1            1
#> 
#> $linfct_interactions
#>              (Intercept) Conditionb.b
#> Conditiona             1            0
#> Conditionb.b           1            1
#> 
xx <- data.frame( Y = 1:10 , Condition = c(rep("a",5), rep("ab",5)) )
m <- lm(Y ~ Condition, data = xx)
linfct_from_model(m)
#> $linfct_factors
#>             (Intercept) Conditionab
#> Conditiona            1           0
#> Conditionab           1           1
#> 
#> $linfct_interactions
#>             (Intercept) Conditionab
#> Conditiona            1           0
#> Conditionab           1           1
#> 
```
