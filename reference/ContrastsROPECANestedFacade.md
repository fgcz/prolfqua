# ROPECA contrast analysis facade for nested input

Encapsulates the pipeline:
[`strategy_lm`](https://fgcz.github.io/prolfqua/reference/strategy.md)
-\>
[`build_model`](https://fgcz.github.io/prolfqua/reference/build_model.md)
-\>
[`Contrasts`](https://fgcz.github.io/prolfqua/reference/Contrasts.md)
-\>
[`ContrastsROPECA`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECA.md).

Input must be nested LFQData (peptide subject_id strictly below the
protein hierarchy).

`ropeca_nested` is a derivation of PECA and ROPECA (Suomi and Elo 2017,
Scientific Reports 7, 5869; Bioconductor package PECA,
[doi:10.18129/B9.bioc.PECA](https://doi.org/10.18129/B9.bioc.PECA) ). It
fits one linear model per peptide, computes the contrast, takes the
median of the signed peptide p-values per protein and corrects it with
the beta distribution of the median of `n` uniform p-values, `n` being
the number of peptides with a p-value, at most 10. The protein fold
change is the median peptide fold change. What differs from PECA and
ROPECA:

- no ROTS (reproducibility-optimized test statistic) step: peptide
  p-values come from ordinary linear models;

- any fixed-effect design and contrast, not only two groups;

- `n` is capped at 10, so proteins with more peptides get larger
  p-values than in PECA;

- with an even number of peptides, the median is taken of the signed
  p-values, not of the t-statistics.

The test is conservative under the null hypothesis; `ropecaV2_nested` is
the calibrated variant.

## Value

An R6 class generator.

## Comparison with PECA

On the IonStar data (`data_ionstar`, dilutions a and e, 767 peptides
observed in all 8 samples, 123 proteins), `ropeca_nested` and
`PECA::PECA_df(test = "t", type = "median")` give the same fold changes
(largest difference 1e-14) and the same p-values for the 58 proteins
with an odd number of peptides up to 10 (largest difference 3e-14). For
the 50 proteins with an even number, the median ratio of the p-values is
1.02 (Spearman correlation 0.99); for the 15 proteins with more than 10
peptides, prolfqua's p-values are larger by a median factor of 1.43
(Spearman correlation 0.87), because of the cap. Code:
`data-raw/compare_peca.R` in the source repository.

## See also

Other modelling:
[`AnovaExtractor`](https://fgcz.github.io/prolfqua/reference/AnovaExtractor.md),
[`Contrasts`](https://fgcz.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsDEqMSFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsDEqMSFacade.md),
[`ContrastsDEqMSVoomFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsDEqMSVoomFacade.md),
[`ContrastsFacadeBase`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.md),
[`ContrastsFirth`](https://fgcz.github.io/prolfqua/reference/ContrastsFirth.md),
[`ContrastsFirthFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsFirthFacade.md),
[`ContrastsHurdleFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsHurdleFacade.md),
[`ContrastsHurdleNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsHurdleNestedFacade.md),
[`ContrastsLMFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLMFacade.md),
[`ContrastsLMImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLMImputeFacade.md),
[`ContrastsLimma`](https://fgcz.github.io/prolfqua/reference/ContrastsLimma.md),
[`ContrastsLimmaFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimmaFacade.md),
[`ContrastsLimmaImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimmaImputeFacade.md),
[`ContrastsLimmaVoomFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimmaVoomFacade.md),
[`ContrastsLimmaVoomImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimmaVoomImputeFacade.md),
[`ContrastsLimpaFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimpaFacade.md),
[`ContrastsLimpaNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLimpaNestedFacade.md),
[`ContrastsLmerNestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsLmerNestedFacade.md),
[`ContrastsModerated`](https://fgcz.github.io/prolfqua/reference/ContrastsModerated.md),
[`ContrastsModeratedDEqMS`](https://fgcz.github.io/prolfqua/reference/ContrastsModeratedDEqMS.md),
[`ContrastsPlotter`](https://fgcz.github.io/prolfqua/reference/ContrastsPlotter.md),
[`ContrastsRLMFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsRLMFacade.md),
[`ContrastsROPECA`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECA.md),
[`ContrastsROPECAV2`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAV2.md),
[`ContrastsROPECAV2NestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAV2NestedFacade.md),
[`ContrastsRfitFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsRfitFacade.md),
[`ContrastsRfitImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsRfitImputeFacade.md),
[`ContrastsTable`](https://fgcz.github.io/prolfqua/reference/ContrastsTable.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`Model`](https://fgcz.github.io/prolfqua/reference/Model.md),
[`ModelFirth`](https://fgcz.github.io/prolfqua/reference/ModelFirth.md),
[`ModelLimma`](https://fgcz.github.io/prolfqua/reference/ModelLimma.md),
[`StrategyLM`](https://fgcz.github.io/prolfqua/reference/StrategyLM.md),
[`StrategyLimma`](https://fgcz.github.io/prolfqua/reference/StrategyLimma.md),
[`StrategyLimpa`](https://fgcz.github.io/prolfqua/reference/StrategyLimpa.md),
[`StrategyLmer`](https://fgcz.github.io/prolfqua/reference/StrategyLmer.md),
[`StrategyLogistf`](https://fgcz.github.io/prolfqua/reference/StrategyLogistf.md),
[`StrategyRLM`](https://fgcz.github.io/prolfqua/reference/StrategyRLM.md),
[`StrategyRfit`](https://fgcz.github.io/prolfqua/reference/StrategyRfit.md),
[`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md),
[`build_model()`](https://fgcz.github.io/prolfqua/reference/build_model.md),
[`build_model_glm_protein()`](https://fgcz.github.io/prolfqua/reference/build_model_glm_protein.md),
[`build_model_impute()`](https://fgcz.github.io/prolfqua/reference/build_model_impute.md),
[`build_model_limma()`](https://fgcz.github.io/prolfqua/reference/build_model_limma.md),
[`build_model_limma_impute()`](https://fgcz.github.io/prolfqua/reference/build_model_limma_impute.md),
[`build_model_limma_voom()`](https://fgcz.github.io/prolfqua/reference/build_model_limma_voom.md),
[`build_model_limma_voom_impute()`](https://fgcz.github.io/prolfqua/reference/build_model_limma_voom_impute.md),
[`build_model_limpa()`](https://fgcz.github.io/prolfqua/reference/build_model_limpa.md),
[`compute_borrowed_variance()`](https://fgcz.github.io/prolfqua/reference/compute_borrowed_variance.md),
[`compute_borrowed_variance_limma()`](https://fgcz.github.io/prolfqua/reference/compute_borrowed_variance_limma.md),
[`df.residual.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/df.residual.rfit_prolfqua.md),
[`group_label()`](https://fgcz.github.io/prolfqua/reference/group_label.md),
[`impute_from_model()`](https://fgcz.github.io/prolfqua/reference/impute_from_model.md),
[`impute_refit_singular()`](https://fgcz.github.io/prolfqua/reference/impute_refit_singular.md),
[`list_facades()`](https://fgcz.github.io/prolfqua/reference/list_facades.md),
[`lookup_facade()`](https://fgcz.github.io/prolfqua/reference/lookup_facade.md),
[`lr_test()`](https://fgcz.github.io/prolfqua/reference/lr_test.md),
[`merge_contrasts_results()`](https://fgcz.github.io/prolfqua/reference/merge_contrasts_results.md),
[`new_imputed_model()`](https://fgcz.github.io/prolfqua/reference/new_imputed_model.md),
[`nobs.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/nobs.rfit_prolfqua.md),
[`pivot_model_contrasts_to_wide()`](https://fgcz.github.io/prolfqua/reference/pivot_model_contrasts_to_wide.md),
[`prolfqua_analysis_methods`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md),
[`register_facade()`](https://fgcz.github.io/prolfqua/reference/register_facade.md),
[`sigma.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/sigma.rfit_prolfqua.md),
[`strategy_limma()`](https://fgcz.github.io/prolfqua/reference/strategy_limma.md),
[`strategy_limpa()`](https://fgcz.github.io/prolfqua/reference/strategy_limpa.md),
[`strategy_logistf()`](https://fgcz.github.io/prolfqua/reference/strategy.md),
[`vcov.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/vcov.rfit_prolfqua.md)

## Super classes

[`ContrastsInterface`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.md)
-\>
[`ContrastsFacadeBase`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.md)
-\> `ContrastsROPECANestedFacade`

## Methods

### Public methods

- [`ContrastsROPECANestedFacade$new()`](#method-ContrastsROPECANestedFacade-initialize)

- [`ContrastsROPECANestedFacade$get_contrasts()`](#method-ContrastsROPECANestedFacade-get_contrasts)

- [`ContrastsROPECANestedFacade$get_Plotter()`](#method-ContrastsROPECANestedFacade-get_Plotter)

- [`ContrastsROPECANestedFacade$to_wide()`](#method-ContrastsROPECANestedFacade-to_wide)

- [`ContrastsROPECANestedFacade$clone()`](#method-ContrastsROPECANestedFacade-clone)

Inherited methods

- [`ContrastsInterface$column_description()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-column_description)
- [`ContrastsInterface$contrast_summary_table()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-contrast_summary_table)
- [`ContrastsInterface$extra_artifacts()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-extra_artifacts)
- [`ContrastsInterface$filter_significant()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-filter_significant)
- [`ContrastsInterface$get_contrast_sides()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_contrast_sides)
- [`ContrastsInterface$get_ora()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_ora)
- [`ContrastsInterface$get_rank()`](https://fgcz.github.io/prolfqua/reference/ContrastsInterface.html#method-get_rank)
- [`ContrastsFacadeBase$.recorded_provenance()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-.recorded_provenance)
- [`ContrastsFacadeBase$get_config()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-get_config)
- [`ContrastsFacadeBase$get_missing()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-get_missing)
- [`ContrastsFacadeBase$provenance()`](https://fgcz.github.io/prolfqua/reference/ContrastsFacadeBase.html#method-provenance)

------------------------------------------------------------------------

### `ContrastsROPECANestedFacade$new()`

initialize

#### Usage

    ContrastsROPECANestedFacade$new(lfqdata, modelstr, contrasts)

#### Arguments

- `lfqdata`:

  nested LFQData (peptide-level)

- `modelstr`:

  model formula string (e.g. "~ group\_")

- `contrasts`:

  named character vector of contrasts

------------------------------------------------------------------------

### `ContrastsROPECANestedFacade$get_contrasts()`

Get contrast results with standardized column names.

ROPECA's `beta.based.significance` is mapped to `p.value` and
`FDR.beta.based.significance` to `FDR`.

Columns not directly produced by ROPECA are derived heuristically:

- `std.error = diff / statistic` (algebraic: t = estimate / SE)

- `sigma = mad.estimate` (MAD of peptide fold changes)

- `df = n_not_na` (number of contributing peptides)

- `conf.low/conf.high` via `diff +/- qt(0.975, df) * |std.error|`

#### Usage

    ContrastsROPECANestedFacade$get_contrasts(all = FALSE)

#### Arguments

- `all`:

  `FALSE` (default) returns the standard columns; `TRUE` appends the
  remaining ROPECA columns

------------------------------------------------------------------------

### `ContrastsROPECANestedFacade$get_Plotter()`

get ContrastsPlotter (uses standardized column names)

#### Usage

    ContrastsROPECANestedFacade$get_Plotter(fc_threshold = 2, fdr_threshold = 0.1)

#### Arguments

- `fc_threshold`:

  fold change threshold

- `fdr_threshold`:

  FDR threshold

------------------------------------------------------------------------

### `ContrastsROPECANestedFacade$to_wide()`

convert results to wide format

#### Usage

    ContrastsROPECANestedFacade$to_wide(columns = c("p.value", "FDR", "statistic"))

#### Arguments

- `columns`:

  value columns to pivot, default `c("p.value", "FDR", "statistic")`

------------------------------------------------------------------------

### `ContrastsROPECANestedFacade$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsROPECANestedFacade$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
lfqdata <- LFQData$new(istar$data, istar$config)
lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#> Transformed response column added: log2_abundance
contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
fa <- ContrastsROPECANestedFacade$new(lfqdata, "~ group_", contrasts)
#> fit (Model): 36 missing responses dropped (e.g. 0EfVhX~0087~ITLb4x1q, 0EfVhX~0087~ahQLlQY7, 0EfVhX~0087~dJkdz7so, 7cbcrd~5725~D5dQ4nKk, 9VUkAq~4703~eIC06D7g, ...): NA response (missing observation) left out of the fit
#> See $provenance().
head(fa$get_contrasts())
#> # A tibble: 6 × 14
#>   modelName   estimate_type protein_Id contrast avgAbd    diff     FDR statistic
#>   <chr>       <chr>         <chr>      <chr>     <dbl>   <dbl>   <dbl>     <dbl>
#> 1 ropeca_nes… observed      0EfVhX~00… A_vs_Ct…   4.27 -0.0742 5.28e-2     -1.75
#> 2 ropeca_nes… observed      7cbcrd~57… A_vs_Ct…   4.51  0.741  9.91e-5      8.79
#> 3 ropeca_nes… observed      9VUkAq~47… A_vs_Ct…   4.47 -0.598  6.91e-6    -12.7 
#> 4 ropeca_nes… observed      BEJI92~52… A_vs_Ct…   4.23  0.277  1.87e-3      3.94
#> 5 ropeca_nes… observed      CGzoYe~21… A_vs_Ct…   4.76 -0.310  3.74e-5     -9.26
#> 6 ropeca_nes… observed      DoWup2~58… A_vs_Ct…   4.43  0.295  1.38e-6     14.7 
#> # ℹ 6 more variables: std.error <dbl>, df <int>, p.value <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, sigma <dbl>
fa$to_wide()
#> # A tibble: 10 × 5
#>    protein_Id diff.A_vs_Ctrl p.value.A_vs_Ctrl FDR.A_vs_Ctrl statistic.A_vs_Ctrl
#>    <chr>               <dbl>             <dbl>         <dbl>               <dbl>
#>  1 0EfVhX~00…        -0.0742       0.0423         0.0528                   -1.75
#>  2 7cbcrd~57…         0.741        0.0000496      0.0000991                 8.79
#>  3 9VUkAq~47…        -0.598        0.00000138     0.00000691              -12.7 
#>  4 BEJI92~52…         0.277        0.00112        0.00187                   3.94
#>  5 CGzoYe~21…        -0.310        0.0000150      0.0000374                -9.26
#>  6 DoWup2~58…         0.295        0.000000138    0.00000138               14.7 
#>  7 Fl4JiV~86…         0.0700       0.740          0.823                     2.08
#>  8 HvIpHG~90…        -0.384        0.00000389     0.0000130                -7.28
#>  9 JcKVfU~96…        -0.0634       0.00218        0.00312                  -1.87
#> 10 SGIVBl~57…        -0.122        0.993          0.993                    -3.88
```
