# R6 class representing modelling result

Holds one fitted model per subject (for example protein) from
[`build_model()`](https://fgcz.github.io/prolfqua/reference/build_model.md),
together with the model table (`model_df`: fit, residual degrees of
freedom, sigma, convergence) and accessors for coefficients, ANOVA
tables and their plots. Contrasts are computed from it with
[Contrasts](https://fgcz.github.io/prolfqua/reference/Contrasts.md).

## Value

An R6 class generator.

## See also

Other modelling:
[`AnovaExtractor`](https://fgcz.github.io/prolfqua/reference/AnovaExtractor.md),
[`Contrasts`](https://fgcz.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsComparison`](https://fgcz.github.io/prolfqua/reference/ContrastsComparison.md),
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
[`ContrastsROPECANestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECANestedFacade.md),
[`ContrastsROPECAV2`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAV2.md),
[`ContrastsROPECAV2NestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECAV2NestedFacade.md),
[`ContrastsRfitFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsRfitFacade.md),
[`ContrastsRfitImputeFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsRfitImputeFacade.md),
[`ContrastsTable`](https://fgcz.github.io/prolfqua/reference/ContrastsTable.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
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

## Super class

[`ModelInterface`](https://fgcz.github.io/prolfqua/reference/ModelInterface.md)
-\> `Model`

## Public fields

- `model_df`:

  data.frame with modelling data and model.

- `model_name`:

  name of model

- `subject_id`:

  e.g. protein_Id

- `model_strategy`:

  function to create the models

- `p.adjust`:

  function to adjust p-values

## Methods

### Public methods

- [`Model$new()`](#method-Model-initialize)

- [`Model$get_coefficients()`](#method-Model-get_coefficients)

- [`Model$get_anova()`](#method-Model-get_anova)

- [`Model$clone()`](#method-Model-clone)

Inherited methods

- [`ModelInterface$anova_histogram()`](https://fgcz.github.io/prolfqua/reference/ModelInterface.html#method-anova_histogram)
- [`ModelInterface$coef_histogram()`](https://fgcz.github.io/prolfqua/reference/ModelInterface.html#method-coef_histogram)
- [`ModelInterface$coef_pairs()`](https://fgcz.github.io/prolfqua/reference/ModelInterface.html#method-coef_pairs)
- [`ModelInterface$coef_volcano()`](https://fgcz.github.io/prolfqua/reference/ModelInterface.html#method-coef_volcano)
- [`ModelInterface$get_provenance_log()`](https://fgcz.github.io/prolfqua/reference/ModelInterface.html#method-get_provenance_log)
- [`ModelInterface$provenance()`](https://fgcz.github.io/prolfqua/reference/ModelInterface.html#method-provenance)

------------------------------------------------------------------------

### `Model$new()`

initialize

#### Usage

    Model$new(
      model_df,
      model_strategy,
      model_name,
      subject_id = "protein_Id",
      p.adjust = prolfqua::adjust_p_values
    )

#### Arguments

- `model_df`:

  dataframe with modelling results

- `model_strategy`:

  model_strategy see
  [`strategy_lmer`](https://fgcz.github.io/prolfqua/reference/strategy.md)

- `model_name`:

  name of model

- `subject_id`:

  subject column name

- `p.adjust`:

  method to adjust p-values

------------------------------------------------------------------------

### `Model$get_coefficients()`

return model coefficient table

#### Usage

    Model$get_coefficients()

------------------------------------------------------------------------

### `Model$get_anova()`

return anova table

#### Usage

    Model$get_anova()

------------------------------------------------------------------------

### `Model$clone()`

The objects of this class are cloneable with this method.

#### Usage

    Model$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r


istar <- sim_lfq_data_peptide_config(n_proteins = 20)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 110 rows with NA response for feature and sample combinations absent from the input.
lfqdata <- LFQData$new(istar$data, istar$config)
lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#> Transformed response column added: log2_abundance
model_name <- "f_condtion_r_peptide"
formula_randomPeptide <-
  strategy_lmer(paste0(lfqdata$response(), " ~ group_ + (1 | peptide_Id)"),
   model_name = model_name)
mod <- prolfqua::build_model(
 lfqdata,
 formula_randomPeptide,
 model_name = model_name)
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> Warning: fit (f_condtion_r_peptide): 7 fits failed (e.g. 0EfVhX~5954, 7cbcrd~8305, CGzoYe~2857, Fl4JiV~7526, JV3Z7t~2956, ...): Error: grouping factors must have > 1 sampled level
#> See $provenance().
#> fit (f_condtion_r_peptide): 4 singular fits (e.g. BEJI92~9143, HC8K98~4958, HvIpHG~4015, R2i6w7~0288): coefficients not estimable or fewer than 2 residual df
#> fit (f_condtion_r_peptide): 102 missing responses dropped (e.g. 0m5WN4~1448, 9VUkAq~4562, At886V~3296, BEJI92~9143, CtOJ9t~2837, ...): NA response (missing observation) left out of the fit
#> See $provenance().

mod$model_df
#> # A tibble: 20 × 12
#> # Groups:   protein_Id [20]
#>    protein_Id  data     linear_model fit_warnings has_model_fit isSingular
#>    <chr>       <list>   <list>       <list>       <lgl>         <lgl>     
#>  1 0EfVhX~5954 <tibble> <chr [1]>    <chr [0]>    FALSE         NA        
#>  2 0m5WN4~1448 <tibble> <lmrMdLmT>   <chr [0]>    TRUE          FALSE     
#>  3 7cbcrd~8305 <tibble> <chr [1]>    <chr [0]>    FALSE         NA        
#>  4 9VUkAq~4562 <tibble> <lmrMdLmT>   <chr [0]>    TRUE          FALSE     
#>  5 At886V~3296 <tibble> <lmrMdLmT>   <chr [0]>    TRUE          FALSE     
#>  6 BEJI92~9143 <tibble> <lmrMdLmT>   <chr [0]>    TRUE          TRUE      
#>  7 CGzoYe~2857 <tibble> <chr [1]>    <chr [0]>    FALSE         NA        
#>  8 CtOJ9t~2837 <tibble> <lmrMdLmT>   <chr [0]>    TRUE          FALSE     
#>  9 DoWup2~2934 <tibble> <lmrMdLmT>   <chr [0]>    TRUE          FALSE     
#> 10 DuwH7n~3402 <tibble> <lmrMdLmT>   <chr [0]>    TRUE          FALSE     
#> 11 Fl4JiV~7526 <tibble> <chr [1]>    <chr [0]>    FALSE         NA        
#> 12 HC8K98~4958 <tibble> <lmrMdLmT>   <chr [0]>    TRUE          TRUE      
#> 13 HvIpHG~4015 <tibble> <lmrMdLmT>   <chr [0]>    TRUE          TRUE      
#> 14 I1Jk2Z~0821 <tibble> <lmrMdLmT>   <chr [0]>    TRUE          FALSE     
#> 15 JV3Z7t~2956 <tibble> <chr [1]>    <chr [0]>    FALSE         NA        
#> 16 JcKVfU~0815 <tibble> <chr [1]>    <chr [0]>    FALSE         NA        
#> 17 JfvT8X~2727 <tibble> <lmrMdLmT>   <chr [0]>    TRUE          FALSE     
#> 18 R2i6w7~0288 <tibble> <lmrMdLmT>   <chr [0]>    TRUE          TRUE      
#> 19 SGIVBl~9558 <tibble> <lmrMdLmT>   <chr [0]>    TRUE          FALSE     
#> 20 r2J0Eh~2687 <tibble> <chr [1]>    <chr [0]>    FALSE         NA        
#> # ℹ 6 more variables: df.residual <dbl>, sigma <dbl>, nr_coef <int>,
#> #   nr_coef_not_NA <int>, converged <lgl>, iterations <int>
aovtable  <- mod$get_anova()
mod$get_coefficients()
#> # A tibble: 39 × 9
#> # Groups:   protein_Id [13]
#>    protein_Id  factor     Estimate Std..Error     df t.value Pr...t.. isSingular
#>    <chr>       <chr>         <dbl>      <dbl>  <dbl>   <dbl>    <dbl> <lgl>     
#>  1 0m5WN4~1448 (Intercep…   4.16       0.0440   4.58  94.6   1.01e- 8 FALSE     
#>  2 0m5WN4~1448 group_B     -0.0614     0.0686  16.2   -0.896 3.83e- 1 FALSE     
#>  3 0m5WN4~1448 group_Ctrl   0.0815     0.0622  16.1    1.31  2.09e- 1 FALSE     
#>  4 9VUkAq~4562 (Intercep…   4.22       0.0449  27.8   94.0   2.49e-36 FALSE     
#>  5 9VUkAq~4562 group_B     -0.0902     0.0395 144.    -2.29  2.38e- 2 FALSE     
#>  6 9VUkAq~4562 group_Ctrl  -0.0699     0.0400 145.    -1.75  8.26e- 2 FALSE     
#>  7 At886V~3296 (Intercep…   4.12       0.0434   7.61  94.9   5.63e-13 FALSE     
#>  8 At886V~3296 group_B      0.0771     0.0393  49.1    1.96  5.55e- 2 FALSE     
#>  9 At886V~3296 group_Ctrl   0.122      0.0383  49.1    3.19  2.52e- 3 FALSE     
#> 10 BEJI92~9143 (Intercep…   4.51       0.123   41.0   36.7   5.60e-33 TRUE      
#> # ℹ 29 more rows
#> # ℹ 1 more variable: nr_coef <int>
mod$coef_histogram()
#> $plot

#> 
#> $name
#> [1] "Coef_Histogram_f_condtion_r_peptide.pdf"
#> 
mod$coef_volcano()
#> $plot

#> 
#> $name
#> [1] "Coef_volcano_plot_f_condtion_r_peptide.pdf"
#> 
mod$coef_pairs()
#> $plot
#> # A tibble: 13 × 4
#>    subject_id  `(Intercept)`  group_B group_Ctrl
#>    <chr>               <dbl>    <dbl>      <dbl>
#>  1 0m5WN4~1448          4.16 -0.0614      0.0815
#>  2 9VUkAq~4562          4.22 -0.0902     -0.0699
#>  3 At886V~3296          4.12  0.0771      0.122 
#>  4 BEJI92~9143          4.51  0.204       0.115 
#>  5 CtOJ9t~2837          4.84 -0.0838     -0.387 
#>  6 DoWup2~2934          4.43  0.316      -0.0374
#>  7 DuwH7n~3402          4.11 -0.00356     0.0990
#>  8 HC8K98~4958          3.89  0.322       0.0752
#>  9 HvIpHG~4015          4.11  0.661      -0.161 
#> 10 I1Jk2Z~0821          3.84  0.219       0.276 
#> 11 JfvT8X~2727          4.48 -0.0844     -0.194 
#> 12 R2i6w7~0288          4.65 -0.410      -0.396 
#> 13 SGIVBl~9558          5.12 -0.158      -0.344 
#> 
#> $name
#> [1] "Coef_Pairsplot_f_condtion_r_peptide.pdf"
#> 
mod$anova_histogram()
#> $plot

#> 
#> $name
#> [1] "Anova_p.values_f_condtion_r_peptide.pdf"
#> 
```
