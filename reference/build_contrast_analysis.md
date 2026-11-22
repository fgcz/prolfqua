# Build a contrast analysis using one of several statistical methods

The factory of prolfqua's modelling: `method` selects a facade class
from the facade registry, and the facade runs the full pipeline from
strategy construction through modelling to contrast computation. Every
facade returns its results through the same interface
(`get_contrasts()`, `get_Plotter()`, `get_missing()`, `to_wide()`), so
the results of different methods have the same columns and can be
compared or combined.

## Usage

``` r
build_contrast_analysis(lfqdata, modelstr, contrasts, method = "lm", ...)
```

## Arguments

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

- modelstr:

  one right-hand-side formula string without the response variable (e.g.
  `"~ group_"`). The response is taken automatically from
  `lfqdata$get_config()$get_response()`.

- contrasts:

  nonempty named character vector of contrast expressions (e.g.
  `c("A_vs_B" = "group_A - group_B")`)

- method:

  a registered facade key. The built-in keys are: linear models `"lm"`,
  `"rlm"`, `"rfit"`; limma `"limma"`, `"limma_voom"`, `"deqms"`,
  `"deqms_voom"`, `"limpa"`, `"limpa_nested"`; limit-of-detection
  imputation `"lm_impute"`, `"limma_impute"`, `"limma_voom_impute"`,
  `"rfit_impute"`; peptide-level models `"lmer_nested"`,
  `"ropeca_nested"`, `"ropecaV2_nested"`; detection models `"firth"`,
  `"binomial_nested"`, `"binomial_firth_nested"`; the hurdle model
  `"hurdle_nested"`, combining intensity and detection. Downstream
  packages may add more via
  [`register_facade`](https://fgcz.github.io/prolfqua/reference/register_facade.md).
  The Keys and column names are case sensitive and must match exactly.
  The authoritative list is
  `names(`[`list_facades`](https://fgcz.github.io/prolfqua/reference/list_facades.md)`())`.
  Defaults to `"lm"`.

- ...:

  additional arguments forwarded to the underlying strategy function
  (e.g. `trend`, `robust` for `strategy_limma`)

## Value

a facade object of the class registered for `method`, see
[`lookup_facade`](https://fgcz.github.io/prolfqua/reference/lookup_facade.md).

## Details

Methods ending in `_nested` take peptide- or precursor-level input and
report protein-level contrasts; the others take input at the level they
report, usually aggregated proteins. `"limpa"` needs input aggregated
with
[`AggregateLimpa`](https://fgcz.github.io/prolfqua/reference/AggregateLimpa.md).

The formula is translated into a design matrix with R's current contrast
options (normally treatment coding). Configured predictors are
categorical labels; their reference level is normally the first sorted
level. Each contrast expression combines group-level coefficient
weights, rather than referring directly to design-matrix coefficient
names. For example, `"group_A - group_Ctrl"` subtracts the control group
mean from the A group mean, including the reference coefficient as
needed. Interaction cell means are named with a colon and quoted with
backticks inside expressions. See the *Under the hood* article for
worked one-factor and `~ A * B` design matrices and linear functions.

Internally, `limpa_nested` and `hurdle_nested` aggregate child features
to proteins. Detection methods (`binomial_nested`,
`binomial_firth_nested`, `firth`, and the detection component of
`hurdle_nested`) encode intensities as detected/not detected. These
aggregation and detection steps have no inverse. The analysis keeps its
input in `analysis$.lfqdata`; for aggregated input this is already at
the modelled hierarchy.
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) expose the
modelled observed responses (inspect each hurdle component separately).
Intensity effects are differences on the modelled scale, hence log2 fold
changes for log2 responses; detection effects are log odds ratios, with
fitted responses reported as probabilities. Transformations performed
before this call retain their input columns in LFQData; this call does
not record a transformation history or infer inverse transformations.

With a fixed design and fixed observations per protein, the per-protein
`lm` and `lmer_nested` fitting loops grow approximately linearly with
the number of proteins and run sequentially. `limma` fits the expression
matrix in one call and has lower per-protein overhead. Single-run
elapsed seconds for 500, 1000 and 2000 proteins were respectively 0.173,
0.288 and 0.592 for `lm`; 0.077, 0.043 and 0.050 for `limma`; and 7.902,
15.314 and 29.667 for `lmer_nested` (R 4.6.2, arm64 macOS, 12 samples,
three peptides per protein for the nested model, no missing responses).
These timings cover facade construction, including fitting, but exclude
lazy contrast extraction. They are illustrative timings, not guarantees:
sample count, design width, missingness and mixed-model convergence
affect cost. Some mixed fits in this run emitted convergence warnings.
Reproduction code and environment details are in
`srr/benchmark-runtime.R` in the source repository.

## Prior art

- Improvements on existing R implementations, applied per protein behind
  one interface: `lm` (stats), `rlm` (MASS), `rfit` (Rfit), `limma` and
  `limma_voom` (limma), `deqms` (DEqMS), `limpa` and `limpa_nested`
  (limpa), `lmer_nested` (lme4, lmerTest), `firth` (logistf).

- Reimplementation in R of a published method: `ropeca_nested`, ROPECA
  (Suomi and Elo 2017, Bioconductor package PECA), on peptide-level
  models of any design.

- First implementations of new methods: the `_impute` methods
  (`lm_impute`, `rfit_impute`, `limma_impute`, `limma_voom_impute`),
  which refit groups without observations at the limit of detection with
  a borrowed variance (Wolski et al. 2023); `ropecaV2_nested`, a median
  peptide z-score with a correlation-aware null; `binomial_nested` and
  `binomial_firth_nested`, per-protein models of peptide detection.

- `hurdle_nested` combines the intensity and detection tests as
  msqrob2's `hypothesisTestHurdle` does, on prolfqua's own components.

The Comparing methods article on the package website compares these
implementations on one dataset. Comparisons with MSstats, msqrob2 and
proDA are in the companion package prolfquabenchmark.

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
[`ContrastsROPECANestedFacade`](https://fgcz.github.io/prolfqua/reference/ContrastsROPECANestedFacade.md),
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
[`build_model()`](https://fgcz.github.io/prolfqua/reference/build_model.md),
[`build_model_glm_protein()`](https://fgcz.github.io/prolfqua/reference/build_model_glm_protein.md),
[`build_model_impute()`](https://fgcz.github.io/prolfqua/reference/build_model_impute.md),
[`build_model_limma()`](https://fgcz.github.io/prolfqua/reference/build_model_limma.md),
[`build_model_limma_impute()`](https://fgcz.github.io/prolfqua/reference/build_model_limma_impute.md),
[`build_model_limma_voom()`](https://fgcz.github.io/prolfqua/reference/build_model_limma_voom.md),
[`build_model_limma_voom_impute()`](https://fgcz.github.io/prolfqua/reference/build_model_limma_voom_impute.md),
[`build_model_limpa()`](https://fgcz.github.io/prolfqua/reference/build_model_limpa.md),
[`build_model_logistf()`](https://fgcz.github.io/prolfqua/reference/build_model_logistf.md),
[`compute_borrowed_variance()`](https://fgcz.github.io/prolfqua/reference/compute_borrowed_variance.md),
[`compute_borrowed_variance_limma()`](https://fgcz.github.io/prolfqua/reference/compute_borrowed_variance_limma.md),
[`compute_contrast()`](https://fgcz.github.io/prolfqua/reference/compute_contrast.md),
[`compute_lmer_contrast()`](https://fgcz.github.io/prolfqua/reference/compute_lmer_contrast.md),
[`df.residual.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/df.residual.rfit_prolfqua.md),
[`get_anova_df()`](https://fgcz.github.io/prolfqua/reference/get_anova_df.md),
[`get_complete_model_fit()`](https://fgcz.github.io/prolfqua/reference/get_complete_model_fit.md),
[`get_p_values_pbeta()`](https://fgcz.github.io/prolfqua/reference/get_p_values_pbeta.md),
[`group_label()`](https://fgcz.github.io/prolfqua/reference/group_label.md),
[`impute_from_model()`](https://fgcz.github.io/prolfqua/reference/impute_from_model.md),
[`impute_refit_singular()`](https://fgcz.github.io/prolfqua/reference/impute_refit_singular.md),
[`linfct_from_model()`](https://fgcz.github.io/prolfqua/reference/linfct_from_model.md),
[`linfct_matrix_contrasts()`](https://fgcz.github.io/prolfqua/reference/linfct_matrix_contrasts.md),
[`list_facades()`](https://fgcz.github.io/prolfqua/reference/list_facades.md),
[`lookup_facade()`](https://fgcz.github.io/prolfqua/reference/lookup_facade.md),
[`lr_test()`](https://fgcz.github.io/prolfqua/reference/lr_test.md),
[`merge_contrasts_results()`](https://fgcz.github.io/prolfqua/reference/merge_contrasts_results.md),
[`model_analyse()`](https://fgcz.github.io/prolfqua/reference/model_analyse.md),
[`model_summary()`](https://fgcz.github.io/prolfqua/reference/model_summary.md),
[`moderated_p_deqms()`](https://fgcz.github.io/prolfqua/reference/moderated_p_deqms.md),
[`moderated_p_deqms_long()`](https://fgcz.github.io/prolfqua/reference/moderated_p_deqms_long.md),
[`moderated_p_limma()`](https://fgcz.github.io/prolfqua/reference/moderated_p_limma.md),
[`new_imputed_model()`](https://fgcz.github.io/prolfqua/reference/new_imputed_model.md),
[`nobs.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/nobs.rfit_prolfqua.md),
[`pivot_model_contrasts_to_wide()`](https://fgcz.github.io/prolfqua/reference/pivot_model_contrasts_to_wide.md),
[`prolfqua_analysis_methods`](https://fgcz.github.io/prolfqua/reference/prolfqua_analysis_methods.md),
[`register_facade()`](https://fgcz.github.io/prolfqua/reference/register_facade.md),
[`sigma.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/sigma.rfit_prolfqua.md),
[`sim_build_models_lm()`](https://fgcz.github.io/prolfqua/reference/sim_build_models_lm.md),
[`sim_build_models_lmer()`](https://fgcz.github.io/prolfqua/reference/sim_build_models_lmer.md),
[`sim_build_models_logistf()`](https://fgcz.github.io/prolfqua/reference/sim_build_models_logistf.md),
[`sim_make_model_lm()`](https://fgcz.github.io/prolfqua/reference/sim_make_model_lm.md),
[`sim_make_model_lmer()`](https://fgcz.github.io/prolfqua/reference/sim_make_model_lmer.md),
[`strategy_limma()`](https://fgcz.github.io/prolfqua/reference/strategy_limma.md),
[`strategy_limpa()`](https://fgcz.github.io/prolfqua/reference/strategy_limpa.md),
[`strategy_logistf()`](https://fgcz.github.io/prolfqua/reference/strategy.md),
[`summary_ropeca_median_p_scaled()`](https://fgcz.github.io/prolfqua/reference/summary_ropeca_median_p_scaled.md),
[`vcov.rfit_prolfqua()`](https://fgcz.github.io/prolfqua/reference/vcov.rfit_prolfqua.md)

## Examples

``` r
istar <- sim_lfq_data_protein_config(n_proteins = 20)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 5 unconfigured columns: average_prot_abundance, sd, mean, N, Replicate.
#> setup_analysis: added 34 rows with NA response for feature and sample combinations absent from the input.
lfqdata <- LFQData$new(istar$data, istar$config)
lfqdata$rename_response("transformedIntensity")
contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

fa_lm    <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "lm")
#> fit (Model): 34 missing responses dropped (e.g. 0m5WN4~1448, 9VUkAq~4562, At886V~3296, BEJI92~9143, CtOJ9t~2837, ...): NA response (missing observation) left out of the fit
#> See $provenance().
head(fa_lm$get_contrasts())
#> # A tibble: 6 × 16
#>   modelName estimate_type protein_Id  contrast     diff avgAbd
#>   <chr>     <chr>         <chr>       <chr>       <dbl>  <dbl>
#> 1 lm        observed      0EfVhX~5954 A_vs_Ctrl  2.72     23.2
#> 2 lm        observed      0m5WN4~1448 A_vs_Ctrl  0.600    17.4
#> 3 lm        observed      7cbcrd~8305 A_vs_Ctrl  2.59     27.0
#> 4 lm        observed      9VUkAq~4562 A_vs_Ctrl  0.0679   19.4
#> 5 lm        observed      At886V~3296 A_vs_Ctrl -1.01     19.1
#> 6 lm        observed      BEJI92~9143 A_vs_Ctrl -0.873    20.9
#> # ℹ 10 more variables: std.error.unmoderated <dbl>, df.unmoderated <int>,
#> #   std.error <dbl>, statistic <dbl>, df <dbl>, p.value <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, sigma <dbl>, FDR <dbl>

fa_limma <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "limma")
#> fit (limma): 34 missing responses dropped (e.g. 0m5WN4~1448, 9VUkAq~4562, At886V~3296, BEJI92~9143, CtOJ9t~2837, ...): NA response (missing observation) left out of the fit
#> See $provenance().
head(fa_limma$get_contrasts())
#> # A tibble: 6 × 16
#>   modelName estimate_type protein_Id  contrast     diff    FDR
#>   <chr>     <chr>         <chr>       <chr>       <dbl>  <dbl>
#> 1 limma     observed      0EfVhX~5954 A_vs_Ctrl  2.72   0.188 
#> 2 limma     observed      0m5WN4~1448 A_vs_Ctrl  0.600  0.623 
#> 3 limma     observed      7cbcrd~8305 A_vs_Ctrl  2.59   0.0271
#> 4 limma     observed      9VUkAq~4562 A_vs_Ctrl  0.0679 0.967 
#> 5 limma     observed      At886V~3296 A_vs_Ctrl -1.01   0.623 
#> 6 limma     observed      BEJI92~9143 A_vs_Ctrl -0.873  0.623 
#> # ℹ 10 more variables: std.error.unmoderated <dbl>, df.unmoderated <dbl>,
#> #   std.error <dbl>, statistic <dbl>, p.value <dbl>, sigma <dbl>, df <dbl>,
#> #   conf.low <dbl>, conf.high <dbl>, avgAbd <dbl>

fa_deqms <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "deqms")
#> fit (Model): 34 missing responses dropped (e.g. 0m5WN4~1448, 9VUkAq~4562, At886V~3296, BEJI92~9143, CtOJ9t~2837, ...): NA response (missing observation) left out of the fit
#> See $provenance().
head(fa_deqms$get_contrasts())
#> # A tibble: 6 × 16
#>   modelName estimate_type contrast  protein_Id     diff avgAbd
#>   <chr>     <chr>         <chr>     <chr>         <dbl>  <dbl>
#> 1 deqms     observed      A_vs_Ctrl 0EfVhX~5954  2.72     23.2
#> 2 deqms     observed      A_vs_Ctrl 0m5WN4~1448  0.600    17.4
#> 3 deqms     observed      A_vs_Ctrl 7cbcrd~8305  2.59     27.0
#> 4 deqms     observed      A_vs_Ctrl 9VUkAq~4562  0.0679   19.4
#> 5 deqms     observed      A_vs_Ctrl At886V~3296 -1.01     19.1
#> 6 deqms     observed      A_vs_Ctrl BEJI92~9143 -0.873    20.9
#> # ℹ 10 more variables: std.error.unmoderated <dbl>, df.unmoderated <int>,
#> #   std.error <dbl>, statistic <dbl>, df <dbl>, p.value <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, sigma <dbl>, FDR <dbl>

istar_pep <- sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
lfqdata_pep <- LFQData$new(istar_pep$data, istar_pep$config)
lfqdata_pep <- lfqdata_pep$get_Transformer()$log2()$lfq
#> Column added : log2_abundance

fa_lmer <- build_contrast_analysis(
  lfqdata_pep,
  "~ group_ + (1 | peptide_Id) + (1 | sampleName)",
  contrasts,
  method = "lmer_nested"
)
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')
#> Warning: fit (Model): 4 fits failed (e.g. 7cbcrd~5725, 9VUkAq~4703, CGzoYe~2147, DoWup2~5896): Error: grouping factors must have > 1 sampled level
#> See $provenance().
#> fit (Model): 5 singular fits (e.g. 0EfVhX~0087, Fl4JiV~8625, HvIpHG~9079, JcKVfU~9653, SGIVBl~5782): coefficients not estimable or fewer than 2 residual df
#> fit (Model): 32 missing responses dropped (e.g. 0EfVhX~0087, BEJI92~5282, Fl4JiV~8625, HvIpHG~9079, SGIVBl~5782): NA response (missing observation) left out of the fit
#> See $provenance().
head(fa_lmer$get_contrasts())
#> # A tibble: 6 × 16
#>   modelName   estimate_type protein_Id  contrast       diff avgAbd
#>   <chr>       <chr>         <chr>       <chr>         <dbl>  <dbl>
#> 1 lmer_nested observed      0EfVhX~0087 A_vs_Ctrl -0.000832   4.34
#> 2 lmer_nested observed      BEJI92~5282 A_vs_Ctrl  0.322      4.22
#> 3 lmer_nested observed      Fl4JiV~8625 A_vs_Ctrl -0.0413     4.38
#> 4 lmer_nested observed      HvIpHG~9079 A_vs_Ctrl -0.372      4.40
#> 5 lmer_nested observed      JcKVfU~9653 A_vs_Ctrl -0.107      5.05
#> 6 lmer_nested observed      SGIVBl~5782 A_vs_Ctrl  0.0307     4.68
#> # ℹ 10 more variables: std.error.unmoderated <dbl>, df.unmoderated <dbl>,
#> #   std.error <dbl>, statistic <dbl>, df <dbl>, p.value <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, sigma <dbl>, FDR <dbl>

fa_ropeca <- build_contrast_analysis(lfqdata_pep, "~ group_", contrasts, method = "ropeca_nested")
#> fit (Model): 36 missing responses dropped (e.g. 0EfVhX~0087~ITLb4x1q, 0EfVhX~0087~ahQLlQY7, 0EfVhX~0087~dJkdz7so, 7cbcrd~5725~D5dQ4nKk, 9VUkAq~4703~eIC06D7g, ...): NA response (missing observation) left out of the fit
#> See $provenance().
head(fa_ropeca$get_contrasts())
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

fa_firth <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "firth")
head(fa_firth$get_contrasts())
#> # A tibble: 6 × 16
#>   modelName estimate_type protein_Id  contrast  sigma    df      diff   FDR
#>   <chr>     <chr>         <chr>       <chr>     <dbl> <int>     <dbl> <dbl>
#> 1 firth     observed      0EfVhX~5954 A_vs_Ctrl     1     9  4.27e-16     1
#> 2 firth     observed      0m5WN4~1448 A_vs_Ctrl     1     9 -2.20e+ 0     1
#> 3 firth     observed      7cbcrd~8305 A_vs_Ctrl     1     9  4.27e-16     1
#> 4 firth     observed      9VUkAq~4562 A_vs_Ctrl     1     9 -1.35e+ 0     1
#> 5 firth     observed      At886V~3296 A_vs_Ctrl     1     9  5.58e-16     1
#> 6 firth     observed      BEJI92~9143 A_vs_Ctrl     1     9 -1.35e+ 0     1
#> # ℹ 8 more variables: std.error <dbl>, statistic <dbl>, p.value <dbl>,
#> #   conf.low <dbl>, conf.high <dbl>, avgAbd <dbl>, std.error.unmoderated <dbl>,
#> #   df.unmoderated <int>
```
