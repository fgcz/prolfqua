# All modelling methods at a glance

## Purpose

[`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md)
selects the model with its `method` argument. This vignette runs every
method listed by
[`list_facades()`](https://fgcz.github.io/prolfqua/reference/list_facades.md)
on one small simulated dataset and shows one line of output per method.
Plots and a discussion of the methods are in the article [Comparing
modelling methods on one
dataset](https://fgcz.github.io/prolfqua/articles/ComparingMethods.html).

![build_contrast_analysis() is a factory: the method key selects the
model, and every model rolls out with the same
interface.](../reference/figures/prolfqua_factory.png)

build_contrast_analysis() is a factory: the method key selects the
model, and every model rolls out with the same interface.

## Data

We simulate peptide intensities for 20 proteins in three groups and
compare group `A` with the control.

``` r
library(prolfqua)
sim <- sim_lfq_data_peptide_config(Nprot = 20, seed = 42)
lfq_peptide <- LFQData$new(sim$data, sim$config)
lfq_peptide <- lfq_peptide$get_Transformer()$log2()$lfq
lfq_protein <- lfq_peptide$get_Aggregator()$aggregate()
lfq_limpa <- AggregateLimpa$new(lfq_peptide, "protein")$aggregate()
contrasts <- c(A_vs_Ctrl = "group_A - group_Ctrl")
```

Each method declares the input it needs:

- `aggregated`: protein-level input, here `lfq_protein` aggregated with
  median polish. `limpa` takes `lfq_limpa`, aggregated by
  `AggregateLimpa`, which adds the per-protein standard errors it needs.
- `nested`: peptide-level input, `lfq_peptide`. The method models the
  peptides and reports protein-level contrasts.

## Fit all methods

``` r
facades <- list_facades()
input_for <- function(method) {
  if (method == "limpa") {
    return(lfq_limpa)
  }
  if (facades[[method]]$needs == "nested") lfq_peptide else lfq_protein
}
model_for <- function(method) {
  if (method == "lmer_nested") "~ group_ + (1 | peptide_Id) + (1 | sampleName)" else "~ group_"
}
results <- lapply(names(facades), function(method) {
  fa <- build_contrast_analysis(input_for(method), model_for(method), contrasts, method = method)
  fa$get_contrasts()
})
names(results) <- names(facades)
```

## One line per method

``` r
about <- c(
  lm = "linear model per protein",
  lm_impute = "lm; groups without observations imputed at the limit of detection (LOD)",
  rlm = "robust linear model (MASS::rlm)",
  rfit = "rank-based regression (Rfit)",
  rfit_impute = "rfit with LOD imputation",
  limma = "limma, empirical-Bayes moderated variances",
  limma_impute = "limma with LOD imputation",
  limma_voom = "limma with vooma precision weights",
  limma_voom_impute = "limma_voom with LOD imputation",
  deqms = "limma with DEqMS peptide-count variance moderation",
  deqms_voom = "deqms with vooma precision weights",
  limpa = "limma on limpa DPC-quantified proteins",
  limpa_nested = "limpa, DPC quantification from peptides",
  firth = "Firth logistic regression, protein detection",
  binomial_nested = "quasibinomial model, peptide detection counts",
  binomial_firth_nested = "binomial_nested with Firth bias reduction",
  hurdle_nested = "lm_impute intensity + binomial_firth_nested detection",
  lmer_nested = "mixed model, peptide and sample random effects",
  ropeca_nested = "peptide fold changes summarised per protein (ROPECA)",
  ropecaV2_nested = "ROPECA with an exact, correlation-aware null for the median peptide z-score"
)
overview <- data.frame(
  method = names(results),
  input = vapply(facades, function(x) x$needs, character(1)),
  proteins = vapply(results, nrow, integer(1)),
  FDR_lt_0.1 = vapply(results, function(x) sum(x$FDR < 0.1, na.rm = TRUE), integer(1)),
  model = about[names(results)]
)
knitr::kable(overview, row.names = FALSE)
```

| method | input | proteins | FDR_lt_0.1 | model |
|:---|:---|---:|---:|:---|
| binomial_firth_nested | nested | 20 | 0 | binomial_nested with Firth bias reduction |
| binomial_nested | nested | 20 | 0 | quasibinomial model, peptide detection counts |
| deqms | aggregated | 20 | 14 | limma with DEqMS peptide-count variance moderation |
| deqms_voom | aggregated | 20 | 13 | deqms with vooma precision weights |
| firth | aggregated | 20 | 0 | Firth logistic regression, protein detection |
| hurdle_nested | nested | 20 | 14 | lm_impute intensity + binomial_firth_nested detection |
| limma | aggregated | 20 | 14 | limma, empirical-Bayes moderated variances |
| limma_impute | aggregated | 20 | 14 | limma with LOD imputation |
| limma_voom | aggregated | 20 | 14 | limma with vooma precision weights |
| limma_voom_impute | aggregated | 20 | 14 | limma_voom with LOD imputation |
| limpa | aggregated | 20 | 16 | limma on limpa DPC-quantified proteins |
| limpa_nested | nested | 20 | 16 | limpa, DPC quantification from peptides |
| lm | aggregated | 20 | 14 | linear model per protein |
| lm_impute | aggregated | 20 | 14 | lm; groups without observations imputed at the limit of detection (LOD) |
| lmer_nested | nested | 12 | 5 | mixed model, peptide and sample random effects |
| rfit | aggregated | 20 | 14 | rank-based regression (Rfit) |
| rfit_impute | aggregated | 20 | 14 | rfit with LOD imputation |
| rlm | aggregated | 20 | 14 | robust linear model (MASS::rlm) |
| ropeca_nested | nested | 20 | 14 | peptide fold changes summarised per protein (ROPECA) |
| ropecaV2_nested | nested | 20 | 16 | ROPECA with an exact, correlation-aware null for the median peptide z-score |

The detection models (`firth`, `binomial_nested`,
`binomial_firth_nested`) test differences in how often a protein or its
peptides are observed. The simulated data has few missing values, so
they find little.

## One protein, all methods

All methods return the same columns. Here is the protein with the
smallest `lm` p-value. For the detection models `diff` is a log odds
ratio, not a log2 fold change. This protein is detected in every sample,
so the detection models report `diff` 0 and a p-value of 1 or `NaN`.

``` r
top <- results$lm$protein_Id[which.min(results$lm$p.value)]
one <- lapply(results, function(x) x[x$protein_Id == top, c("diff", "p.value", "FDR")])
one <- dplyr::bind_rows(one, .id = "method")
knitr::kable(one, digits = 3)
```

| method                |   diff | p.value |   FDR |
|:----------------------|-------:|--------:|------:|
| binomial_firth_nested |  0.000 |   1.000 | 1.000 |
| binomial_nested       |  0.000 |   1.000 | 1.000 |
| deqms                 | -0.469 |   0.000 | 0.000 |
| deqms_voom            | -0.469 |   0.000 | 0.000 |
| firth                 |  0.000 |   1.000 | 1.000 |
| hurdle_nested         | -0.469 |   0.000 | 0.000 |
| limma                 | -0.469 |   0.000 | 0.000 |
| limma_impute          | -0.469 |   0.000 | 0.000 |
| limma_voom            | -0.469 |   0.000 | 0.000 |
| limma_voom_impute     | -0.469 |   0.000 | 0.000 |
| limpa                 | -0.272 |   0.000 | 0.000 |
| limpa_nested          | -0.272 |   0.000 | 0.000 |
| lm                    | -0.469 |   0.000 | 0.000 |
| lm_impute             | -0.469 |   0.000 | 0.000 |
| lmer_nested           | -0.276 |   0.001 | 0.003 |
| rfit                  | -0.463 |   0.000 | 0.000 |
| rfit_impute           | -0.463 |   0.000 | 0.000 |
| rlm                   | -0.469 |   0.000 | 0.000 |
| ropeca_nested         | -0.495 |   0.000 | 0.000 |
| ropecaV2_nested       | -0.495 |   0.000 | 0.000 |

## Session info

``` r
sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 26.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] prolfqua_1.9.4   BiocStyle_2.40.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_1.2.1       dplyr_1.2.1            farver_2.1.2          
#>  [4] S7_0.2.2               fastmap_1.2.0          rpart_4.1.27          
#>  [7] digest_0.6.39          lifecycle_1.0.5        Rfit_0.27.0           
#> [10] survival_3.8-12        statmod_1.5.2          magrittr_2.0.5        
#> [13] compiler_4.6.2         rlang_1.3.0            sass_0.4.10           
#> [16] progress_1.2.3         tools_4.6.2            yaml_2.3.12           
#> [19] data.table_1.18.6.1    limpa_1.4.2            knitr_1.52            
#> [22] prettyunits_1.2.0      htmlwidgets_1.6.4      RColorBrewer_1.1-3    
#> [25] numDeriv_2016.8-1.1    withr_3.0.3            purrr_1.2.2           
#> [28] desc_1.4.3             nnet_7.3-21            grid_4.6.2            
#> [31] jomo_2.7-6             mice_3.19.0            ggplot2_4.0.4         
#> [34] scales_1.4.0           iterators_1.0.14       MASS_7.3-66           
#> [37] cli_3.6.6              rmarkdown_2.32         crayon_1.5.3          
#> [40] ragg_1.5.2             reformulas_0.4.4       generics_0.1.4        
#> [43] otel_0.3.0             minqa_1.2.8            cachem_1.1.0          
#> [46] stringr_1.6.0          operator.tools_1.6.3.1 splines_4.6.2         
#> [49] BiocManager_1.30.27    vctrs_0.7.3            boot_1.3-32           
#> [52] glmnet_5.2             Matrix_1.7-6           jsonlite_2.0.0        
#> [55] bookdown_0.48          hms_1.1.4              mitml_0.4-5           
#> [58] ggrepel_0.9.8          systemfonts_1.3.3      foreach_1.5.2         
#> [61] limma_3.68.5           tidyr_1.3.2            jquerylib_0.1.4       
#> [64] glue_1.8.1             pkgdown_2.2.1          nloptr_2.2.1          
#> [67] pan_2.0                codetools_0.2-20       stringi_1.8.9         
#> [70] gtable_0.3.6           shape_1.4.6.1          lmerTest_3.2-1        
#> [73] lme4_2.1-0             tibble_3.3.1           pillar_1.11.1         
#> [76] htmltools_0.5.9        R6_2.6.1               textshaping_1.0.5     
#> [79] Rdpack_2.6.6           formula.tools_1.7.1    evaluate_1.0.5        
#> [82] lattice_0.23-1         rbibutils_2.4.1        backports_1.5.1       
#> [85] broom_1.0.13           bslib_0.12.0           Rcpp_1.1.2            
#> [88] nlme_3.1-171           mgcv_1.9-4             logistf_1.26.1        
#> [91] xfun_0.62              fs_2.1.0               pkgconfig_2.0.3
```
