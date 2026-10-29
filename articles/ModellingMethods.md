# All modelling methods at a glance

## Purpose

[`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md)
selects the model with its `method` argument. This vignette runs every
method listed by
[`list_facades()`](https://fgcz.github.io/prolfqua/reference/list_facades.md)
on one small simulated dataset and shows one line of output per method.
Plots and a discussion of the methods are in the article [Contrast
facades with parallel
designs](https://fgcz.github.io/prolfqua/articles/ContrastFacades.html).

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

- `same`: protein-level input, here `lfq_protein` aggregated with median
  polish. `limpa` takes `lfq_limpa`, aggregated by `AggregateLimpa`,
  which adds the per-protein standard errors it needs.
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
  firth_nested = "Firth logistic regression, peptide detection",
  binomial_nested = "quasibinomial model, peptide detection counts",
  binomial_firth_nested = "binomial_nested with Firth bias reduction",
  hurdle_nested = "lm_impute intensity + binomial_firth_nested detection",
  lmer_nested = "mixed model, peptide and sample random effects",
  ropeca_nested = "peptide fold changes summarised per protein (ROPECA)"
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

| method                | input  | proteins | FDR_lt_0.1 | model                                                                   |
|:----------------------|:-------|---------:|-----------:|:------------------------------------------------------------------------|
| binomial_firth_nested | nested |       20 |          0 | binomial_nested with Firth bias reduction                               |
| binomial_nested       | nested |       20 |          0 | quasibinomial model, peptide detection counts                           |
| deqms                 | same   |       20 |         14 | limma with DEqMS peptide-count variance moderation                      |
| deqms_voom            | same   |       20 |         13 | deqms with vooma precision weights                                      |
| firth                 | same   |       20 |          0 | Firth logistic regression, protein detection                            |
| firth_nested          | nested |       20 |          0 | Firth logistic regression, peptide detection                            |
| hurdle_nested         | nested |       20 |         14 | lm_impute intensity + binomial_firth_nested detection                   |
| limma                 | same   |       20 |         14 | limma, empirical-Bayes moderated variances                              |
| limma_impute          | same   |       20 |         14 | limma with LOD imputation                                               |
| limma_voom            | same   |       20 |         14 | limma with vooma precision weights                                      |
| limma_voom_impute     | same   |       20 |         14 | limma_voom with LOD imputation                                          |
| limpa                 | same   |       20 |         13 | limma on limpa DPC-quantified proteins                                  |
| limpa_nested          | nested |       20 |         13 | limpa, DPC quantification from peptides                                 |
| lm                    | same   |       20 |         14 | linear model per protein                                                |
| lm_impute             | same   |       20 |         14 | lm; groups without observations imputed at the limit of detection (LOD) |
| lmer_nested           | nested |       12 |          5 | mixed model, peptide and sample random effects                          |
| rfit                  | same   |       20 |         14 | rank-based regression (Rfit)                                            |
| rfit_impute           | same   |       20 |         14 | rfit with LOD imputation                                                |
| rlm                   | same   |       20 |         14 | robust linear model (MASS::rlm)                                         |
| ropeca_nested         | nested |       20 |         14 | peptide fold changes summarised per protein (ROPECA)                    |

The detection models (`firth`, `firth_nested`, `binomial_nested`,
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
| binomial_firth_nested |  0.000 |     NaN |   NaN |
| binomial_nested       |  0.000 |     NaN |   NaN |
| deqms                 | -0.469 |   0.000 | 0.000 |
| deqms_voom            | -0.469 |   0.000 | 0.000 |
| firth                 |  0.000 |   1.000 | 1.000 |
| firth_nested          |  0.000 |   1.000 | 1.000 |
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

## Session info

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
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
#> [1] prolfqua_1.9.0   BiocStyle_2.38.0
#> 
#> loaded via a namespace (and not attached):
#>   [1] Rdpack_2.6.6           gridExtra_2.3.1        rlang_1.3.0           
#>   [4] magrittr_2.0.5         clue_0.3-68            GetoptLong_1.1.1      
#>   [7] otel_0.2.0             matrixStats_1.5.0      compiler_4.5.2        
#>  [10] mgcv_1.9-3             png_0.1-9              systemfonts_1.3.2     
#>  [13] vctrs_0.7.3            pkgconfig_2.0.3        shape_1.4.6.1         
#>  [16] crayon_1.5.3           fastmap_1.2.0          backports_1.5.1       
#>  [19] rmarkdown_2.32         nloptr_2.2.1           ragg_1.5.2            
#>  [22] UpSetR_1.4.1           purrr_1.2.2            xfun_0.62             
#>  [25] glmnet_5.2             jomo_2.7-6             logistf_1.26.1        
#>  [28] cachem_1.1.0           jsonlite_2.0.0         progress_1.2.3        
#>  [31] pan_2.0                Rfit_0.27.0            prettyunits_1.2.0     
#>  [34] broom_1.0.13           parallel_4.5.2         cluster_2.1.8.1       
#>  [37] R6_2.6.1               stringi_1.8.9          bslib_0.12.0          
#>  [40] RColorBrewer_1.1-3     limma_3.66.0           boot_1.3-32           
#>  [43] rpart_4.1.24           numDeriv_2016.8-1.1    jquerylib_0.1.4       
#>  [46] Rcpp_1.1.2             bookdown_0.48          iterators_1.0.14      
#>  [49] knitr_1.52             IRanges_2.44.0         Matrix_1.7-4          
#>  [52] splines_4.5.2          nnet_7.3-20            tidyselect_1.2.1      
#>  [55] yaml_2.3.12            doParallel_1.0.17      codetools_0.2-20      
#>  [58] lmerTest_3.2-1         lattice_0.22-7         tibble_3.3.1          
#>  [61] plyr_1.8.9             withr_3.0.3            S7_0.2.2              
#>  [64] evaluate_1.0.5         desc_1.4.3             survival_3.8-3        
#>  [67] circlize_0.4.18        pillar_1.11.1          BiocManager_1.30.27   
#>  [70] mice_3.19.0            foreach_1.5.2          stats4_4.5.2          
#>  [73] reformulas_0.4.4       plotly_4.12.1          generics_0.1.4        
#>  [76] hms_1.1.4              S4Vectors_0.48.1       ggplot2_4.0.3         
#>  [79] scales_1.4.0           minqa_1.2.8            glue_1.8.1            
#>  [82] tools_4.5.2            data.table_1.18.6.1    lme4_2.1-0            
#>  [85] forcats_1.0.1          fs_2.1.0               grid_4.5.2            
#>  [88] limpa_1.2.5            tidyr_1.3.2            rbibutils_2.4.1       
#>  [91] colorspace_2.1-4       nlme_3.1-168           formula.tools_1.7.1   
#>  [94] cli_3.6.6              textshaping_1.0.5      viridisLite_0.4.3     
#>  [97] ComplexHeatmap_2.26.1  dplyr_1.2.1            gtable_0.3.6          
#> [100] sass_0.4.10            digest_0.6.39          operator.tools_1.6.3.1
#> [103] BiocGenerics_0.56.0    ggrepel_0.9.8          rjson_0.2.23          
#> [106] htmlwidgets_1.6.4      farver_2.1.2           htmltools_0.5.9       
#> [109] pkgdown_2.2.1          lifecycle_1.0.5        httr_1.4.9            
#> [112] GlobalOptions_0.1.4    mitml_0.4-5            statmod_1.5.2         
#> [115] MASS_7.3-65
```
