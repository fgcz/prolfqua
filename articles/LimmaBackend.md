# Limma Backend for prolfqua

## Purpose

This vignette demonstrates the **limma backend** for prolfqua’s
modelling pipeline. The limma backend is a drop-in alternative to the
per-protein `lm` + `squeezeVar` pipeline. Both approaches fit linear
models and perform empirical Bayes variance shrinkage, but:

- **prolfqua default** (`strategy_lm` + `Contrasts` +
  `ContrastsModerated`): fits one
  [`lm()`](https://rdrr.io/r/stats/lm.html) per protein, then applies
  prolfqua’s own `squeezeVarRob` moderation.
- **limma backend** (`strategy_limma` + `ContrastsLimma`): uses limma’s
  matrix-based `lmFit` + `contrasts.fit` + `eBayes` pipeline directly.

The user-facing API is identical: same formula specification, same
contrast specification, same downstream methods (`get_contrasts`,
`get_Plotter`, `to_wide`, `merge_contrasts_results`).

## Data preparation

We use simulated protein-level data with 3 groups (A, B, Ctrl) and some
missing values.

``` r
library(prolfqua)
library(dplyr)

istar <- sim_lfq_data_protein_config(Nprot = 100, weight_missing = 0.3)
lfqdata <- LFQData$new(istar$data, istar$config)
lfqdata$remove_small_intensities()

lt <- lfqdata$get_Transformer()
transformed <- lt$log2()$robscale()$lfq
transformed$rename_response("transformedIntensity")
transformed$factors()
```

    ## # A tibble: 12 × 3
    ##    sample  sampleName group_
    ##    <chr>   <chr>      <chr> 
    ##  1 A_V1    A_V1       A     
    ##  2 A_V2    A_V2       A     
    ##  3 A_V3    A_V3       A     
    ##  4 A_V4    A_V4       A     
    ##  5 B_V1    B_V1       B     
    ##  6 B_V2    B_V2       B     
    ##  7 B_V3    B_V3       B     
    ##  8 B_V4    B_V4       B     
    ##  9 Ctrl_V1 Ctrl_V1    Ctrl  
    ## 10 Ctrl_V2 Ctrl_V2    Ctrl  
    ## 11 Ctrl_V3 Ctrl_V3    Ctrl  
    ## 12 Ctrl_V4 Ctrl_V4    Ctrl

## Define contrasts

The same contrast specification is used for both backends.

``` r
contr_spec <- c(
  "AvsCtrl" = "group_A - group_Ctrl",
  "BvsCtrl" = "group_B - group_Ctrl",
  "AvsB"    = "group_A - group_B"
)
```

## prolfqua default pipeline

``` r
mod_lm <- build_model(transformed, strategy_lm("transformedIntensity ~ group_"))
contr_lm <- prolfqua::Contrasts$new(mod_lm, contr_spec)
contr_moderated <- ContrastsModerated$new(contr_lm)
res_moderated <- contr_moderated$get_contrasts()
```

## Limma backend pipeline

The only difference is using `strategy_limma` + `build_model_limma` +
`ContrastsLimma`.

``` r
strat_limma <- strategy_limma("transformedIntensity ~ group_")
mod_limma <- build_model_limma(transformed, strat_limma)
contr_limma <- ContrastsLimma$new(mod_limma, contr_spec)
res_limma <- contr_limma$get_contrasts()
```

### Model-level inspection

The `ModelLimma` object supports the same inspection methods as `Model`.

``` r
mod_limma$get_anova() |> head()
```

    ## # A tibble: 6 × 5
    ##   protein_Id  F.value       p.value factor                     FDR
    ##   <chr>         <dbl>         <dbl> <chr>                    <dbl>
    ## 1 0EfVhX~3967  6.37   0.00377       group_B+group_Ctrl 0.0133     
    ## 2 0m5WN4~6025 28.8    0.00000000998 group_B+group_Ctrl 0.000000329
    ## 3 0YSKpy~2865  1.24   0.298         group_B+group_Ctrl 0.434      
    ## 4 3QLHfm~8938  7.57   0.00146       group_B+group_Ctrl 0.00805    
    ## 5 3QYop0~7543  0.0688 0.934         group_B+group_Ctrl 0.963      
    ## 6 76k03k~7094  3.12   0.0539        group_B+group_Ctrl 0.133

``` r
mod_limma$get_coefficients() |> head()
```

    ## # A tibble: 6 × 6
    ##   protein_Id  factor      Estimate Std..Error t.value Pr...t..
    ##   <chr>       <chr>          <dbl>      <dbl>   <dbl>    <dbl>
    ## 1 0EfVhX~3967 (Intercept)     4.47     0.0359    124. 6.83e-57
    ## 2 0m5WN4~6025 (Intercept)     4.10     0.0407    101. 6.81e-54
    ## 3 0YSKpy~2865 (Intercept)     4.16     0.0368    113. 4.05e-56
    ## 4 3QLHfm~8938 (Intercept)     4.55     0.0340    134. 2.07e-60
    ## 5 3QYop0~7543 (Intercept)     4.51     0.0350    129. 1.08e-59
    ## 6 76k03k~7094 (Intercept)     4.67     0.0355    132. 4.16e-60

``` r
mod_limma$anova_histogram()
```

    ## $plot

![](LimmaBackend_files/figure-html/model_inspection-1.png)

    ## 
    ## $name
    ## [1] "Anova_p.values_limma.pdf"

## Comparing results

### Fold changes are identical

Both backends estimate the same linear model, so fold changes (log2 FC)
are identical.

``` r
merged <- inner_join(
  select(res_limma, protein_Id, contrast, diff_limma = diff),
  select(res_moderated, protein_Id, contrast, diff_moderated = diff),
  by = c("protein_Id", "contrast")
)

cor(merged$diff_limma, merged$diff_moderated, use = "complete.obs")
```

    ## [1] 1

``` r
plot(merged$diff_limma, merged$diff_moderated,
     xlab = "limma logFC", ylab = "prolfqua moderated logFC",
     pch = 20, col = rgb(0, 0, 0, 0.3))
abline(0, 1, col = "red")
```

![Fold changes: limma vs. prolfqua moderated. Points lie on the
diagonal.](LimmaBackend_files/figure-html/plot_fc-1.png)

Fold changes: limma vs. prolfqua moderated. Points lie on the diagonal.

### P-values differ (different shrinkage)

The p-values differ because limma and prolfqua use different eBayes
implementations, but they are highly correlated.

``` r
merged_p <- inner_join(
  select(res_limma, protein_Id, contrast, p_limma = p.value),
  select(res_moderated, protein_Id, contrast, p_moderated = p.value),
  by = c("protein_Id", "contrast")
)

cor(-log10(merged_p$p_limma), -log10(merged_p$p_moderated), use = "complete.obs")
```

    ## [1] 0.984774

``` r
plot(-log10(merged_p$p_limma), -log10(merged_p$p_moderated),
     xlab = "-log10(p) limma", ylab = "-log10(p) prolfqua moderated",
     pch = 20, col = rgb(0, 0, 0, 0.3))
abline(0, 1, col = "red")
```

![P-values: limma vs. prolfqua moderated. High correlation but not
identical.](LimmaBackend_files/figure-html/plot_pval-1.png)

P-values: limma vs. prolfqua moderated. High correlation but not
identical.

## Volcano plots

``` r
pl_limma <- contr_limma$get_Plotter()
pl_moderated <- contr_moderated$get_Plotter()

gridExtra::grid.arrange(
  pl_limma$volcano()$FDR + ggplot2::ggtitle("limma"),
  pl_moderated$volcano()$FDR + ggplot2::ggtitle("prolfqua moderated"),
  ncol = 1
)
```

![Volcano plots from both
backends.](LimmaBackend_files/figure-html/volcano-1.png)

Volcano plots from both backends.

## Merging with missing value imputation

`ContrastsLimma` integrates seamlessly with `ContrastsMissing` and
`merge_contrasts_results`.

``` r
mC <- ContrastsMissing$new(lfqdata = transformed, contrasts = contr_spec)
merged_result <- merge_contrasts_results(prefer = contr_limma, add = mC)

plotter <- merged_result$merged$get_Plotter()
plotter$volcano()$FDR
```

![](LimmaBackend_files/figure-html/merge_missing-1.png)

## Wide format export

``` r
wide <- contr_limma$to_wide()
head(wide)
```

    ## # A tibble: 6 × 13
    ##   protein_Id diff.AvsCtrl diff.BvsCtrl diff.AvsB p.value.AvsCtrl p.value.BvsCtrl
    ##   <chr>             <dbl>        <dbl>     <dbl>           <dbl>           <dbl>
    ## 1 0EfVhX~39…      0.0875      -0.108      0.195        0.118            0.0555  
    ## 2 0m5WN4~60…     -0.235        0.172     -0.407        0.0000768        0.00258 
    ## 3 0YSKpy~28…      0.00148     -0.0703     0.0718       0.979            0.202   
    ## 4 3QLHfm~89…     -0.0343      -0.177      0.142        0.480            0.000641
    ## 5 3QYop0~75…     -0.0180      -0.00598   -0.0120       0.717            0.904   
    ## 6 76k03k~70…      0.0169      -0.0991     0.116        0.738            0.0544  
    ## # ℹ 7 more variables: p.value.AvsB <dbl>, FDR.AvsCtrl <dbl>, FDR.BvsCtrl <dbl>,
    ## #   FDR.AvsB <dbl>, statistic.AvsCtrl <dbl>, statistic.BvsCtrl <dbl>,
    ## #   statistic.AvsB <dbl>

## Two-factor design

The limma backend handles multi-factor designs in the same way.

``` r
dd <- sim_lfq_data_2Factor_config(Nprot = 100, with_missing = FALSE)
lProt2 <- LFQData$new(dd$data, dd$config)
lProt2$rename_response("transformedIntensity")

strat2 <- strategy_limma("transformedIntensity ~ Treatment + Background")
mod2 <- build_model_limma(lProt2, strat2)

Contr2 <- c(
  "A_vs_B" = "TreatmentA - TreatmentB",
  "X_vs_Z" = "BackgroundX - BackgroundZ"
)
contr2 <- ContrastsLimma$new(mod2, Contr2)

pl2 <- contr2$get_Plotter()
pl2$volcano()$FDR
```

![](LimmaBackend_files/figure-html/two_factor-1.png)

## strategy_limma options

The `strategy_limma` function accepts `trend` and `robust` arguments
that are passed to
[`limma::eBayes`](https://rdrr.io/pkg/limma/man/ebayes.html):

- `trend = TRUE`: allows the prior variance to depend on average
  intensity (useful when variance is intensity-dependent)
- `robust = TRUE`: uses a robust empirical Bayes procedure (resistant to
  outlier variances)

``` r
strat_robust <- strategy_limma(
  "transformedIntensity ~ group_",
  trend = TRUE,
  robust = TRUE
)
mod_robust <- build_model_limma(transformed, strat_robust)
contr_robust <- ContrastsLimma$new(mod_robust, contr_spec)
head(contr_robust$get_contrasts())
```

    ## # A tibble: 6 × 13
    ##   modelName protein_Id  contrast     diff      FDR std.error statistic   p.value
    ##   <chr>     <chr>       <chr>       <dbl>    <dbl>     <dbl>     <dbl>     <dbl>
    ## 1 limma     0EfVhX~3967 AvsCtrl   0.0875  0.400       0.0551    1.59   0.113    
    ## 2 limma     0m5WN4~6025 AvsCtrl  -0.235   0.000672    0.0569   -4.13   0.0000408
    ## 3 limma     0YSKpy~2865 AvsCtrl   0.00148 0.992       0.0641    0.0230 0.982    
    ## 4 limma     3QLHfm~8938 AvsCtrl  -0.0343  0.778       0.0476   -0.721  0.471    
    ## 5 limma     3QYop0~7543 AvsCtrl  -0.0180  0.931       0.0475   -0.380  0.704    
    ## 6 limma     76k03k~7094 AvsCtrl   0.0169  0.931       0.0447    0.378  0.706    
    ## # ℹ 5 more variables: sigma <dbl>, df <dbl>, conf.low <dbl>, conf.high <dbl>,
    ## #   avgAbd <dbl>

## Session Info

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] dplyr_1.2.0    prolfqua_1.5.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] utf8_1.2.6          tidyr_1.3.2         plotly_4.12.0      
    ##  [4] sass_0.4.10         generics_0.1.4      stringi_1.8.7      
    ##  [7] hms_1.1.4           digest_0.6.39       magrittr_2.0.4     
    ## [10] evaluate_1.0.5      grid_4.5.2          RColorBrewer_1.1-3 
    ## [13] fastmap_1.2.0       plyr_1.8.9          jsonlite_2.0.0     
    ## [16] progress_1.2.3      ggrepel_0.9.7       limma_3.66.0       
    ## [19] gridExtra_2.3       httr_1.4.8          purrr_1.2.1        
    ## [22] viridisLite_0.4.3   scales_1.4.0        UpSetR_1.4.0       
    ## [25] lazyeval_0.2.2      textshaping_1.0.4   jquerylib_0.1.4    
    ## [28] cli_3.6.5           crayon_1.5.3        rlang_1.1.7        
    ## [31] withr_3.0.2         cachem_1.1.0        yaml_2.3.12        
    ## [34] otel_0.2.0          tools_4.5.2         ggplot2_4.0.2      
    ## [37] forcats_1.0.1       vctrs_0.7.1         R6_2.6.1           
    ## [40] lifecycle_1.0.5     fs_1.6.6            htmlwidgets_1.6.4  
    ## [43] MASS_7.3-65         ragg_1.5.0          pkgconfig_2.0.3    
    ## [46] desc_1.4.3          pkgdown_2.2.0       pillar_1.11.1      
    ## [49] bslib_0.10.0        gtable_0.3.6        glue_1.8.0         
    ## [52] data.table_1.18.2.1 Rcpp_1.1.1          statmod_1.5.1      
    ## [55] systemfonts_1.3.1   xfun_0.56           tibble_3.3.1       
    ## [58] tidyselect_1.2.1    knitr_1.51          farver_2.1.2       
    ## [61] htmltools_0.5.9     labeling_0.4.3      rmarkdown_2.30     
    ## [64] pheatmap_1.0.13     compiler_4.5.2      prettyunits_1.2.0  
    ## [67] S7_0.2.1
