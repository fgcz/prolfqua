# Filling the gaps: model-based imputation for heatmaps and PCA

## The problem

Heatmaps, PCA and clustering need a complete data matrix. Label-free
data are not complete: a protein is often not observed in some samples,
and often not in any sample of one group. While such proteins are
frequently below the limit of detection (LOD), complete absence can also
reflect true biological on/off expression, which makes them candidates
of high biological interest. PCA drops every protein with a missing
value, and with it exactly these proteins.

``` r
library(prolfqua)
library(ggplot2)
sim <- sim_lfq_data_protein_config(n_proteins = 100, weight_missing = 0.5, seed = 1234)
lfq <- LFQData$new(sim$data, sim$config)

wide <- lfq$data_wide(as.matrix = TRUE)$data
empty_group <- lfq$data_long() |>
  dplyr::group_by(protein_Id, group_) |>
  dplyr::summarise(empty = all(is.na(abundance)), .groups = "drop") |>
  dplyr::filter(empty)
c(
  proteins = nrow(wide),
  complete = sum(stats::complete.cases(wide)),
  with_empty_group = dplyr::n_distinct(empty_group$protein_Id)
)
```

    ##         proteins         complete with_empty_group 
    ##              100               23                5

The proteins with an empty group:

``` r
lfq$data_long() |>
  dplyr::filter(protein_Id %in% empty_group$protein_Id) |>
  dplyr::select(protein_Id, sampleName, abundance) |>
  tidyr::pivot_wider(names_from = sampleName, values_from = abundance) |>
  knitr::kable(digits = 1)
```

| protein_Id  | A_V1 | A_V2 | A_V3 | A_V4 | B_V1 | B_V2 | B_V3 | B_V4 | Ctrl_V1 | Ctrl_V2 | Ctrl_V3 | Ctrl_V4 |
|:------------|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|--------:|--------:|--------:|--------:|
| cnCfON~9769 |   NA | 15.9 | 14.3 | 14.4 |   NA |   NA |   NA |   NA |      NA |    15.8 |    17.5 |      NA |
| DoWup2~8058 | 17.9 |   NA | 18.9 |   NA |   NA |   NA |   NA |   NA |      NA |      NA |      NA |      17 |
| LMQBev~8107 |   NA | 16.0 |   NA |   NA |   NA | 15.8 | 16.0 | 15.8 |      NA |      NA |      NA |      NA |
| umjHSH~2090 |   NA |   NA |   NA |   NA |   NA |   NA | 15.8 | 16.5 |      NA |      NA |    17.4 |      NA |
| ZHBkZv~0394 |   NA |   NA |   NA |   NA |   NA |   NA |   NA | 13.3 |      NA |      NA |    12.6 |      NA |

## Fold changes are a separate question

To test for differences you do not need to impute anything yourself:
`build_contrast_analysis(method = "lm_impute")` also estimates fold
changes for proteins with an empty group, refitting them at the LOD and
borrowing their variance from the other proteins. See [Comparing two
groups](https://fgcz.github.io/prolfqua/articles/Comparing2Groups.html)
and [Comparing modelling methods on one
dataset](https://fgcz.github.io/prolfqua/articles/ComparingMethods.html).
Here we reuse the models of that fit to fill the data matrix.

## Model-based imputation

[`impute_from_model()`](https://fgcz.github.io/prolfqua/reference/impute_from_model.md)
fills every missing cell with the prediction of its protein model and
keeps the observed values. The summary tells which route each protein
took:

- `complete`: no missing values, nothing to fill
- `fitted`: the model fitted on the observed values predicts the missing
  cells from their group
- `lod_refit`: a protein whose fit failed, typically because a group is
  empty, is refitted with its missing values set to the LOD and observed
  values below the LOD raised to it; the empty group gets a value at the
  LOD

``` r
contrasts <- c(AvsC = "group_A - group_Ctrl", BvsC = "group_B - group_Ctrl")
fa <- build_contrast_analysis(lfq, "~ group_", contrasts, method = "lm_impute")
filled <- impute_from_model(fa$model, lfq)
table(filled$summary$route)
```

    ## 
    ##  complete    fitted lod_refit 
    ##        23        70         7

`MissingHelpers` estimates the LOD from the data. Groups in which a
protein is observed in only one sample, or in the fewest samples, are
the closest the data comes to the detection limit; the LOD is the median
of their mean intensities. Pass `lod` to
[`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md)
to use a LOD of your own.

``` r
MissingHelpers$new(lfq$data_long(), lfq$get_config())$get_lod()
```

    ## [1] 15.97769

The proteins with an empty group after imputation; compare with the
table above:

``` r
filled$lfqdata$data_long() |>
  dplyr::filter(protein_Id %in% empty_group$protein_Id) |>
  dplyr::select(protein_Id, sampleName, abundance) |>
  tidyr::pivot_wider(names_from = sampleName, values_from = abundance) |>
  knitr::kable(digits = 1)
```

| protein_Id  | A_V1 | A_V2 | A_V3 | A_V4 | B_V1 | B_V2 | B_V3 | B_V4 | Ctrl_V1 | Ctrl_V2 | Ctrl_V3 | Ctrl_V4 |
|:------------|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|--------:|--------:|--------:|--------:|
| cnCfON~9769 | 16.0 | 15.9 | 14.3 | 14.4 | 16.0 | 16.0 | 16.0 | 16.0 |    16.5 |    15.8 |    17.5 |    16.5 |
| DoWup2~8058 | 17.9 | 17.6 | 18.9 | 17.6 | 16.0 | 16.0 | 16.0 | 16.0 |    16.4 |    16.4 |    16.4 |    17.0 |
| LMQBev~8107 | 16.0 | 16.0 | 16.0 | 16.0 | 16.0 | 15.8 | 16.0 | 15.8 |    16.0 |    16.0 |    16.0 |    16.0 |
| umjHSH~2090 | 16.0 | 16.0 | 16.0 | 16.0 | 16.2 | 16.2 | 15.8 | 16.5 |    16.9 |    16.9 |    17.4 |    16.9 |
| ZHBkZv~0394 | 16.0 | 16.0 | 16.0 | 16.0 | 16.0 | 16.0 | 16.0 | 13.3 |    16.0 |    16.0 |    12.6 |    16.0 |

## Heatmap before and after

We cluster the proteins once, on the imputed data, and draw both
heatmaps with that row order, so each row shows the same protein before
and after.

``` r
after <- filled$lfqdata$data_wide(as.matrix = TRUE)
row_order <- rownames(after$data)[hclust(dist(t(scale(t(after$data)))))$order]
```

``` r
ht_before <- lfq$get_Plotter()$heatmap(row_order = row_order)
ht_before
```

![Before imputation: missing cells are
grey.](ModelBasedImputation_files/figure-html/heatmapBefore-1.png)

Before imputation: missing cells are grey.

``` r
ht_after <- filled$lfqdata$get_Plotter()$heatmap(row_order = row_order)
ht_after
```

![After imputation: every cell holds a
value.](ModelBasedImputation_files/figure-html/heatmapAfter-1.png)

After imputation: every cell holds a value.

## PCA before and after

Before imputation PCA can only use the proteins without missing values;
after imputation it uses all proteins.

``` r
n_complete <- sum(filled$summary$route == "complete")
gridExtra::grid.arrange(
  lfq$get_Plotter()$pca() +
    ggtitle(paste0("Before imputation: ", n_complete, " complete proteins")),
  filled$lfqdata$get_Plotter()$pca() +
    ggtitle(paste0("After imputation: ", nrow(filled$summary), " proteins")),
  ncol = 1
)
```

![PCA before (top) and after (bottom)
imputation.](ModelBasedImputation_files/figure-html/pca-1.png)

PCA before (top) and after (bottom) imputation.

## Caveats

After imputation the groups separate more cleanly, in the heatmap and in
the PCA. Part of that is circular: each filled value is the model’s
estimate for its group, so imputation adds the group structure of the
model to the data. The filled values also carry no noise. Use the filled
matrix to look at the data; for tests use the contrasts of `lm_impute`,
which account for the imputation.

## Session Info

``` r
sessionInfo()
```

    ## R version 4.6.2 (2026-10-29)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 26.04.1 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
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
    ## [1] ggplot2_4.0.4    prolfqua_1.9.9   BiocStyle_2.40.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] shape_1.4.6.1         circlize_0.4.18       gtable_0.3.6         
    ##  [4] rjson_0.2.23          xfun_0.62             bslib_0.12.0         
    ##  [7] htmlwidgets_1.6.4     GlobalOptions_0.1.4   ggrepel_0.9.8        
    ## [10] vctrs_0.7.3           tools_4.6.2           generics_0.1.4       
    ## [13] stats4_4.6.2          parallel_4.6.2        tibble_3.3.1         
    ## [16] cluster_2.1.8.3       pkgconfig_2.0.3       RColorBrewer_1.1-3   
    ## [19] S4Vectors_0.50.3      S7_1.0.0              desc_1.4.3           
    ## [22] lifecycle_1.0.5       compiler_4.6.2        farver_2.1.2         
    ## [25] stringr_1.6.0         textshaping_1.1.0     progress_1.2.3       
    ## [28] statmod_1.5.2         codetools_0.2-20      ComplexHeatmap_2.28.0
    ## [31] clue_0.3-68           htmltools_0.5.9       sass_0.5.0           
    ## [34] yaml_2.3.12           pillar_1.11.1         pkgdown_2.2.1        
    ## [37] crayon_1.5.3          jquerylib_0.1.4       tidyr_1.3.2          
    ## [40] MASS_7.3-66           cachem_1.1.0          limma_3.68.5         
    ## [43] iterators_1.0.14      foreach_1.5.2         tidyselect_1.2.1     
    ## [46] digest_0.6.39         stringi_1.8.9         dplyr_1.2.1          
    ## [49] purrr_1.2.2           bookdown_0.48         labeling_0.4.3       
    ## [52] fastmap_1.2.0         grid_4.6.2            colorspace_2.1-4     
    ## [55] cli_3.6.6             magrittr_2.0.5        withr_3.0.3          
    ## [58] prettyunits_1.2.0     scales_1.4.0          rmarkdown_2.32       
    ## [61] matrixStats_1.5.0     otel_0.3.0            gridExtra_2.3.1      
    ## [64] ragg_1.5.2            png_0.1-9             GetoptLong_1.1.1     
    ## [67] hms_1.1.4             evaluate_1.0.5        knitr_1.52           
    ## [70] IRanges_2.46.0        doParallel_1.0.17     rlang_1.3.0          
    ## [73] Rcpp_1.1.2            glue_1.8.1            BiocManager_1.30.27  
    ## [76] BiocGenerics_0.58.1   jsonlite_2.0.0        R6_2.6.1             
    ## [79] systemfonts_1.3.4     fs_2.1.0
