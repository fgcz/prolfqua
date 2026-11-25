# Comparing Two Groups with prolfqua

## Purpose

This vignette demonstrates how two conditions, e.g., treatment
vs. control, can be compared and the differences statistically tested.
We use simulated protein abundances of an LFQ experiment with two
treatment groups (`A`, `B`) and a control group (`Ctrl`). After first
examining the data using QC plots and then normalizing the data, we
compare each treatment group with the control. The output of the
comparison is the difference in the mean intensities for quantified
proteins (log2 fold-change) in each group, along with statistical
parameters such as degrees of freedom, standard errors, p-value and the
FDR.

## Simulated protein abundances

``` r
library(prolfqua)
```

[`sim_lfq_data_protein_config()`](https://fgcz.github.io/prolfqua/reference/sim_lfq_data_protein_config.md)
simulates protein abundances in long format, one row per protein and
sample, and returns them together with a matching
`AnalysisConfiguration`. To read your own data, see the article
[Importing your
data](https://fgcz.github.io/prolfqua/articles/CreatingConfigurations.html).

``` r
xx <- prolfqua::sim_lfq_data_protein_config(n_proteins = 100, weight_missing = 0.5)
xx
```

    ## $data
    ## # A tibble: 1,200 × 7
    ##    sample sampleName group_ protein_Id  abundance qValue nrPeptides
    ##    <chr>  <chr>      <chr>  <chr>           <dbl>  <dbl>      <dbl>
    ##  1 A_V1   A_V1       A      0EfVhX~3967      20.4      0          2
    ##  2 A_V1   A_V1       A      0m5WN4~6025      16.4      0          7
    ##  3 A_V1   A_V1       A      0YSKpy~2865      NA       NA         NA
    ##  4 A_V1   A_V1       A      3QLHfm~8938      23.0      0          1
    ##  5 A_V1   A_V1       A      3QYop0~7543      23.1      0          7
    ##  6 A_V1   A_V1       A      76k03k~7094      25.5      0          4
    ##  7 A_V1   A_V1       A      7cbcrd~7351      27.1      0          2
    ##  8 A_V1   A_V1       A      7QuTub~1867      18.4      0          4
    ##  9 A_V1   A_V1       A      7soopj~5352      19.5      0          2
    ## 10 A_V1   A_V1       A      7zeekV~7127      20.6      0          8
    ## # ℹ 1,190 more rows
    ## 
    ## $config
    ## <AnalysisConfiguration>
    ##   Public:
    ##     annotation_vars: function () 
    ##     bin_resp: 
    ##     clone: function (deep = FALSE) 
    ##     factor_depth: 1
    ##     factor_keys: function () 
    ##     factor_keys_depth: function () 
    ##     factors: list
    ##     file_name: sample
    ##     get_response: function () 
    ##     hierarchy: list
    ##     hierarchy_depth: 1
    ##     hierarchy_keys: function () 
    ##     hierarchy_keys_depth: function (names = TRUE) 
    ##     id_vars: function () 
    ##     ident_q_value: qValue
    ##     ident_score: 
    ##     is_response_transformed: FALSE
    ##     min_peptides_protein: 2
    ##     norm_value: NULL
    ##     nr_children: nrPeptides
    ##     opt_mz: 
    ##     opt_rt: 
    ##     opt_se: 
    ##     pattern_contaminants: NULL
    ##     pattern_decoys: NULL
    ##     pop_response: function () 
    ##     sample_name: sampleName
    ##     sep: ~
    ##     set_response: function (col_name) 
    ##     value_vars: function () 
    ##     work_intensity: abundance

Create the `LFQData` class instance and remove zero or very small
intensities, which most quantification software reports for missing
values.

``` r
lfqdata <- prolfqua::LFQData$new(xx$data, xx$config)
lfqdata$remove_small_intensities()
lfqdata$factors()
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

You can always convert the data into wide format.

``` r
lfqdata$data_wide()$data[1:3,1:7]
```

    ## # A tibble: 3 × 7
    ##   protein_Id   A_V1  A_V2  A_V3  A_V4  B_V1  B_V2
    ##   <chr>       <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ## 1 0EfVhX~3967  20.4  23.8  22.3  21.5  NA    19.6
    ## 2 0m5WN4~6025  16.4  NA    NA    18.2  24.3  21.7
    ## 3 0YSKpy~2865  NA    17.7  16.7  17.7  18.2  17.5

### Visualization of not normalized data

After this first setting up of the analysis we show now how to normalize
the proteins and the effect of normalization. Furthermore we use some
functions to visualize the missing values in our data.

``` r
lfqplotter <- lfqdata$get_Plotter()
density_nn <- lfqplotter$intensity_distribution_density()
```

#### Visualization of missing data

``` r
lfqplotter$na_heatmap()
```

![Heatmap of missing protein intensities, black - missing, white -
present](Comparing2Groups_files/figure-html/makeMissingHeatmap-1.png)

Heatmap of missing protein intensities, black - missing, white - present

``` r
lfqdata$get_Summariser()$plot_missingness_per_group()
```

![\# of proteins with 0,1,...N missing
values](Comparing2Groups_files/figure-html/missignessPerGroup-1.png)

\# of proteins with 0,1,…N missing values

``` r
lfqplotter$missingness_histogram()
```

![Intensity distribution of proteins depending on \# of missing
values](Comparing2Groups_files/figure-html/missignessHistogram-1.png)

Intensity distribution of proteins depending on \# of missing values

#### Computing standard deviations, mean and CV.

Other important statistics such as coefficient of variation, means and
standard deviations can be easily calculated using the `get_Stats`
function and visualized with a violin plot.

``` r
stats <- lfqdata$get_Stats()
stats$violin()
```

![Violin plots of CVs in the different groups and among all
groups](Comparing2Groups_files/figure-html/PlotCVDistributions-1.png)

Violin plots of CVs in the different groups and among all groups

``` r
prolfqua::table_facade( stats$stats_quantiles()$wide, paste0("quantile of ",stats$stat ))
```

| probs |        A |       All |        B |     Ctrl |
|------:|---------:|----------:|---------:|---------:|
|  0.10 | 2.121145 |  3.416437 | 1.445998 | 1.499604 |
|  0.25 | 2.603169 |  4.316819 | 2.734033 | 2.135081 |
|  0.50 | 3.860492 |  5.077622 | 3.487257 | 3.744767 |
|  0.75 | 5.552086 |  7.229931 | 5.449944 | 5.701197 |
|  0.90 | 7.745671 | 10.211265 | 6.931290 | 8.396967 |

quantile of CV

``` r
stats$density_median()
```

![Distribution of CV's for top 50% and bottom 50% proteins by
intensity.](Comparing2Groups_files/figure-html/plotCVsplitbyMedianIntensity-1.png)

Distribution of CV’s for top 50% and bottom 50% proteins by intensity.

### Normalize protein intensities and show diagnostic plots

We normalize the data by $`\log_2`$ transforming and then $`z-scaling`$.

``` r
lt <- lfqdata$get_Transformer()
transformed <- lt$log2()$robscale()$lfq
transformed$is_transformed()
```

    ## [1] TRUE

``` r
pl <- transformed$get_Plotter()
density_norm <- pl$intensity_distribution_density()
```

``` r
gridExtra::grid.arrange(density_nn, density_norm)
```

![Distribution of intensities before and after
normalization.](Comparing2Groups_files/figure-html/showIntensityDistributions-1.png)

Distribution of intensities before and after normalization.

``` r
pl$pairs_smooth()
```

![Scatterplot
matrix](Comparing2Groups_files/figure-html/plotScatterMatrix-1.png)

Scatterplot matrix

    ## NULL

``` r
p <- pl$heatmap_cor()
```

``` r
p
```

![Heatmap, Rows - proteins, Columns -
samples](Comparing2Groups_files/figure-html/plotHeatmap-1.png)

Heatmap, Rows - proteins, Columns - samples

## Fit the model and compute contrasts

We compare each treatment group with the control.
[`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md)
fits one model per protein and computes the contrasts in one call. We
use the method `lm_impute`: a linear model per protein, with the protein
variances moderated by the empirical Bayes method so that proteins with
few observations borrow strength from all proteins. Proteins without
observations in one group are refitted with the missing values set to
the limit of detection (LOD). The article [Filling the
gaps](https://fgcz.github.io/prolfqua/articles/ModelBasedImputation.html)
uses these models to complete the data for heatmaps and PCA. The
vignette [All modelling methods at a
glance](https://fgcz.github.io/prolfqua/articles/ModellingMethods.md)
lists the other methods.

``` r
contr_spec <- c(
  "AvsC" = "group_A - group_Ctrl",
  "BvsC" = "group_B - group_Ctrl"
)
fa <- build_contrast_analysis(transformed, "~ group_", contr_spec, method = "lm_impute")
contrdf <- fa$get_contrasts()
```

The column `estimate_type` tells how each estimate was obtained:
`observed` from the measured intensities, `lod_imputed` from the refit
at the LOD.

``` r
dplyr::count(contrdf, contrast, estimate_type)
```

    ## # A tibble: 4 × 3
    ##   contrast estimate_type     n
    ##   <chr>    <chr>         <int>
    ## 1 AvsC     lod_imputed       7
    ## 2 AvsC     observed         93
    ## 3 BvsC     lod_imputed       7
    ## 4 BvsC     observed         93

``` r
plotter <- fa$get_Plotter()
plotter$volcano()$FDR
```

![Volcano plots of the two
contrasts.](Comparing2Groups_files/figure-html/plotVolcanos-1.png)

Volcano plots of the two contrasts.

``` r
plotter$ma_plot()
```

![MA plot showing the dependency of mean abundance with respect to the
difference](Comparing2Groups_files/figure-html/showMAplot-1.png)

MA plot showing the dependency of mean abundance with respect to the
difference

The proteins with FDR \< 0.1 in at least one contrast:

``` r
sig <- dplyr::filter(contrdf, FDR < 0.1)
signif <- transformed$get_copy()
signif$set_data(signif$data_long() |> dplyr::filter(protein_Id %in% sig$protein_Id))
hmSig <- signif$get_Plotter()$heatmap()
```

``` r
hmSig
```

![Heatmap of the proteins with FDR \< 0.1 in at least one
contrast](Comparing2Groups_files/figure-html/showSigHeatmap-1.png)

Heatmap of the proteins with FDR \< 0.1 in at least one contrast

``` r
myProteinIDS <- head(unique(contrdf$protein_Id), 3)
dplyr::filter(contrdf, protein_Id %in% myProteinIDS)
```

    ## # A tibble: 6 × 16
    ##   modelName estimate_type protein_Id  contrast    diff avgAbd
    ##   <chr>     <chr>         <chr>       <chr>      <dbl>  <dbl>
    ## 1 lm_impute observed      0EfVhX~3967 AvsC      0.0808   4.43
    ## 2 lm_impute observed      0m5WN4~6025 AvsC     -0.225    4.23
    ## 3 lm_impute observed      0YSKpy~2865 AvsC      0.0192   4.16
    ## 4 lm_impute observed      0EfVhX~3967 BvsC     -0.119    4.33
    ## 5 lm_impute observed      0m5WN4~6025 BvsC      0.162    4.42
    ## 6 lm_impute observed      0YSKpy~2865 BvsC     -0.0640   4.11
    ## # ℹ 10 more variables: std.error.unmoderated <dbl>, df.unmoderated <dbl>,
    ## #   std.error <dbl>, statistic <dbl>, df <dbl>, p.value <dbl>, conf.low <dbl>,
    ## #   conf.high <dbl>, sigma <dbl>, FDR <dbl>

Protein and contrast pairs that `lm_impute` could not estimate, if any:

``` r
fa$get_missing()
```

    ## # A tibble: 0 × 2
    ## # ℹ 2 variables: protein_Id <chr>, contrast <chr>

## Gene Set Enrichment Analysis (GSEA)

After computing contrasts, proteins can be ranked by log2 fold change or
t-statistic and subjected to gene set enrichment analysis. The typical
workflow is:

1.  Extract UniProt IDs from the contrast results using
    [`prolfqua::get_uniprot_id_from_fasta_header()`](https://fgcz.github.io/prolfqua/reference/get_uniprot_id_from_fasta_header.md)
2.  Map UniProt IDs to Entrez Gene IDs using `AnnotationDbi::mapIds()`
    with a species-specific annotation package
3.  Create a named ranked list (statistic values named by gene IDs)
4.  Run GSEA using `clusterProfiler::gseGO()` or
    `clusterProfiler::gseKEGG()`

For details see:

- [clusterProfiler
  book](https://yulab-smu.top/biomedical-knowledge-mining-book/) –
  comprehensive guide to enrichment analysis in R
- [clusterProfiler Bioconductor
  page](https://bioconductor.org/packages/clusterProfiler/)
- [AnnotationDbi](https://bioconductor.org/packages/AnnotationDbi/) –
  for mapping between gene identifier types

The `prolfqua` package is described in (Wolski et al. 2023).

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
    ## [1] prolfqua_1.9.10  dplyr_1.2.1      BiocStyle_2.40.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1      viridisLite_0.4.3     farver_2.1.2         
    ##  [4] S7_1.0.0              fastmap_1.2.0         promises_1.5.0       
    ##  [7] digest_0.6.39         mime_0.13             lifecycle_1.0.5      
    ## [10] cluster_2.1.8.3       statmod_1.5.2         magrittr_2.0.5       
    ## [13] compiler_4.6.2        progress_1.2.3        rlang_1.3.0          
    ## [16] sass_0.5.0            tools_4.6.2           utf8_1.2.6           
    ## [19] yaml_2.3.12           data.table_1.18.6.1   knitr_1.52           
    ## [22] prettyunits_1.2.0     labeling_0.4.3        htmlwidgets_1.6.4    
    ## [25] RColorBrewer_1.1-3    KernSmooth_2.23-27    withr_3.0.3          
    ## [28] purrr_1.2.2           BiocGenerics_0.58.1   desc_1.4.3           
    ## [31] grid_4.6.2            stats4_4.6.2          xtable_1.8-8         
    ## [34] colorspace_2.1-4      ggplot2_4.0.4         scales_1.4.0         
    ## [37] iterators_1.0.14      MASS_7.3-66           cli_3.6.6            
    ## [40] rmarkdown_2.32        crayon_1.5.3          ragg_1.5.2           
    ## [43] generics_0.1.4        otel_0.3.0            httr_1.4.9           
    ## [46] rjson_0.2.23          cachem_1.1.0          stringr_1.6.0        
    ## [49] parallel_4.6.2        BiocManager_1.30.27   matrixStats_1.5.0    
    ## [52] vctrs_0.7.3           jsonlite_2.0.0        bookdown_0.48        
    ## [55] hms_1.1.4             IRanges_2.46.0        GetoptLong_1.1.1     
    ## [58] S4Vectors_0.50.3      ggrepel_0.9.8         clue_0.3-68          
    ## [61] systemfonts_1.3.4     crosstalk_1.2.2       foreach_1.5.2        
    ## [64] limma_3.68.5          plotly_4.12.1         tidyr_1.3.2          
    ## [67] jquerylib_0.1.4       glue_1.8.1            pkgdown_2.2.1        
    ## [70] codetools_0.2-20      stringi_1.8.9         gtable_0.3.6         
    ## [73] shape_1.4.6.1         later_1.4.8           ComplexHeatmap_2.28.0
    ## [76] tibble_3.3.1          pillar_1.11.1         htmltools_0.5.9      
    ## [79] circlize_0.4.18       R6_2.6.1              textshaping_1.1.0    
    ## [82] doParallel_1.0.17     evaluate_1.0.5        shiny_1.15.0         
    ## [85] png_0.1-9             httpuv_1.6.17         bslib_0.12.0         
    ## [88] Rcpp_1.1.2            gridExtra_2.3.1       xfun_0.62            
    ## [91] fs_2.1.0              pkgconfig_2.0.3       GlobalOptions_0.1.4

## References

Wolski, Witold E., Paolo Nanni, Jonas Grossmann, Maria d’Errico, Ralph
Schlapbach, and Christian Panse. 2023. “prolfqua: A Comprehensive
R-Package for Proteomics Differential Expression Analysis.” *Journal of
Proteome Research* 22 (4): 1092–1104.
<https://doi.org/10.1021/acs.jproteome.2c00441>.
