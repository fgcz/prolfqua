# Contrasts from factors

## Purpose

When working with factorial experimental designs (e.g. two conditions
crossed with three time points), specifying all the relevant contrasts
by hand is tedious and error-prone. The `generate_contrasts` family of
functions in prolfqua automates this by generating contrast
specifications from factor levels. These contrast strings can then be
passed directly to the `Contrasts` class for statistical testing.

This vignette demonstrates how to:

- Generate main effect, level-specific, and interaction contrasts for a
  two-factor design
- Use `annotation_add_contrasts` to produce a combined annotation and
  contrast table ready for analysis

## Group labelling convention

All contrast generation functions assume that the group levels in the
fitted model follow the naming convention `G_<primary>_<secondary>`,
which is produced by
[`group_label()`](https://fgcz.github.io/prolfqua/reference/group_label.md):

``` r
library(prolfqua)

group_label("WT", "T0")
```

    ## [1] "G_WT_T0"

``` r
group_label("KO", "T300")
```

    ## [1] "G_KO_T300"

This means that before fitting a model, the data must contain a grouping
column with levels in this format. The `annotation_add_contrasts`
function creates such a column automatically using
[`tidyr::unite`](https://tidyr.tidyverse.org/reference/unite.html).

## Building contrasts step by step

Consider a two-factor design with genotype (WT, KO) and time point (T0,
T150, T300):

``` r
primary_levels <- c("WT", "KO")
secondary_levels <- c("T0", "T150", "T300")
```

### Main effect contrasts

Main effects average across all levels of the secondary factor. For
example, the main effect of KO vs WT is the average difference across
all time points:

``` r
me <- main_effect_contrasts(primary_levels, secondary_levels)
data.frame(ContrastName = names(me), Contrast = unlist(me))
```

    ##          ContrastName
    ## KO_vs_WT     KO_vs_WT
    ##                                                                           Contrast
    ## KO_vs_WT (G_KO_T0 + G_KO_T150 + G_KO_T300)/3 - (G_WT_T0 + G_WT_T150 + G_WT_T300)/3

Swapping the roles of primary and secondary gives main effects for time
points averaged across genotypes:

``` r
me2 <- main_effect_contrasts(secondary_levels, primary_levels)
data.frame(ContrastName = names(me2), Contrast = unlist(me2))
```

    ##              ContrastName                                              Contrast
    ## T150_vs_T0     T150_vs_T0     (G_T150_WT + G_T150_KO)/2 - (G_T0_WT + G_T0_KO)/2
    ## T300_vs_T0     T300_vs_T0     (G_T300_WT + G_T300_KO)/2 - (G_T0_WT + G_T0_KO)/2
    ## T300_vs_T150 T300_vs_T150 (G_T300_WT + G_T300_KO)/2 - (G_T150_WT + G_T150_KO)/2

### Level-specific contrasts

These compare primary factor levels at each individual level of the
secondary factor:

``` r
ls <- level_specific_contrasts(primary_levels, secondary_levels)
data.frame(ContrastName = names(ls), Contrast = unlist(ls))
```

    ##                      ContrastName              Contrast
    ## KO_vs_WT_at_T0     KO_vs_WT_at_T0     G_KO_T0 - G_WT_T0
    ## KO_vs_WT_at_T150 KO_vs_WT_at_T150 G_KO_T150 - G_WT_T150
    ## KO_vs_WT_at_T300 KO_vs_WT_at_T300 G_KO_T300 - G_WT_T300

### Interaction contrasts

Interaction contrasts test whether the difference between primary levels
changes across secondary levels (difference of differences):

``` r
ic <- interaction_contrasts(primary_levels, secondary_levels)
data.frame(ContrastName = names(ic), Contrast = unlist(ic))
```

    ##                                                              ContrastName
    ## interaction_KO_vs_WT_at_T150_vs_T0     interaction_KO_vs_WT_at_T150_vs_T0
    ## interaction_KO_vs_WT_at_T300_vs_T0     interaction_KO_vs_WT_at_T300_vs_T0
    ## interaction_KO_vs_WT_at_T300_vs_T150 interaction_KO_vs_WT_at_T300_vs_T150
    ##                                                                               Contrast
    ## interaction_KO_vs_WT_at_T150_vs_T0       (G_KO_T150 - G_WT_T150) - (G_KO_T0 - G_WT_T0)
    ## interaction_KO_vs_WT_at_T300_vs_T0       (G_KO_T300 - G_WT_T300) - (G_KO_T0 - G_WT_T0)
    ## interaction_KO_vs_WT_at_T300_vs_T150 (G_KO_T300 - G_WT_T300) - (G_KO_T150 - G_WT_T150)

### Single-factor contrasts

For a one-factor design, `generate_contrasts_for_factor` generates all
pairwise comparisons:

``` r
group_levels <- c("CondA", "CondB", "CondC")
sf <- generate_contrasts_for_factor(group_levels)
data.frame(ContrastName = names(sf), Contrast = unlist(sf))
```

    ##                  ContrastName      Contrast
    ## CondB_vs_CondA CondB_vs_CondA CondB - CondA
    ## CondC_vs_CondA CondC_vs_CondA CondC - CondA
    ## CondC_vs_CondB CondC_vs_CondB CondC - CondB

## Generating all contrasts at once

`generate_contrasts` combines main effects, level-specific, and
interaction contrasts into a single data frame:

``` r
all_contrasts <- generate_contrasts(primary_levels, secondary_levels)
knitr::kable(all_contrasts, row.names = FALSE)
```

| ContrastName | Contrast |
|:---|:---|
| KO_vs_WT | (G_KO_T0 + G_KO_T150 + G_KO_T300)/3 - (G_WT_T0 + G_WT_T150 + G_WT_T300)/3 |
| KO_vs_WT_at_T0 | G_KO_T0 - G_WT_T0 |
| KO_vs_WT_at_T150 | G_KO_T150 - G_WT_T150 |
| KO_vs_WT_at_T300 | G_KO_T300 - G_WT_T300 |
| interaction_KO_vs_WT_at_T150_vs_T0 | (G_KO_T150 - G_WT_T150) - (G_KO_T0 - G_WT_T0) |
| interaction_KO_vs_WT_at_T300_vs_T0 | (G_KO_T300 - G_WT_T300) - (G_KO_T0 - G_WT_T0) |
| interaction_KO_vs_WT_at_T300_vs_T150 | (G_KO_T300 - G_WT_T300) - (G_KO_T150 - G_WT_T150) |

To exclude interaction contrasts:

``` r
no_int <- generate_contrasts(primary_levels, secondary_levels, interactions = FALSE)
knitr::kable(no_int, row.names = FALSE)
```

| ContrastName | Contrast |
|:---|:---|
| KO_vs_WT | (G_KO_T0 + G_KO_T150 + G_KO_T300)/3 - (G_WT_T0 + G_WT_T150 + G_WT_T300)/3 |
| KO_vs_WT_at_T0 | G_KO_T0 - G_WT_T0 |
| KO_vs_WT_at_T150 | G_KO_T150 - G_WT_T150 |
| KO_vs_WT_at_T300 | G_KO_T300 - G_WT_T300 |

## Working with annotation tables

In a typical prolfquapp workflow, you start with a sample annotation
table that has columns for the two factors. `annotation_add_contrasts`
creates a united Group column, generates all contrasts, and binds them
alongside the annotation:

``` r
# an annotation table of 36 samples: factor_A (genotype WT / KO),
# factor_B (time T0 / T150 / T300), 6 replicates each
design <- expand.grid(
  rep = 1:6, factor_B = c("T0", "T150", "T300"), factor_A = c("WT", "KO"),
  stringsAsFactors = FALSE
)
annotation <- data.frame(
  Name = paste(design$factor_A, design$factor_B, design$rep, sep = "_"),
  factor_A = design$factor_A,
  factor_B = design$factor_B
)
head(annotation)
```

    ##      Name factor_A factor_B
    ## 1 WT_T0_1       WT       T0
    ## 2 WT_T0_2       WT       T0
    ## 3 WT_T0_3       WT       T0
    ## 4 WT_T0_4       WT       T0
    ## 5 WT_T0_5       WT       T0
    ## 6 WT_T0_6       WT       T0

``` r
result <- annotation_add_contrasts(
  annotation,
  primary_col = "factor_A",
  secondary_col = "factor_B",
  prefix = "primary"
)

# The annotation with united Group and contrast columns
knitr::kable(head(result$annot[, c("Name", "Group", "ContrastName", "Contrast")], 10),
             row.names = FALSE)
```

| Name | Group | ContrastName | Contrast |
|:---|:---|:---|:---|
| WT_T0_1 | WT_T0 | WT_vs_KO | (G_WT_T0 + G_WT_T150 + G_WT_T300)/3 - (G_KO_T0 + G_KO_T150 + G_KO_T300)/3 |
| WT_T0_2 | WT_T0 | WT_vs_KO_at_T0 | G_WT_T0 - G_KO_T0 |
| WT_T0_3 | WT_T0 | WT_vs_KO_at_T150 | G_WT_T150 - G_KO_T150 |
| WT_T0_4 | WT_T0 | WT_vs_KO_at_T300 | G_WT_T300 - G_KO_T300 |
| WT_T0_5 | WT_T0 | interaction_WT_vs_KO_at_T150_vs_T0 | (G_WT_T150 - G_KO_T150) - (G_WT_T0 - G_KO_T0) |
| WT_T0_6 | WT_T0 | interaction_WT_vs_KO_at_T300_vs_T0 | (G_WT_T300 - G_KO_T300) - (G_WT_T0 - G_KO_T0) |
| WT_T150_1 | WT_T150 | interaction_WT_vs_KO_at_T300_vs_T150 | (G_WT_T300 - G_KO_T300) - (G_WT_T150 - G_KO_T150) |
| WT_T150_2 | WT_T150 | NA | NA |
| WT_T150_3 | WT_T150 | NA | NA |
| WT_T150_4 | WT_T150 | NA | NA |

``` r
# Suggested output file name
result$name
```

    ## [1] "DEA_primary_dataset.csv"

Swapping primary and secondary factors gives contrasts from the other
perspective:

``` r
result2 <- annotation_add_contrasts(
  annotation,
  primary_col = "factor_B",
  secondary_col = "factor_A",
  prefix = "secondary"
)
knitr::kable(head(result2$annot[, c("Name", "Group", "ContrastName", "Contrast")], 10),
             row.names = FALSE)
```

| Name | Group | ContrastName | Contrast |
|:---|:---|:---|:---|
| WT_T0_1 | T0_WT | T150_vs_T0 | (G_T150_KO + G_T150_WT)/2 - (G_T0_KO + G_T0_WT)/2 |
| WT_T0_2 | T0_WT | T300_vs_T0 | (G_T300_KO + G_T300_WT)/2 - (G_T0_KO + G_T0_WT)/2 |
| WT_T0_3 | T0_WT | T300_vs_T150 | (G_T300_KO + G_T300_WT)/2 - (G_T150_KO + G_T150_WT)/2 |
| WT_T0_4 | T0_WT | T150_vs_T0_at_KO | G_T150_KO - G_T0_KO |
| WT_T0_5 | T0_WT | T150_vs_T0_at_WT | G_T150_WT - G_T0_WT |
| WT_T0_6 | T0_WT | T300_vs_T0_at_KO | G_T300_KO - G_T0_KO |
| WT_T150_1 | T150_WT | T300_vs_T0_at_WT | G_T300_WT - G_T0_WT |
| WT_T150_2 | T150_WT | T300_vs_T150_at_KO | G_T300_KO - G_T150_KO |
| WT_T150_3 | T150_WT | T300_vs_T150_at_WT | G_T300_WT - G_T150_WT |
| WT_T150_4 | T150_WT | interaction_T150_vs_T0_at_WT_vs_KO | (G_T150_WT - G_T0_WT) - (G_T150_KO - G_T0_KO) |

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
    ## [1] prolfqua_1.9.5   BiocStyle_2.40.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] limma_3.68.5        gtable_0.3.6        jsonlite_2.0.0     
    ##  [4] dplyr_1.2.1         compiler_4.6.2      BiocManager_1.30.27
    ##  [7] Rcpp_1.1.2          tidyselect_1.2.1    tidyr_1.3.2        
    ## [10] jquerylib_0.1.4     systemfonts_1.3.4   scales_1.4.0       
    ## [13] textshaping_1.1.0   statmod_1.5.2       yaml_2.3.12        
    ## [16] fastmap_1.2.0       ggplot2_4.0.4       R6_2.6.1           
    ## [19] generics_0.1.4      knitr_1.52          MASS_7.3-66        
    ## [22] htmlwidgets_1.6.4   ggrepel_0.9.8       tibble_3.3.1       
    ## [25] bookdown_0.48       desc_1.4.3          RColorBrewer_1.1-3 
    ## [28] bslib_0.12.0        pillar_1.11.1       rlang_1.3.0        
    ## [31] cachem_1.1.0        xfun_0.62           S7_1.0.0           
    ## [34] fs_2.1.0            sass_0.5.0          otel_0.3.0         
    ## [37] cli_3.6.6           withr_3.0.3         pkgdown_2.2.1      
    ## [40] magrittr_2.0.5      digest_0.6.39       grid_4.6.2         
    ## [43] lifecycle_1.0.5     vctrs_0.7.3         evaluate_1.0.5     
    ## [46] glue_1.8.1          farver_2.1.2        ragg_1.5.2         
    ## [49] purrr_1.2.2         rmarkdown_2.32      tools_4.6.2        
    ## [52] pkgconfig_2.0.3     htmltools_0.5.9
