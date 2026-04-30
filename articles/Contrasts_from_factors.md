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
[`group_label()`](https://wolski.github.io/prolfqua/reference/group_label.md):

``` r
library(prolfqua)

group_label("MI", "T0")
```

    ## [1] "G_MI_T0"

``` r
group_label("MINOCA", "T300")
```

    ## [1] "G_MINOCA_T300"

This means that before fitting a model, the data must contain a grouping
column with levels in this format. The `annotation_add_contrasts`
function creates such a column automatically using
[`tidyr::unite`](https://tidyr.tidyverse.org/reference/unite.html).

## Building contrasts step by step

Consider a two-factor design with disease type (MI, MINOCA) and time
point (T0, T150, T300):

``` r
primary_levels <- c("MI", "MINOCA")
secondary_levels <- c("T0", "T150", "T300")
```

### Main effect contrasts

Main effects average across all levels of the secondary factor. For
example, the main effect of MINOCA vs MI is the average difference
across all time points:

``` r
me <- main_effect_contrasts(primary_levels, secondary_levels)
data.frame(ContrastName = names(me), Contrast = unlist(me))
```

    ##              ContrastName
    ## MINOCA_vs_MI MINOCA_vs_MI
    ##                                                                                               Contrast
    ## MINOCA_vs_MI ( (G_MINOCA_T0 + G_MINOCA_T150 + G_MINOCA_T300)/3 - (G_MI_T0 + G_MI_T150 + G_MI_T300)/3 )

Swapping the roles of primary and secondary gives main effects for time
points averaged across disease types:

``` r
me2 <- main_effect_contrasts(secondary_levels, primary_levels)
data.frame(ContrastName = names(me2), Contrast = unlist(me2))
```

    ##              ContrastName
    ## T150_vs_T0     T150_vs_T0
    ## T300_vs_T0     T300_vs_T0
    ## T300_vs_T150 T300_vs_T150
    ##                                                                       Contrast
    ## T150_vs_T0       ( (G_T150_MI + G_T150_MINOCA)/2 - (G_T0_MI + G_T0_MINOCA)/2 )
    ## T300_vs_T0       ( (G_T300_MI + G_T300_MINOCA)/2 - (G_T0_MI + G_T0_MINOCA)/2 )
    ## T300_vs_T150 ( (G_T300_MI + G_T300_MINOCA)/2 - (G_T150_MI + G_T150_MINOCA)/2 )

### Level-specific contrasts

These compare primary factor levels at each individual level of the
secondary factor:

``` r
ls <- level_specific_contrasts(primary_levels, secondary_levels)
data.frame(ContrastName = names(ls), Contrast = unlist(ls))
```

    ##                              ContrastName                  Contrast
    ## MINOCA_vs_MI_at_T0     MINOCA_vs_MI_at_T0     G_MINOCA_T0 - G_MI_T0
    ## MINOCA_vs_MI_at_T150 MINOCA_vs_MI_at_T150 G_MINOCA_T150 - G_MI_T150
    ## MINOCA_vs_MI_at_T300 MINOCA_vs_MI_at_T300 G_MINOCA_T300 - G_MI_T300

### Interaction contrasts

Interaction contrasts test whether the difference between primary levels
changes across secondary levels (difference of differences):

``` r
ic <- interaction_contrasts(primary_levels, secondary_levels)
data.frame(ContrastName = names(ic), Contrast = unlist(ic))
```

    ##                                                                      ContrastName
    ## interaction_MINOCA_vs_MI_at_T150_vs_T0     interaction_MINOCA_vs_MI_at_T150_vs_T0
    ## interaction_MINOCA_vs_MI_at_T300_vs_T0     interaction_MINOCA_vs_MI_at_T300_vs_T0
    ## interaction_MINOCA_vs_MI_at_T300_vs_T150 interaction_MINOCA_vs_MI_at_T300_vs_T150
    ##                                                                                           Contrast
    ## interaction_MINOCA_vs_MI_at_T150_vs_T0       (G_MINOCA_T150 - G_MI_T150) - (G_MINOCA_T0 - G_MI_T0)
    ## interaction_MINOCA_vs_MI_at_T300_vs_T0       (G_MINOCA_T300 - G_MI_T300) - (G_MINOCA_T0 - G_MI_T0)
    ## interaction_MINOCA_vs_MI_at_T300_vs_T150 (G_MINOCA_T300 - G_MI_T300) - (G_MINOCA_T150 - G_MI_T150)

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

| ContrastName                             | Contrast                                                                                  |
|:-----------------------------------------|:------------------------------------------------------------------------------------------|
| MINOCA_vs_MI                             | ( (G_MINOCA_T0 + G_MINOCA_T150 + G_MINOCA_T300)/3 - (G_MI_T0 + G_MI_T150 + G_MI_T300)/3 ) |
| MINOCA_vs_MI_at_T0                       | G_MINOCA_T0 - G_MI_T0                                                                     |
| MINOCA_vs_MI_at_T150                     | G_MINOCA_T150 - G_MI_T150                                                                 |
| MINOCA_vs_MI_at_T300                     | G_MINOCA_T300 - G_MI_T300                                                                 |
| interaction_MINOCA_vs_MI_at_T150_vs_T0   | (G_MINOCA_T150 - G_MI_T150) - (G_MINOCA_T0 - G_MI_T0)                                     |
| interaction_MINOCA_vs_MI_at_T300_vs_T0   | (G_MINOCA_T300 - G_MI_T300) - (G_MINOCA_T0 - G_MI_T0)                                     |
| interaction_MINOCA_vs_MI_at_T300_vs_T150 | (G_MINOCA_T300 - G_MI_T300) - (G_MINOCA_T150 - G_MI_T150)                                 |

To exclude interaction contrasts:

``` r
no_int <- generate_contrasts(primary_levels, secondary_levels, interactions = FALSE)
knitr::kable(no_int, row.names = FALSE)
```

| ContrastName         | Contrast                                                                                  |
|:---------------------|:------------------------------------------------------------------------------------------|
| MINOCA_vs_MI         | ( (G_MINOCA_T0 + G_MINOCA_T150 + G_MINOCA_T300)/3 - (G_MI_T0 + G_MI_T150 + G_MI_T300)/3 ) |
| MINOCA_vs_MI_at_T0   | G_MINOCA_T0 - G_MI_T0                                                                     |
| MINOCA_vs_MI_at_T150 | G_MINOCA_T150 - G_MI_T150                                                                 |
| MINOCA_vs_MI_at_T300 | G_MINOCA_T300 - G_MI_T300                                                                 |

## Working with annotation tables

In a typical prolfquapp workflow, you start with a sample annotation
table that has columns for the two factors. `annotation_add_contrasts`
creates a united Group column, generates all contrasts, and binds them
alongside the annotation:

``` r
# x5463yzwer453bbb is a bundled example annotation table
# with factor_A (MI / MINOCA) and factor_B (T0 / T150 / T300)
data("x5463yzwer453bbb", package = "prolfqua")
head(x5463yzwer453bbb[, c("Name", "Group", "factor_A", "factor_B")])
```

    ## # A tibble: 6 × 4
    ##   Name         Group       factor_A factor_B
    ##   <chr>        <chr>       <chr>    <chr>   
    ## 1 MI_150_6     MI_T150     MI       T150    
    ## 2 MINOCA_0_3   MINOCA_T0   MINOCA   T0      
    ## 3 MINOCA_300_5 MINOCA_T300 MINOCA   T300    
    ## 4 MI_150_5     MI_T150     MI       T150    
    ## 5 MI_150_1     MI_T150     MI       T150    
    ## 6 MI_0_3       MI_T0       MI       T0

``` r
result <- annotation_add_contrasts(
  x5463yzwer453bbb,
  primary_col = "factor_A",
  secondary_col = "factor_B",
  prefix = "primary"
)

# The annotation with united Group and contrast columns
knitr::kable(head(result$annot[, c("Name", "Group", "ContrastName", "Contrast")], 10),
             row.names = FALSE)
```

| Name         | Group       | ContrastName                             | Contrast                                                                                  |
|:-------------|:------------|:-----------------------------------------|:------------------------------------------------------------------------------------------|
| MI_150_6     | MI_T150     | MINOCA_vs_MI                             | ( (G_MINOCA_T0 + G_MINOCA_T150 + G_MINOCA_T300)/3 - (G_MI_T0 + G_MI_T150 + G_MI_T300)/3 ) |
| MINOCA_0_3   | MINOCA_T0   | MINOCA_vs_MI_at_T0                       | G_MINOCA_T0 - G_MI_T0                                                                     |
| MINOCA_300_5 | MINOCA_T300 | MINOCA_vs_MI_at_T150                     | G_MINOCA_T150 - G_MI_T150                                                                 |
| MI_150_5     | MI_T150     | MINOCA_vs_MI_at_T300                     | G_MINOCA_T300 - G_MI_T300                                                                 |
| MI_150_1     | MI_T150     | interaction_MINOCA_vs_MI_at_T150_vs_T0   | (G_MINOCA_T150 - G_MI_T150) - (G_MINOCA_T0 - G_MI_T0)                                     |
| MI_0_3       | MI_T0       | interaction_MINOCA_vs_MI_at_T300_vs_T0   | (G_MINOCA_T300 - G_MI_T300) - (G_MINOCA_T0 - G_MI_T0)                                     |
| MINOCA_150_1 | MINOCA_T150 | interaction_MINOCA_vs_MI_at_T300_vs_T150 | (G_MINOCA_T300 - G_MI_T300) - (G_MINOCA_T150 - G_MI_T150)                                 |
| MI_300_2     | MI_T300     | NA                                       | NA                                                                                        |
| MINOCA_150_5 | MINOCA_T150 | NA                                       | NA                                                                                        |
| MI_300_3     | MI_T300     | NA                                       | NA                                                                                        |

``` r
# Suggested output file name
result$name
```

    ## [1] "DEA_primary_dataset.csv"

Swapping primary and secondary factors gives contrasts from the other
perspective:

``` r
result2 <- annotation_add_contrasts(
  x5463yzwer453bbb,
  primary_col = "factor_B",
  secondary_col = "factor_A",
  prefix = "secondary"
)
knitr::kable(head(result2$annot[, c("Name", "Group", "ContrastName", "Contrast")], 10),
             row.names = FALSE)
```

| Name         | Group       | ContrastName                           | Contrast                                                          |
|:-------------|:------------|:---------------------------------------|:------------------------------------------------------------------|
| MI_150_6     | T150_MI     | T150_vs_T0                             | ( (G_T150_MI + G_T150_MINOCA)/2 - (G_T0_MI + G_T0_MINOCA)/2 )     |
| MINOCA_0_3   | T0_MINOCA   | T300_vs_T0                             | ( (G_T300_MI + G_T300_MINOCA)/2 - (G_T0_MI + G_T0_MINOCA)/2 )     |
| MINOCA_300_5 | T300_MINOCA | T300_vs_T150                           | ( (G_T300_MI + G_T300_MINOCA)/2 - (G_T150_MI + G_T150_MINOCA)/2 ) |
| MI_150_5     | T150_MI     | T150_vs_T0_at_MI                       | G_T150_MI - G_T0_MI                                               |
| MI_150_1     | T150_MI     | T150_vs_T0_at_MINOCA                   | G_T150_MINOCA - G_T0_MINOCA                                       |
| MI_0_3       | T0_MI       | T300_vs_T0_at_MI                       | G_T300_MI - G_T0_MI                                               |
| MINOCA_150_1 | T150_MINOCA | T300_vs_T0_at_MINOCA                   | G_T300_MINOCA - G_T0_MINOCA                                       |
| MI_300_2     | T300_MI     | T300_vs_T150_at_MI                     | G_T300_MI - G_T150_MI                                             |
| MINOCA_150_5 | T150_MINOCA | T300_vs_T150_at_MINOCA                 | G_T300_MINOCA - G_T150_MINOCA                                     |
| MI_300_3     | T300_MI     | interaction_T150_vs_T0_at_MINOCA_vs_MI | (G_T150_MINOCA - G_T0_MINOCA) - (G_T150_MI - G_T0_MI)             |

## Session Info

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
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
    ## [1] prolfqua_1.6.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1       viridisLite_0.4.3      dplyr_1.2.1           
    ##  [4] farver_2.1.2           S7_0.2.2               fastmap_1.2.0         
    ##  [7] lazyeval_0.2.3         digest_0.6.39          rpart_4.1.24          
    ## [10] lifecycle_1.0.5        survival_3.8-3         statmod_1.5.1         
    ## [13] magrittr_2.0.5         compiler_4.5.2         rlang_1.2.0           
    ## [16] sass_0.4.10            tools_4.5.2            utf8_1.2.6            
    ## [19] yaml_2.3.12            data.table_1.18.2.1    knitr_1.51            
    ## [22] htmlwidgets_1.6.4      plyr_1.8.9             RColorBrewer_1.1-3    
    ## [25] withr_3.0.2            purrr_1.2.2            desc_1.4.3            
    ## [28] nnet_7.3-20            grid_4.5.2             jomo_2.7-6            
    ## [31] mice_3.19.0            ggplot2_4.0.3          scales_1.4.0          
    ## [34] iterators_1.0.14       MASS_7.3-65            cli_3.6.6             
    ## [37] UpSetR_1.4.0           rmarkdown_2.31         ragg_1.5.2            
    ## [40] reformulas_0.4.4       generics_0.1.4         otel_0.2.0            
    ## [43] httr_1.4.8             minqa_1.2.8            cachem_1.1.0          
    ## [46] operator.tools_1.6.3.1 splines_4.5.2          vctrs_0.7.3           
    ## [49] boot_1.3-32            glmnet_4.1-10          Matrix_1.7-4          
    ## [52] jsonlite_2.0.0         mitml_0.4-5            ggrepel_0.9.8         
    ## [55] systemfonts_1.3.2      foreach_1.5.2          limma_3.66.0          
    ## [58] plotly_4.12.0          tidyr_1.3.2            jquerylib_0.1.4       
    ## [61] glue_1.8.1             pkgdown_2.2.0          nloptr_2.2.1          
    ## [64] pan_1.9                codetools_0.2-20       shape_1.4.6.1         
    ## [67] gtable_0.3.6           lme4_2.0-1             tibble_3.3.1          
    ## [70] pillar_1.11.1          htmltools_0.5.9        R6_2.6.1              
    ## [73] textshaping_1.0.5      Rdpack_2.6.6           formula.tools_1.7.1   
    ## [76] evaluate_1.0.5         lattice_0.22-7         rbibutils_2.4.1       
    ## [79] backports_1.5.1        pheatmap_1.0.13        broom_1.0.12          
    ## [82] bslib_0.10.0           Rcpp_1.1.1-1.1         gridExtra_2.3         
    ## [85] nlme_3.1-168           mgcv_1.9-3             logistf_1.26.1        
    ## [88] xfun_0.57              fs_2.1.0               forcats_1.0.1         
    ## [91] pkgconfig_2.0.3
