# Transforming and calibrating LFQ data

## The active response

`LFQData` holds the measurements and their configuration. Its
`response()` method gives the name of the intensity column used by
downstream operations. `get_Transformer()` creates a transformer with
its own copy of the data. Transformation methods return the transformer
for chaining, and `$lfq` retrieves the result.

This example uses protein-level data from three samples. The two
standards, `StdA` and `StdB`, each have a constant amount across
samples. Their measured signals are doubled in `S2` and halved in `S3`
relative to `S1`.

``` r
library(prolfqua)

samples <- c("S1", "S2", "S3")
quant <- data.frame(
  protein_Id = rep(c("StdA", "StdB", "Target"), each = 3),
  sampleName = rep(samples, 3),
  raw.file = rep(paste0(samples, ".raw"), 3),
  condition = rep(c("Control", "Control", "Treatment"), 3),
  intensity = c(64, 128, 32, 256, 512, 128, 512, 1024, 512),
  nr_children = 1
)

cfg <- AnalysisConfiguration$new()
cfg$file_name <- "raw.file"
cfg$hierarchy <- list(protein_Id = "protein_Id")
cfg$factors <- list(condition = "condition")
cfg$set_response("intensity")
lfq <- LFQData$new(quant, cfg, setup = TRUE)

logged <- lfq$get_Transformer()$log2()$lfq
logged$response()
#> [1] "log2_intensity"
head(logged$data_long()[c("protein_Id", "sampleName", lfq$response(), logged$response())])
#> # A tibble: 6 × 4
#>   protein_Id sampleName intensity log2_intensity
#>   <chr>      <chr>          <dbl>          <dbl>
#> 1 StdA       S1                64              6
#> 2 StdB       S1               256              8
#> 3 Target     S1               512              9
#> 4 StdA       S2               128              7
#> 5 StdB       S2               512              9
#> 6 Target     S2              1024             10
```

[`log2()`](https://rdrr.io/r/base/Log.html) expects positive measured
intensities and carries missing intensities through as `NA`. Convert
zeros that encode non-detection to `NA` before logging. The
`is_transformed()` flag records transformation status. Track the scale
separately when using an operation that requires a particular scale.

| Method | Operation on the active response |
|----|----|
| [`log2()`](https://rdrr.io/r/base/Log.html) | Take the base-2 logarithm of each intensity. |
| `intensity_array(.func)` | Apply a function such as `asinh` to the response vector. |
| `intensity_matrix(.func)` | Apply a matrix function to features in rows and samples in columns. For example, `intensity_matrix(vsn::justvsn)` fits a variance-stabilizing transformation to raw intensities. |
| `robscale()` | Centre each sample by its median and scale by its relative median absolute deviation. The default restores a common level based on the original medians. |
| `robscale_subset(lfqsubset)` | Estimate the medians and median absolute deviations from a subset of features and apply them to the full data. |

## Calibrate against internal standards

Internal standards are reference proteins measured in the same samples
as the targets. Select proteins expected to have stable amounts across
samples. `center_to_reference()` subtracts the median log2 intensity of
the reference proteins from every feature in each sample. The median of
the reference proteins is then zero in each sample with a usable
reference. When all reference intensities for a sample are missing, its
centred intensities are missing.

``` r
reference_ids <- data.frame(protein_Id = c("StdA", "StdB"))
reference <- logged$get_subset(reference_ids)
calibrated_internal <- logged$get_Transformer()$center_to_reference(reference)$lfq
```

In prolfquapp, give the measured protein IDs in the DEA YAML under
`processing_options$internal`. For a protein-level analysis,
`transform: none` means a log2 transform without additional global
normalization. With `internal` configured, prolfquapp computes a median
from the finite log2 reference intensities per sample and calls
`subtract_sample_offsets()`. If a sample has no finite reference
intensity, it stops and lists the affected sample IDs.

``` yaml
processing_options:
  transform: none
  internal: [StdA, StdB]
```

The VSN and robscale paths centre on their respective transformed
scales. For nested analyses, prolfquapp uses the internal standards for
the protein abundances shown in the report. The peptide model uses the
transformed peptide data before that protein-level centring.

## Apply offsets from external calibrants

Use `subtract_sample_offsets()` when the correction comes from
measurements outside the LFQ table. Supply one finite log2 offset per
sample, named by the configured sample ID. For example, suppose an
external calibrant has a reference signal of 128, and you have matched
one calibrant measurement to each sample. Divide each measured signal by
128 to obtain the fold change in response, then take its log2.

``` r
external <- data.frame(
  sampleName = c("S3", "S1", "S2"),
  measured_intensity = c(64, 128, 256)
)
offsets <- setNames(log2(external$measured_intensity / 128), external$sampleName)
offsets
#> S3 S1 S2 
#> -1  0  1

calibrated_external <- logged$get_Transformer()$subtract_sample_offsets(offsets)$lfq
```

The order of the calibrant rows can differ from the order in the LFQ
table. The method matches by sample ID and rejects missing, duplicate or
unknown IDs. The caller is responsible for establishing that the
external measurement tracks a multiplicative effect shared by the
features in that sample.

If the log2 intensity is $`y`$ and the offset is $`d`$, the calibrated
intensity is $`y - d`$. On the original scale this divides the
measurement by $`2^d`$. All features in a sample receive the same shift,
which preserves within-sample log2 differences. Between-sample
differences change by the difference between their offsets.

## Compare the two calibrations

``` r
target_values <- function(lfq) {
  rows <- lfq$data_long()
  target <- rows[rows$protein_Id == "Target", ]
  setNames(target[[lfq$response()]], target$sampleName)[samples]
}

knitr::kable(
  data.frame(
    Sample = samples,
    `Original log2` = unname(target_values(logged)),
    `Internal standards` = unname(target_values(calibrated_internal)),
    `External calibrant` = unname(target_values(calibrated_external)),
    check.names = FALSE
  ),
  caption = "Target log2 intensity before and after calibration."
)
```

| Sample | Original log2 | Internal standards | External calibrant |
|:-------|--------------:|-------------------:|-------------------:|
| S1     |             9 |                  2 |                  9 |
| S2     |            10 |                  2 |                  9 |
| S3     |             9 |                  3 |                 10 |

Target log2 intensity before and after calibration.

Both calibrations recover a target that is equal in `S1` and `S2` and
twice as abundant in `S3`. The absolute log2 levels differ because
internal centring puts the reference median at zero, while the external
offsets are anchored to the calibrant’s reference signal of 128. The
example treats them as two estimates of the same sample effect, so it
applies each to a separate copy of the logged data.

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
#> [1] prolfqua_1.9.5   BiocStyle_2.40.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] limma_3.68.5        gtable_0.3.6        jsonlite_2.0.0     
#>  [4] dplyr_1.2.1         compiler_4.6.2      BiocManager_1.30.27
#>  [7] Rcpp_1.1.2          tidyselect_1.2.1    tidyr_1.3.2        
#> [10] jquerylib_0.1.4     systemfonts_1.3.4   scales_1.4.0       
#> [13] textshaping_1.1.0   statmod_1.5.2       yaml_2.3.12        
#> [16] fastmap_1.2.0       ggplot2_4.0.4       R6_2.6.1           
#> [19] generics_0.1.4      knitr_1.52          MASS_7.3-66        
#> [22] htmlwidgets_1.6.4   ggrepel_0.9.8       tibble_3.3.1       
#> [25] bookdown_0.48       desc_1.4.3          RColorBrewer_1.1-3 
#> [28] bslib_0.12.0        pillar_1.11.1       rlang_1.3.0        
#> [31] utf8_1.2.6          cachem_1.1.0        xfun_0.62          
#> [34] S7_1.0.0            fs_2.1.0            sass_0.5.0         
#> [37] otel_0.3.0          cli_3.6.6           withr_3.0.3        
#> [40] pkgdown_2.2.1       magrittr_2.0.5      digest_0.6.39      
#> [43] grid_4.6.2          lifecycle_1.0.5     vctrs_0.7.3        
#> [46] evaluate_1.0.5      glue_1.8.1          farver_2.1.2       
#> [49] ragg_1.5.2          purrr_1.2.2         rmarkdown_2.32     
#> [52] tools_4.6.2         pkgconfig_2.0.3     htmltools_0.5.9
```
