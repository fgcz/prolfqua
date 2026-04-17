# Scale data using a subset of the data

this should reduce the overall variance.

## Usage

``` r
scale_with_subset(lfqdata, lfqsubset, preserve_mean = FALSE, get_scales = TRUE)
```

## Arguments

- lfqdata:

  LFQData object with full dataset

- lfqsubset:

  LFQData object with subset for computing scales

- preserve_mean:

  default FALSE - sets mean to zero

- get_scales:

  return a list of transformed data and the scaling parameters

## Value

list with data, scales, and colname

## See also

Other preprocessing:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`apply_to_response_matrix()`](https://wolski.github.io/prolfqua/reference/apply_to_response_matrix.md),
[`filter_proteins_by_peptide_count()`](https://wolski.github.io/prolfqua/reference/filter_proteins_by_peptide_count.md),
[`get_robscales()`](https://wolski.github.io/prolfqua/reference/get_robscales.md),
[`robust_scale()`](https://wolski.github.io/prolfqua/reference/robust_scale.md)

## Examples

``` r
bb <- sim_lfq_data_peptide_config(Nprot = 100)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfqdata <- LFQData$new(bb$data, bb$config)
lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#> Column added : log2_abundance
s1 <- get_robscales(lfqdata)
res <- scale_with_subset(lfqdata, lfqdata)
#> Joining with `by = join_by(sampleName, isotopeLabel, protein_Id, peptide_Id)`
cfg <- lfqdata$get_config()$clone(deep = TRUE)
cfg$set_response(res$colname)
lfqres <- LFQData$new(res$data, cfg)
s2 <- get_robscales(lfqres)
stopifnot(abs(mean(s1$mads) - mean(s2$mads)) < 1e-6)
```
