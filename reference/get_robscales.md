# compute median and standard deviation for each sample

compute median and standard deviation for each sample

## Usage

``` r
get_robscales(lfqdata)
```

## See also

Other preprocessing:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`apply_to_response_matrix()`](https://fgcz.github.io/prolfqua/reference/apply_to_response_matrix.md),
[`filter_proteins_by_peptide_count()`](https://fgcz.github.io/prolfqua/reference/filter_proteins_by_peptide_count.md),
[`nr_B_in_A()`](https://fgcz.github.io/prolfqua/reference/nr_B_in_A.md),
[`robust_scale()`](https://fgcz.github.io/prolfqua/reference/robust_scale.md),
[`scale_with_subset()`](https://fgcz.github.io/prolfqua/reference/scale_with_subset.md)

## Examples

``` r

bb <- prolfqua::sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
lfqdata <- prolfqua::LFQData$new(bb$data, bb$config)
lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#> Column added : log2_abundance
s1 <- get_robscales(lfqdata)

res <- scale_with_subset(lfqdata, lfqdata)
lfqres <- prolfqua::LFQData$new(res$data, lfqdata$get_config()$clone(deep = TRUE))
s2 <- get_robscales(lfqres)
abs(mean(s1$mads) - mean(s2$mads)) < 0.1
#> [1] TRUE

```
