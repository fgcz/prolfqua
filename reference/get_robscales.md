# compute median and standard deviation for each sample

compute median and standard deviation for each sample

## Usage

``` r
get_robscales(lfqdata)
```

## See also

Other preprocessing:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`apply_to_response_matrix()`](https://wolski.github.io/prolfqua/reference/apply_to_response_matrix.md),
[`filter_proteins_by_peptide_count()`](https://wolski.github.io/prolfqua/reference/filter_proteins_by_peptide_count.md),
[`robust_scale()`](https://wolski.github.io/prolfqua/reference/robust_scale.md),
[`scale_with_subset()`](https://wolski.github.io/prolfqua/reference/scale_with_subset.md)

## Examples

``` r

bb <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfqdata <- prolfqua::LFQData$new(bb$data, bb$config)
lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#> Column added : log2_abundance
s1 <- get_robscales(lfqdata)

res <- scale_with_subset(lfqdata, lfqdata)
#> Joining with `by = join_by(sampleName, isotopeLabel, protein_Id, peptide_Id)`
lfqres <- prolfqua::LFQData$new(res$data, lfqdata$get_config()$clone(deep = TRUE))
s2 <- get_robscales(lfqres)
abs(mean(s1$mads) - mean(s2$mads)) < 0.1
#> [1] TRUE

```
