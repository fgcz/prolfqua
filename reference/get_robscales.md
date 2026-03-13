# compute median and standard deviation for each sample

compute median and standard deviation for each sample

## Usage

``` r
get_robscales(data, config)
```

## See also

Other preprocessing:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`apply_to_response_matrix()`](https://wolski.github.io/prolfqua/reference/apply_to_response_matrix.md),
[`filter_proteins_by_peptide_count()`](https://wolski.github.io/prolfqua/reference/filter_proteins_by_peptide_count.md),
[`normalize_log2_robscale()`](https://wolski.github.io/prolfqua/reference/normalize_log2_robscale.md),
[`robust_scale()`](https://wolski.github.io/prolfqua/reference/robust_scale.md),
[`scale_with_subset()`](https://wolski.github.io/prolfqua/reference/scale_with_subset.md),
[`scale_with_subset_by_factors()`](https://wolski.github.io/prolfqua/reference/scale_with_subset_by_factors.md)

## Examples

``` r

bb <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
conf <- bb$config
sample_analysis <- bb$data
pepIntensityNormalized <- transform_work_intensity(sample_analysis, conf, log2)
#> Column added : log2_abundance
s1 <- get_robscales(pepIntensityNormalized, conf)

res <- scale_with_subset(pepIntensityNormalized, pepIntensityNormalized, conf)
#> Joining with `by = join_by(sampleName, isotopeLabel, protein_Id, peptide_Id)`
s2 <- get_robscales(res$data, conf)
abs(mean(s1$mads) - mean(s2$mads)) < 0.1
#> [1] TRUE

```
