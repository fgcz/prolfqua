# normalize data by log2 and robust scaling

normalize data by log2 and robust scaling

## Usage

``` r
normalize_log2_robscale(pdata, config)
```

## Arguments

- pdata:

  data.frame

- config:

  AnalysisConfiguration

## Value

list with data.frame (data) and updated config (config)

## See also

Other preprocessing:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`apply_to_response_matrix()`](https://wolski.github.io/prolfqua/reference/apply_to_response_matrix.md),
[`filter_proteins_by_peptide_count()`](https://wolski.github.io/prolfqua/reference/filter_proteins_by_peptide_count.md),
[`get_robscales()`](https://wolski.github.io/prolfqua/reference/get_robscales.md),
[`robust_scale()`](https://wolski.github.io/prolfqua/reference/robust_scale.md),
[`scale_with_subset()`](https://wolski.github.io/prolfqua/reference/scale_with_subset.md),
[`scale_with_subset_by_factors()`](https://wolski.github.io/prolfqua/reference/scale_with_subset_by_factors.md)

## Examples

``` r

bb <- sim_lfq_data_peptide_config(Nprot = 100)
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
xx <- normalize_log2_robscale(bb$data, bb$config)
#> Column added : log2_abundance
#> Joining with `by = join_by(sampleName, isotopeLabel,
#> protein_Id, peptide_Id)`
xx$config$workIntensity
#> [1] "abundance"            "log2_abundance"       "transformedIntensity"
```
