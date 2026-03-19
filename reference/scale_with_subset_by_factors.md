# Scale using a subset of the data, within factor levels (e.g. use for pulldown data)

This method reduces the variance within the group.

## Usage

``` r
scale_with_subset_by_factors(data, subset, config, preserveMean = FALSE)
```

## Arguments

- data:

  tibble with data

- subset:

  tibble with subset of the data which will be used to derive scaling
  parameters

- config:

  configuration

- preserveMean:

  default FALSE then set mean to 0

## See also

Other preprocessing:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`apply_to_response_matrix()`](https://wolski.github.io/prolfqua/reference/apply_to_response_matrix.md),
[`filter_proteins_by_peptide_count()`](https://wolski.github.io/prolfqua/reference/filter_proteins_by_peptide_count.md),
[`get_robscales()`](https://wolski.github.io/prolfqua/reference/get_robscales.md),
[`normalize_log2_robscale()`](https://wolski.github.io/prolfqua/reference/normalize_log2_robscale.md),
[`robust_scale()`](https://wolski.github.io/prolfqua/reference/robust_scale.md),
[`scale_with_subset()`](https://wolski.github.io/prolfqua/reference/scale_with_subset.md)

## Examples

``` r


bb <- sim_lfq_data_peptide_config(Nprot = 100)
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
conf <- bb$config$clone(deep=TRUE)
sample_analysis <- bb$data

res <- transform_work_intensity(sample_analysis, conf, log2)
#> Column added : log2_abundance
res <- scale_with_subset_by_factors(res, res, conf)
#> Joining with `by = join_by(sampleName, isotopeLabel,
#> protein_Id, peptide_Id)`
#> Joining with `by = join_by(sampleName, isotopeLabel,
#> protein_Id, peptide_Id)`
#> Joining with `by = join_by(sampleName, isotopeLabel,
#> protein_Id, peptide_Id)`

```
