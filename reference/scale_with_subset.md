# Scale data using a subset of the data

this should reduce the overall variance.

## Usage

``` r
scale_with_subset(
  data,
  subset,
  config,
  preserveMean = FALSE,
  get_scales = TRUE
)
```

## Arguments

- data:

  the whole dataset

- subset:

  a subset of the dataset

- config:

  configuration

- preserveMean:

  default FALSE - sets mean to zero

- get_scales:

  return a list of transformed data and the scaling parameters

## See also

Other preprocessing:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`apply_to_response_matrix()`](https://wolski.github.io/prolfqua/reference/apply_to_response_matrix.md),
[`filter_proteins_by_peptide_count()`](https://wolski.github.io/prolfqua/reference/filter_proteins_by_peptide_count.md),
[`get_robscales()`](https://wolski.github.io/prolfqua/reference/get_robscales.md),
[`normalize_log2_robscale()`](https://wolski.github.io/prolfqua/reference/normalize_log2_robscale.md),
[`robust_scale()`](https://wolski.github.io/prolfqua/reference/robust_scale.md),
[`scale_with_subset_by_factors()`](https://wolski.github.io/prolfqua/reference/scale_with_subset_by_factors.md)

## Examples

``` r


bb <-sim_lfq_data_peptide_config(Nprot = 100)
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
conf <- bb$config$clone(deep=TRUE)
sample_analysis <- bb$data

res <- transform_work_intensity(sample_analysis, conf, log2)
#> Column added : log2_abundance
s1 <- get_robscales(res, conf)
res <- scale_with_subset(res, res, conf)
#> Warning: Expected 2 pieces. Additional pieces discarded in 4200 rows [1, 2, 3, 4, 5, 6,
#> 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].
#> Joining with `by = join_by(sampleName, protein_Id, peptide_Id)`
s2 <- get_robscales(res$data, conf)
stopifnot(abs(mean(s1$mads) - mean(s2$mads)) < 1e-6)
```
