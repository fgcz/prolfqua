# Apply function requiring a matrix to tidy table

Apply function requiring a matrix to tidy table

## Usage

``` r
apply_to_response_matrix(data, config, .func, .funcname = NULL)
```

## Arguments

- data:

  data.frame

- config:

  AnalysisConfiguration

- .func:

  function

- .funcname:

  name of function (used for creating new column)

## See also

Other preprocessing:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`filter_proteins_by_peptide_count()`](https://wolski.github.io/prolfqua/reference/filter_proteins_by_peptide_count.md),
[`get_robscales()`](https://wolski.github.io/prolfqua/reference/get_robscales.md),
[`normalize_log2_robscale()`](https://wolski.github.io/prolfqua/reference/normalize_log2_robscale.md),
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
data <- bb$data
conf <- bb$config
res <- apply_to_response_matrix(data, conf, .func = base::scale)
#> Joining with `by = join_by(sampleName, isotopeLabel,
#> protein_Id, peptide_Id)`

stopifnot("abundance_base..scale" %in% colnames(res))
stopifnot("abundance_base..scale" == conf$get_response())
conf <- bb$config$clone(deep=TRUE)
conf$workIntensity <- "abundance"
res <- apply_to_response_matrix(data, conf$clone(deep=TRUE), .func = robust_scale)
#> Joining with `by = join_by(sampleName, isotopeLabel,
#> protein_Id, peptide_Id)`

# Normalize data using the vsn method from bioconductor

if( require("vsn")){
 res <- apply_to_response_matrix(data, conf$clone(deep=TRUE), .func = vsn::justvsn)
}
#> Joining with `by = join_by(sampleName, isotopeLabel,
#> protein_Id, peptide_Id)`
```
