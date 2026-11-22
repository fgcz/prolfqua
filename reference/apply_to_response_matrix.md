# Apply function requiring a matrix to tidy table

Apply function requiring a matrix to tidy table

## Usage

``` r
apply_to_response_matrix(lfqdata, .func, .funcname = NULL)
```

## Arguments

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

- .func:

  function taking and returning a matrix

- .funcname:

  name of function (used for creating new column)

## Value

list with `data` (data.frame) and `colname` (new column name)

## See also

Other preprocessing:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`filter_proteins_by_peptide_count()`](https://fgcz.github.io/prolfqua/reference/filter_proteins_by_peptide_count.md),
[`get_robscales()`](https://fgcz.github.io/prolfqua/reference/get_robscales.md),
[`nr_B_in_A()`](https://fgcz.github.io/prolfqua/reference/nr_B_in_A.md),
[`robust_scale()`](https://fgcz.github.io/prolfqua/reference/robust_scale.md),
[`scale_with_subset()`](https://fgcz.github.io/prolfqua/reference/scale_with_subset.md)

## Examples

``` r
bb <- sim_lfq_data_peptide_config(n_proteins = 100)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 528 rows with NA response for feature and sample combinations absent from the input.
lfqdata <- LFQData$new(bb$data, bb$config)
res <- apply_to_response_matrix(lfqdata, .func = base::scale)
stopifnot("abundance_base..scale" %in% colnames(res$data))

res <- apply_to_response_matrix(lfqdata, .func = robust_scale)

if (requireNamespace("vsn", quietly = TRUE)) {
  res <- apply_to_response_matrix(lfqdata, .func = vsn::justvsn)
}
```
