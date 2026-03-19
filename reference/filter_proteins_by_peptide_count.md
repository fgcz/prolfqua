# Keep only those proteins with 2 IDENTIFIED peptides

Keep only those proteins with 2 IDENTIFIED peptides

## Usage

``` r
filter_proteins_by_peptide_count(pdata, config)
```

## Arguments

- pdata:

  data.frame

- config:

  AnalysisConfiguration

## Value

list with data.frame (data) and name of new column (name)

## See also

Other preprocessing:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`apply_to_response_matrix()`](https://wolski.github.io/prolfqua/reference/apply_to_response_matrix.md),
[`get_robscales()`](https://wolski.github.io/prolfqua/reference/get_robscales.md),
[`normalize_log2_robscale()`](https://wolski.github.io/prolfqua/reference/normalize_log2_robscale.md),
[`robust_scale()`](https://wolski.github.io/prolfqua/reference/robust_scale.md),
[`scale_with_subset()`](https://wolski.github.io/prolfqua/reference/scale_with_subset.md),
[`scale_with_subset_by_factors()`](https://wolski.github.io/prolfqua/reference/scale_with_subset_by_factors.md)

## Examples

``` r
istar <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
filterPep <- prolfqua::filter_proteins_by_peptide_count( istar$data ,  istar$config )
#> Column added : nr_peptide_Id_IN_protein_Id
x <- prolfqua::summarize_hierarchy(filterPep$data , istar$config)
stopifnot(x$peptide_Id_n >= istar$config$min_peptides_protein)
```
