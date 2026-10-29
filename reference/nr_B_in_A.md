# Count the children of each parent in the hierarchy

The number of distinct \`level_b\` entries, e.g. peptides, within each
\`level_a\` entry, e.g. protein.

## Usage

``` r
nr_B_in_A(
  pdata,
  level_a,
  level_b,
  name = paste0("nr_", paste(level_b, collapse = "_"), "_IN_", paste(level_a, collapse =
    "_"))
)
```

## Arguments

- pdata:

  data.frame

- level_a:

  parent key columns

- level_b:

  child key columns

- name:

  name of the count column

## Value

data.frame of the \`level_a\` keys and the count column

## See also

Other preprocessing:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`apply_to_response_matrix()`](https://fgcz.github.io/prolfqua/reference/apply_to_response_matrix.md),
[`filter_proteins_by_peptide_count()`](https://fgcz.github.io/prolfqua/reference/filter_proteins_by_peptide_count.md),
[`get_robscales()`](https://fgcz.github.io/prolfqua/reference/get_robscales.md),
[`robust_scale()`](https://fgcz.github.io/prolfqua/reference/robust_scale.md),
[`scale_with_subset()`](https://fgcz.github.io/prolfqua/reference/scale_with_subset.md)

## Examples

``` r
istar <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
head(prolfqua:::nr_B_in_A(istar$data, "protein_Id", "peptide_Id"))
#> # A tibble: 6 × 2
#>   protein_Id  nr_peptide_Id_IN_protein_Id
#>   <chr>                             <int>
#> 1 0EfVhX~0087                           3
#> 2 7cbcrd~5725                           1
#> 3 9VUkAq~4703                           1
#> 4 BEJI92~5282                           2
#> 5 CGzoYe~2147                           1
#> 6 DoWup2~5896                           1
```
