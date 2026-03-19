# how many peptides per protein in each sample

how many peptides per protein in each sample

## Usage

``` r
nr_B_in_A_per_sample(data, config, nested = TRUE)
```

## See also

Other summary:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`hierarchy_counts()`](https://wolski.github.io/prolfqua/reference/hierarchy_counts.md),
[`hierarchy_counts_sample()`](https://wolski.github.io/prolfqua/reference/hierarchy_counts_sample.md),
[`summarize_hierarchy()`](https://wolski.github.io/prolfqua/reference/summarize_hierarchy.md)

## Examples

``` r

bb <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
nr_B_in_A_per_sample(bb$data, bb$config, nested =FALSE)
#> completing cases
#> # A tibble: 120 × 7
#>    sample sampleName protein_Id  group_ isotopeLabel nrPep present
#>    <chr>  <chr>      <chr>       <chr>  <chr>        <int>   <dbl>
#>  1 A_V1   A_V1       0EfVhX~0087 A      light            3       2
#>  2 A_V1   A_V1       7cbcrd~5725 A      light            1       1
#>  3 A_V1   A_V1       9VUkAq~4703 A      light            1       1
#>  4 A_V1   A_V1       BEJI92~5282 A      light            2       2
#>  5 A_V1   A_V1       CGzoYe~2147 A      light            1       1
#>  6 A_V1   A_V1       DoWup2~5896 A      light            1       1
#>  7 A_V1   A_V1       Fl4JiV~8625 A      light            4       3
#>  8 A_V1   A_V1       HvIpHG~9079 A      light            2       1
#>  9 A_V1   A_V1       JcKVfU~9653 A      light            7       7
#> 10 A_V1   A_V1       SGIVBl~5782 A      light            6       4
#> # ℹ 110 more rows
bb <- prolfqua::sim_lfq_data_protein_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
nr_B_in_A_per_sample(bb$data, bb$config, nested=FALSE)
#> Warning: here is no B in A
#> completing cases
#> # A tibble: 120 × 7
#>    sample sampleName protein_Id  group_ isotopeLabel nrPep present
#>    <chr>  <chr>      <chr>       <chr>  <chr>        <int>   <dbl>
#>  1 A_V1   A_V1       0EfVhX~0087 A      light            1       1
#>  2 A_V1   A_V1       7cbcrd~5725 A      light            1       1
#>  3 A_V1   A_V1       9VUkAq~4703 A      light            1       1
#>  4 A_V1   A_V1       BEJI92~5282 A      light            1       1
#>  5 A_V1   A_V1       CGzoYe~2147 A      light            1       1
#>  6 A_V1   A_V1       DoWup2~5896 A      light            1       0
#>  7 A_V1   A_V1       Fl4JiV~8625 A      light            1       1
#>  8 A_V1   A_V1       HvIpHG~9079 A      light            1       1
#>  9 A_V1   A_V1       JcKVfU~9653 A      light            1       1
#> 10 A_V1   A_V1       SGIVBl~5782 A      light            1       1
#> # ℹ 110 more rows
```
