# Summarize hierarchy counts

E.g compute number of peptides for each protein

## Usage

``` r
summarize_hierarchy(
  pdata,
  hierarchy_keys,
  isotope_label = "isotopeLabel",
  hierarchy = hierarchy_keys,
  factors = character()
)
```

## Arguments

- pdata:

  data.frame

- hierarchy_keys:

  character vector — all hierarchy column names

- isotope_label:

  character — isotope label column name

- hierarchy:

  character vector — hierarchy level to group by (default
  hierarchy_keys)

- factors:

  character vector — factor columns to include

## See also

Other summary:
[`HierarchyCountsSample`](https://wolski.github.io/prolfqua/reference/HierarchyCountsSample.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`hierarchy_counts()`](https://wolski.github.io/prolfqua/reference/hierarchy_counts.md),
[`hierarchy_counts_sample()`](https://wolski.github.io/prolfqua/reference/hierarchy_counts_sample.md)

## Examples

``` r
bb <- sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq <- LFQData$new(bb$data, bb$config)
summarize_hierarchy(lfq$data_long(), lfq$hierarchy_keys(), lfq$isotope_label())
#> # A tibble: 28 × 3
#> # Groups:   protein_Id [10]
#>    protein_Id  peptide_Id isotopeLabel_n
#>    <chr>       <chr>               <int>
#>  1 0EfVhX~0087 ITLb4x1q                1
#>  2 0EfVhX~0087 ahQLlQY7                1
#>  3 0EfVhX~0087 dJkdz7so                1
#>  4 7cbcrd~5725 D5dQ4nKk                1
#>  5 9VUkAq~4703 eIC06D7g                1
#>  6 BEJI92~5282 HBkZvdhT                1
#>  7 BEJI92~5282 qQ1GK8Un                1
#>  8 CGzoYe~2147 mjHSHhoe                1
#>  9 DoWup2~5896 KVUnZ6oZ                1
#> 10 Fl4JiV~8625 GsUIOl6Q                1
#> # ℹ 18 more rows
summarize_hierarchy(lfq$data_long(), lfq$hierarchy_keys(), lfq$isotope_label(),
 factors = character())
#> # A tibble: 28 × 3
#> # Groups:   protein_Id [10]
#>    protein_Id  peptide_Id isotopeLabel_n
#>    <chr>       <chr>               <int>
#>  1 0EfVhX~0087 ITLb4x1q                1
#>  2 0EfVhX~0087 ahQLlQY7                1
#>  3 0EfVhX~0087 dJkdz7so                1
#>  4 7cbcrd~5725 D5dQ4nKk                1
#>  5 9VUkAq~4703 eIC06D7g                1
#>  6 BEJI92~5282 HBkZvdhT                1
#>  7 BEJI92~5282 qQ1GK8Un                1
#>  8 CGzoYe~2147 mjHSHhoe                1
#>  9 DoWup2~5896 KVUnZ6oZ                1
#> 10 Fl4JiV~8625 GsUIOl6Q                1
#> # ℹ 18 more rows
summarize_hierarchy(lfq$data_long(), lfq$hierarchy_keys(), lfq$isotope_label(),
 hierarchy = lfq$relevant_hierarchy_keys())
#> # A tibble: 10 × 3
#>    protein_Id  isotopeLabel_n peptide_Id_n
#>    <chr>                <int>        <int>
#>  1 0EfVhX~0087              1            3
#>  2 7cbcrd~5725              1            1
#>  3 9VUkAq~4703              1            1
#>  4 BEJI92~5282              1            2
#>  5 CGzoYe~2147              1            1
#>  6 DoWup2~5896              1            1
#>  7 Fl4JiV~8625              1            4
#>  8 HvIpHG~9079              1            2
#>  9 JcKVfU~9653              1            7
#> 10 SGIVBl~5782              1            6
```
