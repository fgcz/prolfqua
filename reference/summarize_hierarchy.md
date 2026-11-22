# Summarize hierarchy counts

E.g compute number of peptides for each protein: for each entry of
`hierarchy` (and of `factors`), the number of distinct entries of the
other hierarchy levels, in columns named `<level>_n`. Counts
experiment-wide, whether observed or not.

## Usage

``` r
summarize_hierarchy(
  pdata,
  hierarchy_keys,
  hierarchy = hierarchy_keys,
  factors = character()
)
```

## Arguments

- pdata:

  data.frame

- hierarchy_keys:

  character vector - all hierarchy column names

- hierarchy:

  character vector - hierarchy level to group by (default
  hierarchy_keys)

- factors:

  character vector - factor columns to include

## Value

a tibble with one row per combination of the `factors` and `hierarchy`
columns and, for each other hierarchy level, a column `<level>_n` with
the number of its distinct entries

## See also

Other summary:
[`HierarchyCountsSample`](https://fgcz.github.io/prolfqua/reference/HierarchyCountsSample.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`hierarchy_counts()`](https://fgcz.github.io/prolfqua/reference/hierarchy_counts.md)

Other hierarchy counting:
[`HierarchyCountsSample`](https://fgcz.github.io/prolfqua/reference/HierarchyCountsSample.md),
[`hierarchy_counts()`](https://fgcz.github.io/prolfqua/reference/hierarchy_counts.md),
[`nr_children_experiment()`](https://fgcz.github.io/prolfqua/reference/nr_children_experiment.md),
[`nr_features_experiment()`](https://fgcz.github.io/prolfqua/reference/nr_features_experiment.md),
[`nr_obs_sample()`](https://fgcz.github.io/prolfqua/reference/nr_obs_sample.md)

## Examples

``` r
bb <- sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
lfq <- LFQData$new(bb$data, bb$config)
summarize_hierarchy(lfq$data_long(), lfq$hierarchy_keys())
#> # A tibble: 28 × 2
#>    protein_Id  peptide_Id
#>    <chr>       <chr>     
#>  1 0EfVhX~0087 ITLb4x1q  
#>  2 0EfVhX~0087 ahQLlQY7  
#>  3 0EfVhX~0087 dJkdz7so  
#>  4 7cbcrd~5725 D5dQ4nKk  
#>  5 9VUkAq~4703 eIC06D7g  
#>  6 BEJI92~5282 HBkZvdhT  
#>  7 BEJI92~5282 qQ1GK8Un  
#>  8 CGzoYe~2147 mjHSHhoe  
#>  9 DoWup2~5896 KVUnZ6oZ  
#> 10 Fl4JiV~8625 GsUIOl6Q  
#> # ℹ 18 more rows
summarize_hierarchy(lfq$data_long(), lfq$hierarchy_keys(), factors = character())
#> # A tibble: 28 × 2
#>    protein_Id  peptide_Id
#>    <chr>       <chr>     
#>  1 0EfVhX~0087 ITLb4x1q  
#>  2 0EfVhX~0087 ahQLlQY7  
#>  3 0EfVhX~0087 dJkdz7so  
#>  4 7cbcrd~5725 D5dQ4nKk  
#>  5 9VUkAq~4703 eIC06D7g  
#>  6 BEJI92~5282 HBkZvdhT  
#>  7 BEJI92~5282 qQ1GK8Un  
#>  8 CGzoYe~2147 mjHSHhoe  
#>  9 DoWup2~5896 KVUnZ6oZ  
#> 10 Fl4JiV~8625 GsUIOl6Q  
#> # ℹ 18 more rows
summarize_hierarchy(lfq$data_long(), lfq$hierarchy_keys(),
 hierarchy = lfq$relevant_hierarchy_keys())
#> # A tibble: 10 × 2
#>    protein_Id  peptide_Id_n
#>    <chr>              <int>
#>  1 0EfVhX~0087            3
#>  2 7cbcrd~5725            1
#>  3 9VUkAq~4703            1
#>  4 BEJI92~5282            2
#>  5 CGzoYe~2147            1
#>  6 DoWup2~5896            1
#>  7 Fl4JiV~8625            4
#>  8 HvIpHG~9079            2
#>  9 JcKVfU~9653            7
#> 10 SGIVBl~5782            6
```
