# Create an interaction column from factors

Create an interaction column from factors

## Usage

``` r
make_interaction_column(data, columns, sep = ".")
```

## Arguments

- data:

  data frame

- columns:

  factor columns to combine

- sep:

  separator between the levels in the combined labels

## Value

`data` with an added factor column `interaction` combining the values of
`columns`, each prefixed with its column name and joined by `sep`

## See also

Other configuration:
[`AnalysisConfiguration`](https://fgcz.github.io/prolfqua/reference/AnalysisConfiguration.md),
[`ContrastConfiguration`](https://fgcz.github.io/prolfqua/reference/ContrastConfiguration.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`complete_cases()`](https://fgcz.github.io/prolfqua/reference/complete_cases.md),
[`list_to_contrast_configuration()`](https://fgcz.github.io/prolfqua/reference/list_to_contrast_configuration.md),
[`make_reduced_hierarchy_config()`](https://fgcz.github.io/prolfqua/reference/make_reduced_hierarchy_config.md),
[`r6_extract_values()`](https://fgcz.github.io/prolfqua/reference/r6_extract_values.md),
[`sample_subset()`](https://fgcz.github.io/prolfqua/reference/sample_subset.md),
[`separate_hierarchy()`](https://fgcz.github.io/prolfqua/reference/separate_hierarchy.md),
[`setup_analysis()`](https://fgcz.github.io/prolfqua/reference/setup_analysis.md),
[`table_factors()`](https://fgcz.github.io/prolfqua/reference/table_factors.md),
[`table_factors_size()`](https://fgcz.github.io/prolfqua/reference/table_factors_size.md)

## Examples

``` r
xx <- data.frame(A = c("a","a","a"), B = c("d","d","e"))
x <- make_interaction_column(xx, c("B","A"))
x <- make_interaction_column(xx, c("A"))
bb <- prolfqua::sim_lfq_data_protein_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 5 unconfigured columns: average_prot_abundance, sd, mean, N, Replicate.
#> setup_analysis: added 16 rows with NA response for feature and sample combinations absent from the input.
config <- bb$config
analysis <- bb$data

config$factor_depth <- 1
make_interaction_column(analysis,
   config$factor_keys_depth())
#> # A tibble: 120 × 8
#>    sample sampleName group_ protein_Id  abundance qValue nrPeptides interaction
#>    <chr>  <chr>      <chr>  <chr>           <dbl>  <dbl>      <dbl> <fct>      
#>  1 A_V1   A_V1       A      0EfVhX~0087      20.1      0          3 group_A    
#>  2 A_V1   A_V1       A      7cbcrd~5725      22.0      0          1 group_A    
#>  3 A_V1   A_V1       A      9VUkAq~4703      19.8      0          1 group_A    
#>  4 A_V1   A_V1       A      BEJI92~5282      21.2      0          2 group_A    
#>  5 A_V1   A_V1       A      CGzoYe~2147      29.4      0          1 group_A    
#>  6 A_V1   A_V1       A      DoWup2~5896      NA       NA         NA group_A    
#>  7 A_V1   A_V1       A      Fl4JiV~8625      20.1      0          4 group_A    
#>  8 A_V1   A_V1       A      HvIpHG~9079      21.7      0          2 group_A    
#>  9 A_V1   A_V1       A      JcKVfU~9653      34.5      0          7 group_A    
#> 10 A_V1   A_V1       A      SGIVBl~5782      23.9      0          6 group_A    
#> # ℹ 110 more rows
```
