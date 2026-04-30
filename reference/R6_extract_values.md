# Extract all value slots in an R6 object

Extract all value slots in an R6 object

## Usage

``` r
R6_extract_values(r6class)
```

## Arguments

- r6class:

  r6 class

## Value

The computed result.

## See also

Other configuration:
[`AnalysisConfiguration`](https://wolski.github.io/prolfqua/reference/AnalysisConfiguration.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`complete_cases()`](https://wolski.github.io/prolfqua/reference/complete_cases.md),
[`make_interaction_column()`](https://wolski.github.io/prolfqua/reference/make_interaction_column.md),
[`make_reduced_hierarchy_config()`](https://wolski.github.io/prolfqua/reference/make_reduced_hierarchy_config.md),
[`sample_subset()`](https://wolski.github.io/prolfqua/reference/sample_subset.md),
[`separate_hierarchy()`](https://wolski.github.io/prolfqua/reference/separate_hierarchy.md),
[`setup_analysis()`](https://wolski.github.io/prolfqua/reference/setup_analysis.md),
[`table_factors()`](https://wolski.github.io/prolfqua/reference/table_factors.md),
[`table_factors_size()`](https://wolski.github.io/prolfqua/reference/table_factors_size.md)

## Examples

``` r
bb <- sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
values <- R6_extract_values(bb$config)
names(values)
#>  [1] "min_peptides_protein"    "hierarchy_depth"        
#>  [3] "hierarchy"               "factor_depth"           
#>  [5] "factors"                 "bin_resp"               
#>  [7] "is_response_transformed" "work_intensity"         
#>  [9] "nr_children"             "opt_se"                 
#> [11] "opt_mz"                  "opt_rt"                 
#> [13] "ident_score"             "ident_q_value"          
#> [15] "isotope_label"           "sample_name"            
#> [17] "file_name"               "sep"                    
```
