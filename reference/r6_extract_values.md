# Extract all value slots in an R6 object

Turns an R6 object, such as an
[AnalysisConfiguration](https://fgcz.github.io/prolfqua/reference/AnalysisConfiguration.md),
into a plain list that can be written to YAML;
[`list_to_analysis_configuration()`](https://fgcz.github.io/prolfqua/reference/list_to_analysis_configuration.md)
reads it back.

## Usage

``` r
r6_extract_values(r6class)
```

## Arguments

- r6class:

  r6 class

## Value

a named list with the value fields of `r6class`; fields holding other R6
objects become nested lists, functions and environments are dropped.

## See also

Other configuration:
[`AnalysisConfiguration`](https://fgcz.github.io/prolfqua/reference/AnalysisConfiguration.md),
[`ContrastConfiguration`](https://fgcz.github.io/prolfqua/reference/ContrastConfiguration.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`list_to_contrast_configuration()`](https://fgcz.github.io/prolfqua/reference/list_to_contrast_configuration.md),
[`setup_analysis()`](https://fgcz.github.io/prolfqua/reference/setup_analysis.md)

## Examples

``` r
bb <- sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
values <- r6_extract_values(bb$config)
names(values)
#>  [1] "min_peptides_protein"    "hierarchy_depth"        
#>  [3] "hierarchy"               "factor_depth"           
#>  [5] "factors"                 "bin_resp"               
#>  [7] "is_response_transformed" "work_intensity"         
#>  [9] "nr_children"             "opt_se"                 
#> [11] "opt_mz"                  "opt_rt"                 
#> [13] "ident_score"             "ident_q_value"          
#> [15] "sample_name"             "file_name"              
#> [17] "sep"                    
```
