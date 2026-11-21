# Separates hierarchy columns into starting columns

Separates hierarchy columns into starting columns

## Usage

``` r
separate_hierarchy(data, config)
```

## Arguments

- data:

  data.frame

- config:

  AnalysisConfiguration

## Value

`data` with each hierarchy key at the current depth split at
`config$sep` into its component columns, which are added after the key
column

## See also

Other configuration:
[`AnalysisConfiguration`](https://fgcz.github.io/prolfqua/reference/AnalysisConfiguration.md),
[`ContrastConfiguration`](https://fgcz.github.io/prolfqua/reference/ContrastConfiguration.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`complete_cases()`](https://fgcz.github.io/prolfqua/reference/complete_cases.md),
[`list_to_contrast_configuration()`](https://fgcz.github.io/prolfqua/reference/list_to_contrast_configuration.md),
[`make_interaction_column()`](https://fgcz.github.io/prolfqua/reference/make_interaction_column.md),
[`make_reduced_hierarchy_config()`](https://fgcz.github.io/prolfqua/reference/make_reduced_hierarchy_config.md),
[`r6_extract_values()`](https://fgcz.github.io/prolfqua/reference/r6_extract_values.md),
[`sample_subset()`](https://fgcz.github.io/prolfqua/reference/sample_subset.md),
[`setup_analysis()`](https://fgcz.github.io/prolfqua/reference/setup_analysis.md),
[`table_factors()`](https://fgcz.github.io/prolfqua/reference/table_factors.md),
[`table_factors_size()`](https://fgcz.github.io/prolfqua/reference/table_factors_size.md)

## Examples

``` r

bb <- sim_lfq_data_protein_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 5 unconfigured columns: average_prot_abundance, sd, mean, N, Replicate.
#> setup_analysis: added 16 rows with NA response for feature and sample combinations absent from the input.
dt <- separate_hierarchy(bb$data, bb$config)
base::setdiff(colnames(dt) ,colnames(bb$data))
#> [1] "proteinID" "idtype2"  
stopifnot(ncol(dt) >= ncol(bb$data))
```
