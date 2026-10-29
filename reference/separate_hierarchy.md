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

## See also

Other configuration:
[`AnalysisConfiguration`](https://fgcz.github.io/prolfqua/reference/AnalysisConfiguration.md),
[`ContrastConfiguration`](https://fgcz.github.io/prolfqua/reference/ContrastConfiguration.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`R6_extract_values()`](https://fgcz.github.io/prolfqua/reference/R6_extract_values.md),
[`complete_cases()`](https://fgcz.github.io/prolfqua/reference/complete_cases.md),
[`list_to_ContrastConfiguration()`](https://fgcz.github.io/prolfqua/reference/list_to_ContrastConfiguration.md),
[`make_interaction_column()`](https://fgcz.github.io/prolfqua/reference/make_interaction_column.md),
[`make_reduced_hierarchy_config()`](https://fgcz.github.io/prolfqua/reference/make_reduced_hierarchy_config.md),
[`sample_subset()`](https://fgcz.github.io/prolfqua/reference/sample_subset.md),
[`setup_analysis()`](https://fgcz.github.io/prolfqua/reference/setup_analysis.md),
[`table_factors()`](https://fgcz.github.io/prolfqua/reference/table_factors.md),
[`table_factors_size()`](https://fgcz.github.io/prolfqua/reference/table_factors_size.md)

## Examples

``` r

bb <- sim_lfq_data_protein_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
dt <- separate_hierarchy(bb$data, bb$config)
base::setdiff(colnames(dt) ,colnames(bb$data))
#> [1] "proteinID" "idtype2"  
stopifnot(ncol(dt) >= ncol(bb$data))
```
