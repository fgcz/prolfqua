# Make reduced hierarchy configuration

Make reduced hierarchy configuration

## Usage

``` r
make_reduced_hierarchy_config(config, work_intensity, hierarchy)
```

## Arguments

- config:

  AnalysisConfiguration

- work_intensity:

  work intensity column

- hierarchy:

  new reduced hierarchy

## Value

AnalysisConfiguration with reduced hierarchy

## See also

Other configuration:
[`AnalysisConfiguration`](https://fgcz.github.io/prolfqua/reference/AnalysisConfiguration.md),
[`ContrastConfiguration`](https://fgcz.github.io/prolfqua/reference/ContrastConfiguration.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`R6_extract_values()`](https://fgcz.github.io/prolfqua/reference/R6_extract_values.md),
[`complete_cases()`](https://fgcz.github.io/prolfqua/reference/complete_cases.md),
[`list_to_ContrastConfiguration()`](https://fgcz.github.io/prolfqua/reference/list_to_ContrastConfiguration.md),
[`make_interaction_column()`](https://fgcz.github.io/prolfqua/reference/make_interaction_column.md),
[`sample_subset()`](https://fgcz.github.io/prolfqua/reference/sample_subset.md),
[`separate_hierarchy()`](https://fgcz.github.io/prolfqua/reference/separate_hierarchy.md),
[`setup_analysis()`](https://fgcz.github.io/prolfqua/reference/setup_analysis.md),
[`table_factors()`](https://fgcz.github.io/prolfqua/reference/table_factors.md),
[`table_factors_size()`](https://fgcz.github.io/prolfqua/reference/table_factors_size.md)

## Examples

``` r
bb <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done

red <- make_reduced_hierarchy_config(bb$config,
 "testintensity",
 bb$config$hierarchy[1])
stopifnot(red$get_response() == "testintensity")
stopifnot(length(red$hierarchy) == 1)
```
