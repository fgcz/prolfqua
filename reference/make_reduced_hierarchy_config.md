# Make reduced hierarchy configuration

Make reduced hierarchy configuration

## Usage

``` r
make_reduced_hierarchy_config(config, workIntensity, hierarchy)
```

## Arguments

- config:

  AnalysisConfiguration

- workIntensity:

  work intensity column

- hierarchy:

  new reduced hierarchy

## Value

AnalysisConfiguration with reduced hieararchy

## See also

Other configuration:
[`AnalysisConfiguration`](https://wolski.github.io/prolfqua/reference/AnalysisConfiguration.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`R6_extract_values()`](https://wolski.github.io/prolfqua/reference/R6_extract_values.md),
[`complete_cases()`](https://wolski.github.io/prolfqua/reference/complete_cases.md),
[`concrete_AnalysisConfiguration`](https://wolski.github.io/prolfqua/reference/concrete_AnalysisConfiguration.md),
[`make_interaction_column()`](https://wolski.github.io/prolfqua/reference/make_interaction_column.md),
[`sample_subset()`](https://wolski.github.io/prolfqua/reference/sample_subset.md),
[`separate_hierarchy()`](https://wolski.github.io/prolfqua/reference/separate_hierarchy.md),
[`setup_analysis()`](https://wolski.github.io/prolfqua/reference/setup_analysis.md),
[`table_factors()`](https://wolski.github.io/prolfqua/reference/table_factors.md),
[`table_factors_size()`](https://wolski.github.io/prolfqua/reference/table_factors_size.md)

## Examples

``` r
bb <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done

red <- make_reduced_hierarchy_config(bb$config,
 "testintensity",
 bb$config$hierarchy[1])
stopifnot(red$get_response() == "testintensity")
stopifnot(length(red$hierarchy) == 1)
```
