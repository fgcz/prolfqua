# Complete cases

The tidy table does not need to contain missing data. This function
re-establishes the missing observations in a sample.

## Usage

``` r
complete_cases(lfqdata)
```

## Arguments

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

## Value

a tibble of the long-format data with one row for every combination of
sample and hierarchy element; added rows have `NA` intensities

## See also

Other configuration:
[`AnalysisConfiguration`](https://fgcz.github.io/prolfqua/reference/AnalysisConfiguration.md),
[`ContrastConfiguration`](https://fgcz.github.io/prolfqua/reference/ContrastConfiguration.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`list_to_contrast_configuration()`](https://fgcz.github.io/prolfqua/reference/list_to_contrast_configuration.md),
[`make_interaction_column()`](https://fgcz.github.io/prolfqua/reference/make_interaction_column.md),
[`make_reduced_hierarchy_config()`](https://fgcz.github.io/prolfqua/reference/make_reduced_hierarchy_config.md),
[`r6_extract_values()`](https://fgcz.github.io/prolfqua/reference/r6_extract_values.md),
[`sample_subset()`](https://fgcz.github.io/prolfqua/reference/sample_subset.md),
[`separate_hierarchy()`](https://fgcz.github.io/prolfqua/reference/separate_hierarchy.md),
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
lfq <- LFQData$new(bb$data, bb$config)
xx <- complete_cases(lfq)
stopifnot(nrow(bb$data) <= nrow(xx))
```
