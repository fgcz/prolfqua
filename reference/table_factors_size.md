# Table of distinct factors with group sizes

Table of distinct factors with group sizes

## Usage

``` r
table_factors_size(
  pdata,
  file_name,
  sample_name,
  factor_keys,
  factor_keys_depth
)
```

## Arguments

- pdata:

  data.frame

- file_name:

  character - file name column

- sample_name:

  character - sample name column

- factor_keys:

  character vector - all factor column names

- factor_keys_depth:

  character vector - factor columns at current depth

## Value

a tibble with one row per combination of the `factor_keys_depth` levels
and a column `n` with the number of samples

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
[`separate_hierarchy()`](https://fgcz.github.io/prolfqua/reference/separate_hierarchy.md),
[`setup_analysis()`](https://fgcz.github.io/prolfqua/reference/setup_analysis.md),
[`table_factors()`](https://fgcz.github.io/prolfqua/reference/table_factors.md)

## Examples

``` r
istar <- sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
lfq <- LFQData$new(istar$data, istar$config)
xx <- table_factors_size(lfq$data_long(), lfq$file_name(),
  lfq$sample_name(), lfq$factor_keys(), lfq$relevant_factor_keys())
stopifnot(all(xx$n == 4))
```
