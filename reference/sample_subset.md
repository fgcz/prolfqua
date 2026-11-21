# Sample subset of proteins/peptides/precursors

Sample subset of proteins/peptides/precursors

## Usage

``` r
sample_subset(size, pdata, hierarchy_keys_depth)
```

## Arguments

- size:

  size of sample

- pdata:

  tidy table

- hierarchy_keys_depth:

  character vector - hierarchy columns at current depth

## Value

a tibble with all rows of `pdata` for `size` randomly drawn hierarchy
elements

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
[`separate_hierarchy()`](https://fgcz.github.io/prolfqua/reference/separate_hierarchy.md),
[`setup_analysis()`](https://fgcz.github.io/prolfqua/reference/setup_analysis.md),
[`table_factors()`](https://fgcz.github.io/prolfqua/reference/table_factors.md),
[`table_factors_size()`](https://fgcz.github.io/prolfqua/reference/table_factors_size.md)

## Examples

``` r
bb <- sim_lfq_data_peptide_config(n_proteins = 5)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 39 rows with NA response for feature and sample combinations absent from the input.
subset <- sample_subset(2, bb$data, bb$config$hierarchy_keys_depth())
length(unique(subset[[bb$config$hierarchy_keys_depth()[1]]]))
#> [1] 2
```
