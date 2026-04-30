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

  character vector — hierarchy columns at current depth

## See also

Other configuration:
[`AnalysisConfiguration`](https://wolski.github.io/prolfqua/reference/AnalysisConfiguration.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`R6_extract_values()`](https://wolski.github.io/prolfqua/reference/R6_extract_values.md),
[`complete_cases()`](https://wolski.github.io/prolfqua/reference/complete_cases.md),
[`make_interaction_column()`](https://wolski.github.io/prolfqua/reference/make_interaction_column.md),
[`make_reduced_hierarchy_config()`](https://wolski.github.io/prolfqua/reference/make_reduced_hierarchy_config.md),
[`separate_hierarchy()`](https://wolski.github.io/prolfqua/reference/separate_hierarchy.md),
[`setup_analysis()`](https://wolski.github.io/prolfqua/reference/setup_analysis.md),
[`table_factors()`](https://wolski.github.io/prolfqua/reference/table_factors.md),
[`table_factors_size()`](https://wolski.github.io/prolfqua/reference/table_factors_size.md)

## Examples

``` r
bb <- sim_lfq_data_peptide_config(Nprot = 5)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
subset <- sample_subset(2, bb$data, bb$config$hierarchy_keys_depth())
#> Sampling 2protein_Id
#> Joining with `by = join_by(protein_Id)`
length(unique(subset[[bb$config$hierarchy_keys_depth()[1]]]))
#> [1] 2
```
