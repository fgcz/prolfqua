# Count distinct elements for each level of hierarchy and istope

E.g. number of proteins, peptides, precursors in the dataset

## Usage

``` r
hierarchy_counts(pdata, config)
```

## Arguments

- pdata:

  data.frame

- config:

  AnalysisConfiguration

## See also

Other summary:
[`HierarchyCountsSample`](https://wolski.github.io/prolfqua/reference/HierarchyCountsSample.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`hierarchy_counts_sample()`](https://wolski.github.io/prolfqua/reference/hierarchy_counts_sample.md),
[`summarize_hierarchy()`](https://wolski.github.io/prolfqua/reference/summarize_hierarchy.md)

## Examples

``` r
bb <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done

config <- bb$config$clone(deep=TRUE)
data <- bb$data

x <- hierarchy_counts(data, config)
x$protein_Id
#> [1] 10
stopifnot(ncol(x) == length(config$hierarchy_keys()) + 1)
# select non existing protein
data0 <- data |> dplyr::filter( protein_Id == "XYZ")
tmp <- hierarchy_counts(data0, config)
stopifnot(nrow(tmp) == 0)
```
