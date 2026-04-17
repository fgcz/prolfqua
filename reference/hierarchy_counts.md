# Count distinct elements for each level of hierarchy and isotope

E.g. number of proteins, peptides, precursors in the dataset

## Usage

``` r
hierarchy_counts(pdata, hierarchy_keys, isotope_label = "isotopeLabel")
```

## Arguments

- pdata:

  data.frame

- hierarchy_keys:

  character vector — all hierarchy column names

- isotope_label:

  character — isotope label column name

## See also

Other summary:
[`HierarchyCountsSample`](https://wolski.github.io/prolfqua/reference/HierarchyCountsSample.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`hierarchy_counts_sample()`](https://wolski.github.io/prolfqua/reference/hierarchy_counts_sample.md),
[`summarize_hierarchy()`](https://wolski.github.io/prolfqua/reference/summarize_hierarchy.md)

## Examples

``` r
bb <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq <- LFQData$new(bb$data, bb$config)

x <- hierarchy_counts(lfq$data_long(), lfq$hierarchy_keys(), lfq$isotope_label())
x$protein_Id
#> [1] 10
stopifnot(ncol(x) == length(lfq$hierarchy_keys()) + 1)
# select non existing protein
data0 <- lfq$data_long() |> dplyr::filter(protein_Id == "XYZ")
tmp <- hierarchy_counts(data0, lfq$hierarchy_keys(), lfq$isotope_label())
stopifnot(nrow(tmp) == 0)
```
