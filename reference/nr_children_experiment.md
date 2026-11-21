# Max nr_children across samples per hierarchy unit

Aggregates nr_children per sample (via nr_obs_sample), then takes the
max per hierarchy unit. Counts experiment-wide, observed values only.

## Usage

``` r
nr_children_experiment(
  data,
  response,
  hierarchy_keys_depth,
  file_name,
  nr_children_col,
  name_nr_child = "nr_child_exp"
)
```

## Arguments

- data:

  data.frame

- response:

  character - intensity column name

- hierarchy_keys_depth:

  character vector - hierarchy columns at current depth

- file_name:

  character - file name column

- nr_children_col:

  character - nr_children column name

- name_nr_child:

  character - output column name

## Value

The computed result.

## See also

Other hierarchy counting:
[`HierarchyCountsSample`](https://fgcz.github.io/prolfqua/reference/HierarchyCountsSample.md),
[`hierarchy_counts()`](https://fgcz.github.io/prolfqua/reference/hierarchy_counts.md),
[`nr_features_experiment()`](https://fgcz.github.io/prolfqua/reference/nr_features_experiment.md),
[`nr_obs_sample()`](https://fgcz.github.io/prolfqua/reference/nr_obs_sample.md),
[`summarize_hierarchy()`](https://fgcz.github.io/prolfqua/reference/summarize_hierarchy.md)

## Examples

``` r
dd <- prolfqua::sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
lfq <- LFQData$new(dd$data, dd$config)
xd <- nr_children_experiment(lfq$data_long(), lfq$response(),
  lfq$relevant_hierarchy_keys(), lfq$file_name(), lfq$nr_children_col())
stopifnot(min(xd$nr_child_exp) == 1)
```
