# Count distinct child features per hierarchy unit

Counts the number of distinct entries of the hierarchy levels below
\`hierarchy_keys_depth\` (e.g. peptides per protein). Counts
experiment-wide, whether observed or not.

## Usage

``` r
nr_features_experiment(
  data,
  hierarchy_keys,
  hierarchy_keys_depth,
  name_nr_child = "nr_features_exp"
)
```

## Arguments

- data:

  data.frame

- hierarchy_keys:

  character vector — all hierarchy column names

- hierarchy_keys_depth:

  character vector — hierarchy columns at current depth

- name_nr_child:

  character — output column name

## Value

The computed result.

## See also

Other hierarchy counting:
[`HierarchyCountsSample`](https://fgcz.github.io/prolfqua/reference/HierarchyCountsSample.md),
[`hierarchy_counts()`](https://fgcz.github.io/prolfqua/reference/hierarchy_counts.md),
[`nr_children_experiment()`](https://fgcz.github.io/prolfqua/reference/nr_children_experiment.md),
[`nr_obs_sample()`](https://fgcz.github.io/prolfqua/reference/nr_obs_sample.md),
[`summarize_hierarchy()`](https://fgcz.github.io/prolfqua/reference/summarize_hierarchy.md)

## Examples

``` r
dd <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
lfq <- LFQData$new(dd$data, dd$config)
xd <- nr_features_experiment(lfq$data_long(), lfq$hierarchy_keys(),
  lfq$relevant_hierarchy_keys())
stopifnot(min(xd$nr_features_exp) == 1)
```
