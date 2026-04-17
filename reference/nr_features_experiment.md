# Count distinct child features per hierarchy unit

Counts the number of distinct child-level entries (e.g. peptides per
protein).

## Usage

``` r
nr_features_experiment(
  data,
  hierarchy_keys,
  hierarchy_keys_depth,
  name_nr_child = "nr_child_exp"
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

## Examples

``` r
dd <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq <- LFQData$new(dd$data, dd$config)
xd <- nr_features_experiment(lfq$data_long(), lfq$hierarchy_keys(),
  lfq$relevant_hierarchy_keys())
stopifnot(min(xd$nr_child_exp) == 1)
```
