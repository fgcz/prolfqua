# Max nr_children across samples per hierarchy unit

Aggregates nr_children per sample (via nr_obs_sample), then takes the
max per hierarchy unit.

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

  character — intensity column name

- hierarchy_keys_depth:

  character vector — hierarchy columns at current depth

- file_name:

  character — file name column

- nr_children_col:

  character — nr_children column name

- name_nr_child:

  character — output column name

## Value

The computed result.

## Examples

``` r
dd <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq <- LFQData$new(dd$data, dd$config)
xd <- nr_children_experiment(lfq$data_long(), lfq$response(),
  lfq$relevant_hierarchy_keys(), lfq$file_name(), lfq$nr_children_col())
stopifnot(min(xd$nr_child_exp) == 1)
```
