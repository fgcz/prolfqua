# Count observations per experiment (max nr_children across samples)

Count observations per experiment (max nr_children across samples)

## Usage

``` r
nr_obs_experiment(
  data,
  hierarchy_keys,
  hierarchy_keys_depth,
  nr_children_col,
  response = NULL,
  file_name = NULL,
  from_children = TRUE,
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

- nr_children_col:

  character — nr_children column name

- response:

  character — intensity column name (needed when from_children = TRUE)

- file_name:

  character — file name column (needed when from_children = TRUE)

- from_children:

  TRUE compute from nr_children column, FALSE count distinct hierarchy
  entries

- name_nr_child:

  how to name output column

## Examples

``` r
dd <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq <- LFQData$new(dd$data, dd$config)
xd <- nr_obs_experiment(lfq$get_data(), lfq$hierarchy_keys(),
  lfq$relevant_hierarchy_keys(), lfq$nr_children_col(),
  response = lfq$response(), file_name = lfq$file_name())
stopifnot(min(xd$nr_child_exp) == 1)
xd2 <- nr_obs_experiment(lfq$get_data(), lfq$hierarchy_keys(),
  lfq$relevant_hierarchy_keys(), lfq$nr_children_col(), from_children = FALSE)
```
