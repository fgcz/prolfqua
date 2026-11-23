# Count distinct elements for each level of hierarchy per sample

Provides wide, long, and plot views of hierarchy counts. Counts per
sample, observed values only (non-missing response and at least
`nr_children` children).

## Value

An R6 class generator.

## See also

Other summary:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md)

Other hierarchy counting:
[`nr_children_experiment()`](https://fgcz.github.io/prolfqua/reference/nr_children_experiment.md),
[`nr_features_experiment()`](https://fgcz.github.io/prolfqua/reference/nr_features_experiment.md),
[`nr_obs_sample()`](https://fgcz.github.io/prolfqua/reference/nr_obs_sample.md)

## Public fields

- `.summary`:

  summarised data frame

- `.configuration`:

  AnalysisConfiguration object

## Methods

### Public methods

- [`HierarchyCountsSample$new()`](#method-HierarchyCountsSample-initialize)

- [`HierarchyCountsSample$wide()`](#method-HierarchyCountsSample-wide)

- [`HierarchyCountsSample$long()`](#method-HierarchyCountsSample-long)

- [`HierarchyCountsSample$plot()`](#method-HierarchyCountsSample-plot)

- [`HierarchyCountsSample$clone()`](#method-HierarchyCountsSample-clone)

------------------------------------------------------------------------

### `HierarchyCountsSample$new()`

Create a new HierarchyCountsSample

#### Usage

    HierarchyCountsSample$new(pdata, configuration, nr_children = 1)

#### Arguments

- `pdata`:

  data frame

- `configuration`:

  AnalysisConfiguration

- `nr_children`:

  minimum number of children to include

------------------------------------------------------------------------

### `HierarchyCountsSample$wide()`

Return wide-format summary

#### Usage

    HierarchyCountsSample$wide()

------------------------------------------------------------------------

### `HierarchyCountsSample$long()`

Return long-format summary

#### Usage

    HierarchyCountsSample$long()

------------------------------------------------------------------------

### `HierarchyCountsSample$plot()`

Return barplot of hierarchy counts

#### Usage

    HierarchyCountsSample$plot()

------------------------------------------------------------------------

### `HierarchyCountsSample$clone()`

The objects of this class are cloneable with this method.

#### Usage

    HierarchyCountsSample$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
bb <- sim_lfq_data_protein_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 5 unconfigured columns: average_prot_abundance, sd, mean, N, Replicate.
#> setup_analysis: added 16 rows with NA response for feature and sample combinations absent from the input.
counts <- HierarchyCountsSample$new(bb$data, bb$config)
head(counts$wide())
#> # A tibble: 6 × 2
#>   sampleName protein_Id
#>   <chr>           <int>
#> 1 A_V1                9
#> 2 A_V2                9
#> 3 A_V3                7
#> 4 A_V4               10
#> 5 B_V1                9
#> 6 B_V2                9
```
