# Count distinct elements for each level of hierarchy per sample

Count distinct elements for each level of hierarchy per sample

Count distinct elements for each level of hierarchy per sample

## Value

An R6 class generator.

## Details

Provides wide, long, and plot views of hierarchy counts.

## See also

Other summary:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`hierarchy_counts()`](https://wolski.github.io/prolfqua/reference/hierarchy_counts.md),
[`hierarchy_counts_sample()`](https://wolski.github.io/prolfqua/reference/hierarchy_counts_sample.md),
[`summarize_hierarchy()`](https://wolski.github.io/prolfqua/reference/summarize_hierarchy.md)

## Public fields

- `.summary`:

  summarised data frame

- `.configuration`:

  AnalysisConfiguration object

## Methods

### Public methods

- [`HierarchyCountsSample$new()`](#method-HierarchyCountsSample-new)

- [`HierarchyCountsSample$wide()`](#method-HierarchyCountsSample-wide)

- [`HierarchyCountsSample$long()`](#method-HierarchyCountsSample-long)

- [`HierarchyCountsSample$plot()`](#method-HierarchyCountsSample-plot)

- [`HierarchyCountsSample$clone()`](#method-HierarchyCountsSample-clone)

------------------------------------------------------------------------

### Method `new()`

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

### Method `wide()`

Return wide-format summary

#### Usage

    HierarchyCountsSample$wide()

------------------------------------------------------------------------

### Method `long()`

Return long-format summary

#### Usage

    HierarchyCountsSample$long()

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Return barplot of hierarchy counts

#### Usage

    HierarchyCountsSample$plot()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    HierarchyCountsSample$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
bb <- sim_lfq_data_protein_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
counts <- HierarchyCountsSample$new(bb$data, bb$config)
head(counts$wide())
#> # A tibble: 6 × 3
#> # Groups:   isotopeLabel [1]
#>   isotopeLabel sampleName protein_Id
#>   <chr>        <chr>           <int>
#> 1 light        A_V1                9
#> 2 light        A_V2                9
#> 3 light        A_V3                7
#> 4 light        A_V4               10
#> 5 light        B_V1                9
#> 6 light        B_V2                9
```
