# Hierarchy counts per sample

Hierarchy counts per sample

## Usage

``` r
hierarchy_counts_sample(pdata, configuration, nr_children = 1)
```

## Arguments

- pdata:

  data.frame

- configuration:

  AnalysisConfiguration

- nr_children:

  minimum number of children

## Value

[`HierarchyCountsSample`](https://wolski.github.io/prolfqua/reference/HierarchyCountsSample.md)
R6 object

## See also

Other summary:
[`HierarchyCountsSample`](https://wolski.github.io/prolfqua/reference/HierarchyCountsSample.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`hierarchy_counts()`](https://wolski.github.io/prolfqua/reference/hierarchy_counts.md),
[`summarize_hierarchy()`](https://wolski.github.io/prolfqua/reference/summarize_hierarchy.md)

## Examples

``` r
bb <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done

config <- bb$config
data <- bb$data
res <- hierarchy_counts_sample(data, config, nr_children = 1)
x <- res$long()
# filters on peptide level
res <- hierarchy_counts_sample(data, config, nr_children = 2)
x2 <- res$long()
# filters on protein level based on peptide count
bb <- prolfqua::sim_lfq_data_protein_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
res <- hierarchy_counts_sample(bb$data, bb$config, nr_children = 2)
x1 <- res$wide()
res <- hierarchy_counts_sample(bb$data, bb$config, nr_children = 1)
x2 <- res$wide()
x1$nr_children <- 2
x2$nr_children <- 1
xl <- dplyr::bind_rows(x1, x2)

xl$nr_children |> table()
#> 
#>  1  2 
#> 12 12 
nudgeval <-  -mean(xl$protein_Id) * 0.05
ggplot2::ggplot(xl,
  ggplot2::aes(x = sampleName, y = protein_Id, fill = as.character(nr_children))) +
 ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge())

```
