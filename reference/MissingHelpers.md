# Group statistics and limit of detection

Group statistics and limit of detection

Group statistics and limit of detection

## Value

An R6 class generator.

## Details

Computes per-group statistics of the abundances and, from them, the
limit of detection (LOD) used by
[`build_model_impute`](https://fgcz.github.io/prolfqua/reference/build_model_impute.md)
to impute groups without observations.

## Public fields

- `data`:

  data

- `config`:

  config

- `prob`:

  quantile of groups with one observed value to estimate LOD

- `stats`:

  data.frame with group statistics

## Methods

### Public methods

- [`MissingHelpers$new()`](#method-MissingHelpers-new)

- [`MissingHelpers$get_stats()`](#method-MissingHelpers-get_stats)

- [`MissingHelpers$get_lod()`](#method-MissingHelpers-get_lod)

- [`MissingHelpers$clone()`](#method-MissingHelpers-clone)

------------------------------------------------------------------------

### Method [`new()`](https://rdrr.io/r/methods/new.html)

initialize

#### Usage

    MissingHelpers$new(data, config, prob = 0.5)

#### Arguments

- `data`:

  data

- `config`:

  config

- `prob`:

  default 0.5, median of groups with one observed value

------------------------------------------------------------------------

### Method `get_stats()`

get data.frame with statistics

#### Usage

    MissingHelpers$get_stats()

#### Returns

data.frame

------------------------------------------------------------------------

### Method `get_lod()`

determine limit of detection computes quantile of abundances in groups
with a single observation; without such groups, of the partly observed
groups with the fewest observations, and when every group is complete or
empty, the lowest group mean

#### Usage

    MissingHelpers$get_lod()

#### Returns

numeric LOD

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    MissingHelpers$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
dd <- prolfqua::sim_lfq_data_protein_config(Nprot = 100,weight_missing = 2)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
mh <- prolfqua::MissingHelpers$new(dd$data, dd$config, prob = 0.8)
xx <- mh$get_stats()
xx <- mh$get_lod()

dd <- prolfqua::sim_lfq_data_2factor_config(Nprot = 100,weight_missing = 0.1)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
mh <- prolfqua::MissingHelpers$new(dd$data, dd$config, prob = 0.8)
mh$get_stats()$interaction |> table()
#> 
#> TreatmentA:BackgroundX TreatmentB:BackgroundX TreatmentA:BackgroundZ 
#>                    100                    100                    100 
#> TreatmentB:BackgroundZ             TreatmentA             TreatmentB 
#>                    100                    100                    100 
#>            BackgroundX            BackgroundZ 
#>                    100                    100 
mh$get_lod()
#> [1] 15.24351
```
