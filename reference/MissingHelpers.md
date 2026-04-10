# compute group mean by LOD

compute group mean by LOD

compute group mean by LOD

## Details

weight lod by nr of NA's \$(LOD \* nrNas + meanAbundance
\*nrObs)/(nrMeasured)\$

## Public fields

- `data`:

  data

- `config`:

  config

- `prob`:

  quantile of groups with one observed value to estimate LOD

- `stats`:

  data.frame with group statistics

- `weighted`:

  should we weight the LOD

## Methods

### Public methods

- [`MissingHelpers$new()`](#method-MissingHelpers-new)

- [`MissingHelpers$get_stats()`](#method-MissingHelpers-get_stats)

- [`MissingHelpers$get_lod()`](#method-MissingHelpers-get_lod)

- [`MissingHelpers$impute_weighted_lod()`](#method-MissingHelpers-impute_weighted_lod)

- [`MissingHelpers$impute_lod()`](#method-MissingHelpers-impute_lod)

- [`MissingHelpers$get_poolvar()`](#method-MissingHelpers-get_poolvar)

- [`MissingHelpers$get_contrast_estimates()`](#method-MissingHelpers-get_contrast_estimates)

- [`MissingHelpers$get_contrasts()`](#method-MissingHelpers-get_contrasts)

- [`MissingHelpers$clone()`](#method-MissingHelpers-clone)

------------------------------------------------------------------------

### Method [`new()`](https://rdrr.io/r/methods/new.html)

initialize

#### Usage

    MissingHelpers$new(data, config, prob = 0.5, weighted = TRUE)

#### Arguments

- `data`:

  data

- `config`:

  config

- `prob`:

  default 0.5, median of groups with one observed value

- `weighted`:

  should group average be computed used weighting, default TRUE.

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
with a single observation

#### Usage

    MissingHelpers$get_lod()

#### Returns

integer LOD

------------------------------------------------------------------------

### Method `impute_weighted_lod()`

compute group averages using weighted lod

#### Usage

    MissingHelpers$impute_weighted_lod()

------------------------------------------------------------------------

### Method `impute_lod()`

if group average absent substitute with LOD

#### Usage

    MissingHelpers$impute_lod()

------------------------------------------------------------------------

### Method `get_poolvar()`

compute pooled var per protein

#### Usage

    MissingHelpers$get_poolvar(prob = 0.75)

#### Arguments

- `prob`:

  prob of sd from proteins where it can be computed

------------------------------------------------------------------------

### Method `get_contrast_estimates()`

get contrast estimates

#### Usage

    MissingHelpers$get_contrast_estimates(Contrasts)

#### Arguments

- `Contrasts`:

  named array with contrasts

------------------------------------------------------------------------

### Method `get_contrasts()`

compute contrasts

#### Usage

    MissingHelpers$get_contrasts(Contrasts, confint = 0.95, all = FALSE)

#### Arguments

- `Contrasts`:

  vector with contrasts

- `confint`:

  compute confint

- `all`:

  return all columns, default FALSE

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
Contrasts <- c("group.b-a" = "group_A - group_B", "group.a-ctrl" = "group_A - group_Ctrl")
dd <- prolfqua::sim_lfq_data_protein_config(Nprot = 100,weight_missing = 2)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
mh <- prolfqua::MissingHelpers$new(dd$data, dd$config, prob = 0.8,weighted = TRUE)
xx <- mh$get_stats()
xx <- mh$get_lod()
xx <- mh$impute_weighted_lod()
xx <- mh$impute_lod()
xx <- mh$get_poolvar()
bb <- mh$get_contrast_estimates(Contrasts)
#> group.b-a=group_A - group_B
#> group.a-ctrl=group_A - group_Ctrl
#> group.b-a=group_A - group_B
#> group.a-ctrl=group_A - group_Ctrl
#> group.b-a=group_A - group_B
#> group.a-ctrl=group_A - group_Ctrl
mh$get_contrasts(Contrasts)
#> group.b-a=group_A - group_B
#> group.a-ctrl=group_A - group_Ctrl
#> group.b-a=group_A - group_B
#> group.a-ctrl=group_A - group_Ctrl
#> group.b-a=group_A - group_B
#> group.a-ctrl=group_A - group_Ctrl
#> # A tibble: 200 × 18
#>    protein_Id  meanAbundanceImp_group_1 meanAbundanceImp_group_2 estimate
#>    <chr>                          <dbl>                    <dbl>    <dbl>
#>  1 0EfVhX~3967                     22.0                     19.2    2.84 
#>  2 0YSKpy~2865                     19.5                     20.0    0    
#>  3 0m5WN4~6025                     19.5                     23.1   -3.59 
#>  4 3QLHfm~8938                     22.2                     21.8    0.365
#>  5 3QYop0~7543                     22.5                     22.1    0.433
#>  6 76k03k~7094                     25.3                     24.3    1.02 
#>  7 7QuTub~1867                     20.4                     24.6   -4.12 
#>  8 7cbcrd~7351                     26.7                     23.6    3.17 
#>  9 7soopj~5352                     20.1                     20.0    0.113
#> 10 7zeekV~7127                     20.6                     20.8   -0.196
#> # ℹ 190 more rows
#> # ℹ 14 more variables: group_1_name <chr>, group_2_name <chr>, contrast <chr>,
#> #   avgAbd <dbl>, indic <dbl>, nrMeasured_group_1 <int>,
#> #   nrMeasured_group_2 <int>, df <dbl>, sd <dbl>, sdT <dbl>, statistic <dbl>,
#> #   p.value <dbl>, conf.low <dbl>, conf.high <dbl>

dd <- prolfqua::sim_lfq_data_2factor_config(Nprot = 100,weight_missing = 0.1)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done

Contrasts <- c("c1" = "TreatmentA - TreatmentB",
               "C2" = "BackgroundX- BackgroundZ",
               "c3" = "`TreatmentA:BackgroundX` - `TreatmentA:BackgroundZ`",
               "c4" = "`TreatmentB:BackgroundX` - `TreatmentB:BackgroundZ`"
               )
mh <- prolfqua::MissingHelpers$new(dd$data, dd$config, prob = 0.8,weighted = TRUE)
mh$get_stats()$interaction |> table()
#> 
#> TreatmentA:BackgroundX TreatmentB:BackgroundX TreatmentA:BackgroundZ 
#>                    100                    100                    100 
#> TreatmentB:BackgroundZ             TreatmentA             TreatmentB 
#>                    100                    100                    100 
#>            BackgroundX            BackgroundZ 
#>                    100                    100 
mh$get_contrast_estimates(Contrasts)
#> c1=TreatmentA - TreatmentB
#> C2=BackgroundX- BackgroundZ
#> c3=`TreatmentA:BackgroundX` - `TreatmentA:BackgroundZ`
#> c4=`TreatmentB:BackgroundX` - `TreatmentB:BackgroundZ`
#> c1=TreatmentA - TreatmentB
#> C2=BackgroundX- BackgroundZ
#> c3=`TreatmentA:BackgroundX` - `TreatmentA:BackgroundZ`
#> c4=`TreatmentB:BackgroundX` - `TreatmentB:BackgroundZ`
#> c1=TreatmentA - TreatmentB
#> C2=BackgroundX- BackgroundZ
#> c3=`TreatmentA:BackgroundX` - `TreatmentA:BackgroundZ`
#> c4=`TreatmentB:BackgroundX` - `TreatmentB:BackgroundZ`
#> # A tibble: 400 × 11
#>    protein_Id  meanAbundanceImp_group_1 meanAbundanceImp_group_2 estimate
#>    <chr>                          <dbl>                    <dbl>    <dbl>
#>  1 0EfVhX~3967                     21.1                     20.5    0.608
#>  2 0YSKpy~2865                     18.3                     17.1    1.15 
#>  3 0m5WN4~6025                     18.1                     22.2   -4.09 
#>  4 3QLHfm~8938                     23.2                     22.5    0.712
#>  5 3QYop0~7543                     23.3                     23.7   -0.310
#>  6 76k03k~7094                     24.3                     24.6   -0.343
#>  7 7QuTub~1867                     21.6                     25.7   -4.14 
#>  8 7cbcrd~7351                     25.6                     24.8    0.804
#>  9 7soopj~5352                     19.7                     21.6   -1.90 
#> 10 7zeekV~7127                     21.1                     20.3    0.720
#> # ℹ 390 more rows
#> # ℹ 7 more variables: group_1_name <chr>, group_2_name <chr>, contrast <chr>,
#> #   avgAbd <dbl>, indic <dbl>, nrMeasured_group_1 <int>,
#> #   nrMeasured_group_2 <int>
```
