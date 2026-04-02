# compute contrasts

compute contrasts

## Usage

``` r
contrasts_linfct_firth(models, subject_Id = "protein_Id")
```

## Examples

``` r
mod3 <- sim_build_models_logistf(model = "parallel3", weight_missing = 1, peptide=TRUE)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
#> completing cases
#> Joining with `by = join_by(protein_Id)`
#> Joining with `by = join_by(protein_Id)`
contrasts <- c(Avs = "group_A - group_B", AvsCtrl = "group_A - group_Ctrl")
ctrpep <- ContrastsFirth$new(mod3,contrasts)
ctrpep$get_contrast_sides()
#> # A tibble: 2 × 3
#>   contrast group_1 group_2   
#>   <chr>    <chr>   <chr>     
#> 1 Avs      group_A group_B   
#> 2 AvsCtrl  group_A group_Ctrl

xx <- ctrpep$get_linfct()
models1 <- xx$models$models$models1
tmp1 <- contrasts_linfct_firth(models1)
#> contrasts_linfct_firth
models2 <- xx$models$models$models2
tmp2 <- contrasts_linfct_firth(models2)
#> contrasts_linfct_firth
stopifnot(all(dim(tmp1) > 10))
stopifnot(all(dim(tmp2) > 10))
```
