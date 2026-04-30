# Registry of available contrast facade classes

A named list mapping short names to facade class names and their data
requirements. Each entry has:

- class:

  Character string naming the R6 facade class

- needs:

  One of `"aggregated"`, `"nested"`, or `"either"`

## Usage

``` r
FACADE_REGISTRY
```

## Format

An object of class `list` of length 14.

## Examples

``` r
names(FACADE_REGISTRY)
#>  [1] "lm"                "lm_impute"         "lm_missing"       
#>  [4] "limma"             "limma_impute"      "limma_voom"       
#>  [7] "limma_voom_impute" "limpa"             "rlm"              
#> [10] "deqms"             "deqms_voom"        "firth"            
#> [13] "lmer"              "ropeca"           
FACADE_REGISTRY$limma$class
#> [1] "ContrastsLimmaFacade"
```
