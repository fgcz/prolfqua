# Registry of available contrast facade classes

Read-only snapshot of the built-in prolfqua facade registry, taken at
package load time. It is derived from the single source of truth
(`.seed_facade_registry()`); use
[`register_facade`](https://wolski.github.io/prolfqua/reference/register_facade.md)
to add entries from downstream packages,
[`lookup_facade`](https://wolski.github.io/prolfqua/reference/lookup_facade.md)
to resolve a single entry, and
[`list_facades()`](https://wolski.github.io/prolfqua/reference/list_facades.md)
to enumerate the live registry (which reflects downstream
registrations).

## Usage

``` r
FACADE_REGISTRY
```

## Format

An object of class `facade_registry` (inherits from `list`) of length
20.

## Value

A named list of built-in facade registry entries.

## Details

Each entry has fields `class`, `needs`, `package`, and
`needs_saint_annotation`.

## Examples

``` r
names(FACADE_REGISTRY)
#>  [1] "binomial_firth_nested" "binomial_nested"       "deqms"                
#>  [4] "deqms_voom"            "firth"                 "firth_nested"         
#>  [7] "hurdle_nested"         "limma"                 "limma_impute"         
#> [10] "limma_voom"            "limma_voom_impute"     "limpa"                
#> [13] "limpa_nested"          "lm"                    "lm_impute"            
#> [16] "lmer_nested"           "rfit"                  "rfit_impute"          
#> [19] "rlm"                   "ropeca_nested"        
FACADE_REGISTRY$limma$class
#> [1] "ContrastsLimmaFacade"
```
