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
19.

## Details

Each entry has fields `class`, `needs`, `package`, and
`needs_saint_annotation`.

## Examples

``` r
names(FACADE_REGISTRY)
#>  [1] "binomial_nested"   "deqms"             "deqms_voom"       
#>  [4] "firth"             "firth_nested"      "limma"            
#>  [7] "limma_impute"      "limma_voom"        "limma_voom_impute"
#> [10] "limpa"             "limpa_nested"      "lm"               
#> [13] "lm_impute"         "lm_missing"        "lmer_nested"      
#> [16] "rfit"              "rfit_impute"       "rlm"              
#> [19] "ropeca_nested"    
FACADE_REGISTRY$limma$class
#> [1] "ContrastsLimmaFacade"
```
