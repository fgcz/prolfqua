# Registry of available contrast facade classes

Read-only snapshot of the built-in prolfqua facade registry, taken at
package load time. It is derived from the single source of truth
(`.seed_facade_registry()`); use
[`register_facade`](https://fgcz.github.io/prolfqua/reference/register_facade.md)
to add entries from downstream packages,
[`lookup_facade`](https://fgcz.github.io/prolfqua/reference/lookup_facade.md)
to resolve a single entry, and
[`list_facades()`](https://fgcz.github.io/prolfqua/reference/list_facades.md)
to enumerate the live registry (which reflects downstream
registrations).

## Usage

``` r
FACADE_REGISTRY
```

## Value

A named list of built-in facade registry entries.

## Details

Each entry has fields `class`, `needs`, `package`, and
`needs_saint_annotation`.

## Examples

``` r
names(FACADE_REGISTRY)
#>  [1] "binomial_firth_nested"  "binomial_nested"        "deqms"                 
#>  [4] "deqms_voom"             "firth"                  "hurdle_nested"         
#>  [7] "limma"                  "limma_impute"           "limma_voom"            
#> [10] "limma_voom_impute"      "limpa"                  "limpa_nested"          
#> [13] "lm"                     "lm_impute"              "lmer_nested"           
#> [16] "rfit"                   "rfit_impute"            "rlm"                   
#> [19] "ropecaV2_impute_nested" "ropecaV2_nested"        "ropeca_impute_nested"  
#> [22] "ropeca_nested"         
FACADE_REGISTRY$limma$class
#> [1] "ContrastsLimmaFacade"
```
