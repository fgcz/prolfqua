# Registry of available contrast facade classes

Read-only snapshot of the prolfqua facade registry. Use
[`register_facade`](https://wolski.github.io/prolfqua/reference/register_facade.md)
to add entries from downstream packages,
[`lookup_facade`](https://wolski.github.io/prolfqua/reference/lookup_facade.md)
to resolve a single entry, and
[`list_facades()`](https://wolski.github.io/prolfqua/reference/list_facades.md)
to enumerate the current registry.

## Usage

``` r
FACADE_REGISTRY
```

## Format

An object of class `facade_registry` (inherits from `list`) of length
16.

## Details

Each entry has fields `class`, `needs`, `package`, and
`needs_saint_annotation`.

## Examples

``` r
names(FACADE_REGISTRY)
#>  [1] "lm"                "lm_impute"         "lm_missing"       
#>  [4] "limma"             "limma_impute"      "limma_voom"       
#>  [7] "limma_voom_impute" "limpa"             "limpa_nested"     
#> [10] "rlm"               "deqms"             "deqms_voom"       
#> [13] "firth"             "firth_nested"      "lmer_nested"      
#> [16] "ropeca_nested"    
FACADE_REGISTRY$limma$class
#> [1] "ContrastsLimmaFacade"
```
