# Main-effect contrasts of a two-factor design

For every pair of primary levels, the difference of their means averaged
over all secondary levels.

## Usage

``` r
main_effect_contrasts(primary_levels, secondary_levels)
```

## Arguments

- primary_levels:

  character vector of primary factor levels

- secondary_levels:

  character vector of secondary factor levels

## Value

A named list of contrast strings, one per pair of primary levels, named
`<high>_vs_<low>`; groups are named by
[`group_label()`](https://fgcz.github.io/prolfqua/reference/group_label.md).

## See also

Other contrasts:
[`annotation_add_contrasts()`](https://fgcz.github.io/prolfqua/reference/annotation_add_contrasts.md),
[`generate_contrasts()`](https://fgcz.github.io/prolfqua/reference/generate_contrasts.md),
[`generate_contrasts_for_factor()`](https://fgcz.github.io/prolfqua/reference/generate_contrasts_for_factor.md),
[`interaction_contrasts()`](https://fgcz.github.io/prolfqua/reference/interaction_contrasts.md),
[`level_specific_contrasts()`](https://fgcz.github.io/prolfqua/reference/level_specific_contrasts.md)

## Examples

``` r
primary_levels <- c("KO", "WT")
secondary_levels <- c("T150", "T0", "T300")
main_effect_contrasts(primary_levels, secondary_levels)
#> $WT_vs_KO
#> [1] "(G_WT_T150 + G_WT_T0 + G_WT_T300)/3 - (G_KO_T150 + G_KO_T0 + G_KO_T300)/3"
#> 
main_effect_contrasts(secondary_levels, primary_levels)
#> $T0_vs_T150
#> [1] "(G_T0_KO + G_T0_WT)/2 - (G_T150_KO + G_T150_WT)/2"
#> 
#> $T300_vs_T150
#> [1] "(G_T300_KO + G_T300_WT)/2 - (G_T150_KO + G_T150_WT)/2"
#> 
#> $T300_vs_T0
#> [1] "(G_T300_KO + G_T300_WT)/2 - (G_T0_KO + G_T0_WT)/2"
#> 
```
