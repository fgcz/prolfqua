# Interaction contrasts of a two-factor design

For every pair of primary levels and every pair of secondary levels, the
difference of differences: does the primary difference change between
the two secondary levels?

## Usage

``` r
interaction_contrasts(primary_levels, secondary_levels)
```

## Arguments

- primary_levels:

  character vector of primary factor levels

- secondary_levels:

  character vector of secondary factor levels

## Value

A named list of contrast strings, named
`interaction_<high>_vs_<low>_at_<s2>_vs_<s1>`.

## See also

Other contrasts:
[`annotation_add_contrasts()`](https://fgcz.github.io/prolfqua/reference/annotation_add_contrasts.md),
[`generate_contrasts()`](https://fgcz.github.io/prolfqua/reference/generate_contrasts.md),
[`generate_contrasts_for_factor()`](https://fgcz.github.io/prolfqua/reference/generate_contrasts_for_factor.md),
[`level_specific_contrasts()`](https://fgcz.github.io/prolfqua/reference/level_specific_contrasts.md),
[`main_effect_contrasts()`](https://fgcz.github.io/prolfqua/reference/main_effect_contrasts.md)

## Examples

``` r
primary_levels <- c("KO", "WT")
secondary_levels <- c("T0", "T150", "T300")
interaction_contrasts(primary_levels, secondary_levels)
#> $interaction_WT_vs_KO_at_T150_vs_T0
#> [1] "(G_WT_T150 - G_KO_T150) - (G_WT_T0 - G_KO_T0)"
#> 
#> $interaction_WT_vs_KO_at_T300_vs_T0
#> [1] "(G_WT_T300 - G_KO_T300) - (G_WT_T0 - G_KO_T0)"
#> 
#> $interaction_WT_vs_KO_at_T300_vs_T150
#> [1] "(G_WT_T300 - G_KO_T300) - (G_WT_T150 - G_KO_T150)"
#> 
interaction_contrasts(secondary_levels, primary_levels)
#> $interaction_T150_vs_T0_at_WT_vs_KO
#> [1] "(G_T150_WT - G_T0_WT) - (G_T150_KO - G_T0_KO)"
#> 
#> $interaction_T300_vs_T0_at_WT_vs_KO
#> [1] "(G_T300_WT - G_T0_WT) - (G_T300_KO - G_T0_KO)"
#> 
#> $interaction_T300_vs_T150_at_WT_vs_KO
#> [1] "(G_T300_WT - G_T150_WT) - (G_T300_KO - G_T150_KO)"
#> 
```
