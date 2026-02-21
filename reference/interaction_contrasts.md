# Interaction contrasts (difference of differences)

Interaction contrasts (difference of differences)

## Usage

``` r
interaction_contrasts(primary_levels, secondary_levels)
```

## Arguments

- primary_levels:

  character vector of primary factor levels

- secondary_levels:

  character vector of secondary factor levels

## See also

Other contrasts:
[`annotation_add_contrasts()`](https://wolski.github.io/prolfqua/reference/annotation_add_contrasts.md),
[`generate_contrasts()`](https://wolski.github.io/prolfqua/reference/generate_contrasts.md),
[`generate_contrasts_for_factor()`](https://wolski.github.io/prolfqua/reference/generate_contrasts_for_factor.md),
[`level_specific_contrasts()`](https://wolski.github.io/prolfqua/reference/level_specific_contrasts.md),
[`main_effect_contrasts()`](https://wolski.github.io/prolfqua/reference/main_effect_contrasts.md)

## Examples

``` r
primary_levels <- c("MI", "MINOCAM")
secondary_levels <- c("T0", "T150", "T300")
interaction_contrasts(primary_levels, secondary_levels)
#> $interaction_MINOCAM_vs_MI_at_T150_vs_T0
#> [1] "(G_MINOCAM_T150 - G_MI_T150) - (G_MINOCAM_T0 - G_MI_T0)"
#> 
#> $interaction_MINOCAM_vs_MI_at_T300_vs_T0
#> [1] "(G_MINOCAM_T300 - G_MI_T300) - (G_MINOCAM_T0 - G_MI_T0)"
#> 
#> $interaction_MINOCAM_vs_MI_at_T300_vs_T150
#> [1] "(G_MINOCAM_T300 - G_MI_T300) - (G_MINOCAM_T150 - G_MI_T150)"
#> 
interaction_contrasts(secondary_levels, primary_levels)
#> $interaction_T150_vs_T0_at_MINOCAM_vs_MI
#> [1] "(G_T150_MINOCAM - G_T0_MINOCAM) - (G_T150_MI - G_T0_MI)"
#> 
#> $interaction_T300_vs_T0_at_MINOCAM_vs_MI
#> [1] "(G_T300_MINOCAM - G_T0_MINOCAM) - (G_T300_MI - G_T0_MI)"
#> 
#> $interaction_T300_vs_T150_at_MINOCAM_vs_MI
#> [1] "(G_T300_MINOCAM - G_T150_MINOCAM) - (G_T300_MI - G_T150_MI)"
#> 
```
