# main effects contrasts

main effects contrasts

## Usage

``` r
main_effect_contrasts(primary_levels, secondary_levels)
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
[`interaction_contrasts()`](https://wolski.github.io/prolfqua/reference/interaction_contrasts.md),
[`level_specific_contrasts()`](https://wolski.github.io/prolfqua/reference/level_specific_contrasts.md)

## Examples

``` r
primary_levels <- c("MI", "MINOCAM")
secondary_levels <- c("T150", "T0", "T300")
main_effect_contrasts(primary_levels, secondary_levels)
#> $MINOCAM_vs_MI
#> [1] "( (G_MINOCAM_T150 + G_MINOCAM_T0 + G_MINOCAM_T300)/3 - (G_MI_T150 + G_MI_T0 + G_MI_T300)/3 )"
#> 
main_effect_contrasts(secondary_levels, primary_levels)
#> $T0_vs_T150
#> [1] "( (G_T0_MI + G_T0_MINOCAM)/2 - (G_T150_MI + G_T150_MINOCAM)/2 )"
#> 
#> $T300_vs_T150
#> [1] "( (G_T300_MI + G_T300_MINOCAM)/2 - (G_T150_MI + G_T150_MINOCAM)/2 )"
#> 
#> $T300_vs_T0
#> [1] "( (G_T300_MI + G_T300_MINOCAM)/2 - (G_T0_MI + G_T0_MINOCAM)/2 )"
#> 
```
