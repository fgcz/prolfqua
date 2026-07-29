# Combined generate_contrasts

Combined generate_contrasts

## Usage

``` r
generate_contrasts(primary_levels, secondary_levels, interactions = TRUE)
```

## Arguments

- primary_levels:

  character vector of primary factor levels

- secondary_levels:

  character vector of secondary factor levels

- interactions:

  logical, if TRUE include interaction contrasts

## Value

Contrast definitions or contrast results.

## See also

Other contrasts:
[`annotation_add_contrasts()`](https://wolski.github.io/prolfqua/reference/annotation_add_contrasts.md),
[`generate_contrasts_for_factor()`](https://wolski.github.io/prolfqua/reference/generate_contrasts_for_factor.md),
[`interaction_contrasts()`](https://wolski.github.io/prolfqua/reference/interaction_contrasts.md),
[`level_specific_contrasts()`](https://wolski.github.io/prolfqua/reference/level_specific_contrasts.md),
[`main_effect_contrasts()`](https://wolski.github.io/prolfqua/reference/main_effect_contrasts.md)

## Examples

``` r
primary_levels <- c("MI", "MINOCAM")
secondary_levels <- c("T0", "T150", "T300")
generate_contrasts(primary_levels, secondary_levels)
#>                                                                        ContrastName
#> MINOCAM_vs_MI                                                         MINOCAM_vs_MI
#> MINOCAM_vs_MI_at_T0                                             MINOCAM_vs_MI_at_T0
#> MINOCAM_vs_MI_at_T150                                         MINOCAM_vs_MI_at_T150
#> MINOCAM_vs_MI_at_T300                                         MINOCAM_vs_MI_at_T300
#> interaction_MINOCAM_vs_MI_at_T150_vs_T0     interaction_MINOCAM_vs_MI_at_T150_vs_T0
#> interaction_MINOCAM_vs_MI_at_T300_vs_T0     interaction_MINOCAM_vs_MI_at_T300_vs_T0
#> interaction_MINOCAM_vs_MI_at_T300_vs_T150 interaction_MINOCAM_vs_MI_at_T300_vs_T150
#>                                                                                                                           Contrast
#> MINOCAM_vs_MI                             (G_MINOCAM_T0 + G_MINOCAM_T150 + G_MINOCAM_T300)/3 - (G_MI_T0 + G_MI_T150 + G_MI_T300)/3
#> MINOCAM_vs_MI_at_T0                                                                                         G_MINOCAM_T0 - G_MI_T0
#> MINOCAM_vs_MI_at_T150                                                                                   G_MINOCAM_T150 - G_MI_T150
#> MINOCAM_vs_MI_at_T300                                                                                   G_MINOCAM_T300 - G_MI_T300
#> interaction_MINOCAM_vs_MI_at_T150_vs_T0                                    (G_MINOCAM_T150 - G_MI_T150) - (G_MINOCAM_T0 - G_MI_T0)
#> interaction_MINOCAM_vs_MI_at_T300_vs_T0                                    (G_MINOCAM_T300 - G_MI_T300) - (G_MINOCAM_T0 - G_MI_T0)
#> interaction_MINOCAM_vs_MI_at_T300_vs_T150                              (G_MINOCAM_T300 - G_MI_T300) - (G_MINOCAM_T150 - G_MI_T150)
generate_contrasts(secondary_levels, primary_levels)
#>                                                                        ContrastName
#> T150_vs_T0                                                               T150_vs_T0
#> T300_vs_T0                                                               T300_vs_T0
#> T300_vs_T150                                                           T300_vs_T150
#> T150_vs_T0_at_MI                                                   T150_vs_T0_at_MI
#> T150_vs_T0_at_MINOCAM                                         T150_vs_T0_at_MINOCAM
#> T300_vs_T0_at_MI                                                   T300_vs_T0_at_MI
#> T300_vs_T0_at_MINOCAM                                         T300_vs_T0_at_MINOCAM
#> T300_vs_T150_at_MI                                               T300_vs_T150_at_MI
#> T300_vs_T150_at_MINOCAM                                     T300_vs_T150_at_MINOCAM
#> interaction_T150_vs_T0_at_MINOCAM_vs_MI     interaction_T150_vs_T0_at_MINOCAM_vs_MI
#> interaction_T300_vs_T0_at_MINOCAM_vs_MI     interaction_T300_vs_T0_at_MINOCAM_vs_MI
#> interaction_T300_vs_T150_at_MINOCAM_vs_MI interaction_T300_vs_T150_at_MINOCAM_vs_MI
#>                                                                                                  Contrast
#> T150_vs_T0                                    (G_T150_MI + G_T150_MINOCAM)/2 - (G_T0_MI + G_T0_MINOCAM)/2
#> T300_vs_T0                                    (G_T300_MI + G_T300_MINOCAM)/2 - (G_T0_MI + G_T0_MINOCAM)/2
#> T300_vs_T150                              (G_T300_MI + G_T300_MINOCAM)/2 - (G_T150_MI + G_T150_MINOCAM)/2
#> T150_vs_T0_at_MI                                                                      G_T150_MI - G_T0_MI
#> T150_vs_T0_at_MINOCAM                                                       G_T150_MINOCAM - G_T0_MINOCAM
#> T300_vs_T0_at_MI                                                                      G_T300_MI - G_T0_MI
#> T300_vs_T0_at_MINOCAM                                                       G_T300_MINOCAM - G_T0_MINOCAM
#> T300_vs_T150_at_MI                                                                  G_T300_MI - G_T150_MI
#> T300_vs_T150_at_MINOCAM                                                   G_T300_MINOCAM - G_T150_MINOCAM
#> interaction_T150_vs_T0_at_MINOCAM_vs_MI           (G_T150_MINOCAM - G_T0_MINOCAM) - (G_T150_MI - G_T0_MI)
#> interaction_T300_vs_T0_at_MINOCAM_vs_MI           (G_T300_MINOCAM - G_T0_MINOCAM) - (G_T300_MI - G_T0_MI)
#> interaction_T300_vs_T150_at_MINOCAM_vs_MI     (G_T300_MINOCAM - G_T150_MINOCAM) - (G_T300_MI - G_T150_MI)
```
