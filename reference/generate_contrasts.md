# All contrasts of a two-factor design

Combines
[`main_effect_contrasts()`](https://fgcz.github.io/prolfqua/reference/main_effect_contrasts.md),
[`level_specific_contrasts()`](https://fgcz.github.io/prolfqua/reference/level_specific_contrasts.md)
and, optionally,
[`interaction_contrasts()`](https://fgcz.github.io/prolfqua/reference/interaction_contrasts.md).

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

A data.frame with the columns `ContrastName` and `Contrast` (the
contrast string), one row per contrast.

## See also

Other contrasts:
[`annotation_add_contrasts()`](https://fgcz.github.io/prolfqua/reference/annotation_add_contrasts.md),
[`generate_contrasts_for_factor()`](https://fgcz.github.io/prolfqua/reference/generate_contrasts_for_factor.md),
[`interaction_contrasts()`](https://fgcz.github.io/prolfqua/reference/interaction_contrasts.md),
[`level_specific_contrasts()`](https://fgcz.github.io/prolfqua/reference/level_specific_contrasts.md),
[`main_effect_contrasts()`](https://fgcz.github.io/prolfqua/reference/main_effect_contrasts.md)

## Examples

``` r
primary_levels <- c("KO", "WT")
secondary_levels <- c("T0", "T150", "T300")
generate_contrasts(primary_levels, secondary_levels)
#>                                                              ContrastName
#> WT_vs_KO                                                         WT_vs_KO
#> WT_vs_KO_at_T0                                             WT_vs_KO_at_T0
#> WT_vs_KO_at_T150                                         WT_vs_KO_at_T150
#> WT_vs_KO_at_T300                                         WT_vs_KO_at_T300
#> interaction_WT_vs_KO_at_T150_vs_T0     interaction_WT_vs_KO_at_T150_vs_T0
#> interaction_WT_vs_KO_at_T300_vs_T0     interaction_WT_vs_KO_at_T300_vs_T0
#> interaction_WT_vs_KO_at_T300_vs_T150 interaction_WT_vs_KO_at_T300_vs_T150
#>                                                                                                       Contrast
#> WT_vs_KO                             (G_WT_T0 + G_WT_T150 + G_WT_T300)/3 - (G_KO_T0 + G_KO_T150 + G_KO_T300)/3
#> WT_vs_KO_at_T0                                                                               G_WT_T0 - G_KO_T0
#> WT_vs_KO_at_T150                                                                         G_WT_T150 - G_KO_T150
#> WT_vs_KO_at_T300                                                                         G_WT_T300 - G_KO_T300
#> interaction_WT_vs_KO_at_T150_vs_T0                               (G_WT_T150 - G_KO_T150) - (G_WT_T0 - G_KO_T0)
#> interaction_WT_vs_KO_at_T300_vs_T0                               (G_WT_T300 - G_KO_T300) - (G_WT_T0 - G_KO_T0)
#> interaction_WT_vs_KO_at_T300_vs_T150                         (G_WT_T300 - G_KO_T300) - (G_WT_T150 - G_KO_T150)
generate_contrasts(secondary_levels, primary_levels)
#>                                                              ContrastName
#> T150_vs_T0                                                     T150_vs_T0
#> T300_vs_T0                                                     T300_vs_T0
#> T300_vs_T150                                                 T300_vs_T150
#> T150_vs_T0_at_KO                                         T150_vs_T0_at_KO
#> T150_vs_T0_at_WT                                         T150_vs_T0_at_WT
#> T300_vs_T0_at_KO                                         T300_vs_T0_at_KO
#> T300_vs_T0_at_WT                                         T300_vs_T0_at_WT
#> T300_vs_T150_at_KO                                     T300_vs_T150_at_KO
#> T300_vs_T150_at_WT                                     T300_vs_T150_at_WT
#> interaction_T150_vs_T0_at_WT_vs_KO     interaction_T150_vs_T0_at_WT_vs_KO
#> interaction_T300_vs_T0_at_WT_vs_KO     interaction_T300_vs_T0_at_WT_vs_KO
#> interaction_T300_vs_T150_at_WT_vs_KO interaction_T300_vs_T150_at_WT_vs_KO
#>                                                                                   Contrast
#> T150_vs_T0                               (G_T150_KO + G_T150_WT)/2 - (G_T0_KO + G_T0_WT)/2
#> T300_vs_T0                               (G_T300_KO + G_T300_WT)/2 - (G_T0_KO + G_T0_WT)/2
#> T300_vs_T150                         (G_T300_KO + G_T300_WT)/2 - (G_T150_KO + G_T150_WT)/2
#> T150_vs_T0_at_KO                                                       G_T150_KO - G_T0_KO
#> T150_vs_T0_at_WT                                                       G_T150_WT - G_T0_WT
#> T300_vs_T0_at_KO                                                       G_T300_KO - G_T0_KO
#> T300_vs_T0_at_WT                                                       G_T300_WT - G_T0_WT
#> T300_vs_T150_at_KO                                                   G_T300_KO - G_T150_KO
#> T300_vs_T150_at_WT                                                   G_T300_WT - G_T150_WT
#> interaction_T150_vs_T0_at_WT_vs_KO           (G_T150_WT - G_T0_WT) - (G_T150_KO - G_T0_KO)
#> interaction_T300_vs_T0_at_WT_vs_KO           (G_T300_WT - G_T0_WT) - (G_T300_KO - G_T0_KO)
#> interaction_T300_vs_T150_at_WT_vs_KO     (G_T300_WT - G_T150_WT) - (G_T300_KO - G_T150_KO)
```
