# Single-factor contrasts (pairwise comparisons)

Single-factor contrasts (pairwise comparisons)

## Usage

``` r
generate_contrasts_for_factor(levels)
```

## Arguments

- levels:

  character vector of factor levels

## Value

Contrast definitions or contrast results.

## See also

Other contrasts:
[`annotation_add_contrasts()`](https://wolski.github.io/prolfqua/reference/annotation_add_contrasts.md),
[`generate_contrasts()`](https://wolski.github.io/prolfqua/reference/generate_contrasts.md),
[`interaction_contrasts()`](https://wolski.github.io/prolfqua/reference/interaction_contrasts.md),
[`level_specific_contrasts()`](https://wolski.github.io/prolfqua/reference/level_specific_contrasts.md),
[`main_effect_contrasts()`](https://wolski.github.io/prolfqua/reference/main_effect_contrasts.md)

## Examples

``` r
# example code
primary_levels <- c("MI", "MINOCAM")
secondary_levels <- c("T0", "T150", "T300")
x2 <- level_specific_contrasts(primary_levels, secondary_levels)
x3 <- level_specific_contrasts(secondary_levels, primary_levels)
generate_contrasts_for_factor(names(x2))
#> $MINOCAM_vs_MI_at_T150_vs_MINOCAM_vs_MI_at_T0
#> [1] "MINOCAM_vs_MI_at_T150 - MINOCAM_vs_MI_at_T0"
#> 
#> $MINOCAM_vs_MI_at_T300_vs_MINOCAM_vs_MI_at_T0
#> [1] "MINOCAM_vs_MI_at_T300 - MINOCAM_vs_MI_at_T0"
#> 
#> $MINOCAM_vs_MI_at_T300_vs_MINOCAM_vs_MI_at_T150
#> [1] "MINOCAM_vs_MI_at_T300 - MINOCAM_vs_MI_at_T150"
#> 
```
