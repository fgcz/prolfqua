# Level-specific contrasts (per secondary level)

Level-specific contrasts (per secondary level)

## Usage

``` r
level_specific_contrasts(primary_levels, secondary_levels)
```

## Arguments

- primary_levels:

  character vector of primary factor levels

- secondary_levels:

  character vector of secondary factor levels

## Value

Contrast definitions or contrast results.

## See also

Other contrasts:
[`annotation_add_contrasts()`](https://wolski.github.io/prolfqua/reference/annotation_add_contrasts.md),
[`generate_contrasts()`](https://wolski.github.io/prolfqua/reference/generate_contrasts.md),
[`generate_contrasts_for_factor()`](https://wolski.github.io/prolfqua/reference/generate_contrasts_for_factor.md),
[`interaction_contrasts()`](https://wolski.github.io/prolfqua/reference/interaction_contrasts.md),
[`main_effect_contrasts()`](https://wolski.github.io/prolfqua/reference/main_effect_contrasts.md)

## Examples

``` r
# example code

primary_levels <- c("MI", "MINOCAM")
secondary_levels <- c("T0", "T150", "T300")
x2 <- level_specific_contrasts(primary_levels, secondary_levels)
x3 <- level_specific_contrasts(secondary_levels, primary_levels)
```
