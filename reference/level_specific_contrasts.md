# Level-specific contrasts of a two-factor design

For every pair of primary levels, their difference within each secondary
level.

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

A named list of contrast strings, one per pair of primary levels and
secondary level, named `<high>_vs_<low>_at_<secondary>`.

## See also

Other contrasts:
[`annotation_add_contrasts()`](https://fgcz.github.io/prolfqua/reference/annotation_add_contrasts.md),
[`generate_contrasts()`](https://fgcz.github.io/prolfqua/reference/generate_contrasts.md),
[`generate_contrasts_for_factor()`](https://fgcz.github.io/prolfqua/reference/generate_contrasts_for_factor.md),
[`interaction_contrasts()`](https://fgcz.github.io/prolfqua/reference/interaction_contrasts.md),
[`main_effect_contrasts()`](https://fgcz.github.io/prolfqua/reference/main_effect_contrasts.md)

## Examples

``` r
primary_levels <- c("KO", "WT")
secondary_levels <- c("T0", "T150", "T300")
x2 <- level_specific_contrasts(primary_levels, secondary_levels)
x3 <- level_specific_contrasts(secondary_levels, primary_levels)
```
