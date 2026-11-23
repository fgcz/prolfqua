# Pairwise contrasts of one factor

One contrast for every pair of levels: the later level minus the earlier
one.

## Usage

``` r
generate_contrasts_for_factor(levels)
```

## Arguments

- levels:

  character vector of factor levels

## Value

A named list of contrast strings `"<high> - <low>"`, named
`<high>_vs_<low>`.

## See also

Other contrasts:
[`annotation_add_contrasts()`](https://fgcz.github.io/prolfqua/reference/annotation_add_contrasts.md),
[`generate_contrasts()`](https://fgcz.github.io/prolfqua/reference/generate_contrasts.md),
[`interaction_contrasts()`](https://fgcz.github.io/prolfqua/reference/interaction_contrasts.md),
[`level_specific_contrasts()`](https://fgcz.github.io/prolfqua/reference/level_specific_contrasts.md),
[`main_effect_contrasts()`](https://fgcz.github.io/prolfqua/reference/main_effect_contrasts.md)

## Examples

``` r
primary_levels <- c("KO", "WT")
secondary_levels <- c("T0", "T150", "T300")
x2 <- level_specific_contrasts(primary_levels, secondary_levels)
x3 <- level_specific_contrasts(secondary_levels, primary_levels)
generate_contrasts_for_factor(names(x2))
#> $WT_vs_KO_at_T150_vs_WT_vs_KO_at_T0
#> [1] "WT_vs_KO_at_T150 - WT_vs_KO_at_T0"
#> 
#> $WT_vs_KO_at_T300_vs_WT_vs_KO_at_T0
#> [1] "WT_vs_KO_at_T300 - WT_vs_KO_at_T0"
#> 
#> $WT_vs_KO_at_T300_vs_WT_vs_KO_at_T150
#> [1] "WT_vs_KO_at_T300 - WT_vs_KO_at_T150"
#> 
```
