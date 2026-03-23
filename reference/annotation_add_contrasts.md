# DRY function: process and export annotated contrasts

DRY function: process and export annotated contrasts

## Usage

``` r
annotation_add_contrasts(
  df,
  primary_col,
  secondary_col,
  prefix = "G_",
  dataset_id = "dataset",
  decreasing = FALSE,
  interactions = TRUE
)
```

## Arguments

- df:

  data.frame with annotation data

- primary_col:

  character, name of the primary factor column

- secondary_col:

  character, name of the secondary factor column

- prefix:

  character, prefix for output naming (default "G\_")

- dataset_id:

  character, dataset identifier for output naming

- decreasing:

  logical, if TRUE sort factor levels in decreasing order

- interactions:

  logical, if TRUE include interaction contrasts

## See also

Other contrasts:
[`generate_contrasts()`](https://wolski.github.io/prolfqua/reference/generate_contrasts.md),
[`generate_contrasts_for_factor()`](https://wolski.github.io/prolfqua/reference/generate_contrasts_for_factor.md),
[`interaction_contrasts()`](https://wolski.github.io/prolfqua/reference/interaction_contrasts.md),
[`level_specific_contrasts()`](https://wolski.github.io/prolfqua/reference/level_specific_contrasts.md),
[`main_effect_contrasts()`](https://wolski.github.io/prolfqua/reference/main_effect_contrasts.md)

## Examples

``` r

annotation_add_contrasts(prolfqua::x5463yzwer453bbb, "factor_A", "factor_B", "primary")
#> $annot
#> # A tibble: 36 × 5
#>    `Relative Path`                             Name  Group ContrastName Contrast
#>  * <chr>                                       <chr> <chr> <chr>        <chr>   
#>  1 20250429_019_C38312-12_S935004_MI_150_6_S1… MI_1… MI_T… MINOCA_vs_MI ( (G_MI…
#>  2 20250429_003_C38312-21_S935013_MINOCA_0_3_… MINO… MINO… MINOCA_vs_M… G_MINOC…
#>  3 20250429_015_C38312-35_S935027_MINOCA_300_… MINO… MINO… MINOCA_vs_M… G_MINOC…
#>  4 20250429_016_C38312-11_S935003_MI_150_5_S1… MI_1… MI_T… MINOCA_vs_M… G_MINOC…
#>  5 20250429_014_C38312-7_S934999_MI_150_1_S1-… MI_1… MI_T… interaction… (G_MINO…
#>  6 20250429_035_C38312-3_S934995_MI_0_3_S1-C1… MI_0… MI_T0 interaction… (G_MINO…
#>  7 20250429_022_C38312-25_S935017_MINOCA_150_… MINO… MINO… interaction… (G_MINO…
#>  8 20250429_043_C38312-14_S935006_MI_300_2_S1… MI_3… MI_T… NA           NA      
#>  9 20250429_004_C38312-29_S935021_MINOCA_150_… MINO… MINO… NA           NA      
#> 10 20250429_034_C38312-15_S935007_MI_300_3_S1… MI_3… MI_T… NA           NA      
#> # ℹ 26 more rows
#> 
#> $name
#> [1] "DEA_primary_dataset.csv"
#> 
annotation_add_contrasts(prolfqua::x5463yzwer453bbb, "factor_B", "factor_A", "secondary")
#> $annot
#> # A tibble: 36 × 5
#>    `Relative Path`                             Name  Group ContrastName Contrast
#>  * <chr>                                       <chr> <chr> <chr>        <chr>   
#>  1 20250429_019_C38312-12_S935004_MI_150_6_S1… MI_1… T150… T150_vs_T0   ( (G_T1…
#>  2 20250429_003_C38312-21_S935013_MINOCA_0_3_… MINO… T0_M… T300_vs_T0   ( (G_T3…
#>  3 20250429_015_C38312-35_S935027_MINOCA_300_… MINO… T300… T300_vs_T150 ( (G_T3…
#>  4 20250429_016_C38312-11_S935003_MI_150_5_S1… MI_1… T150… T150_vs_T0_… G_T150_…
#>  5 20250429_014_C38312-7_S934999_MI_150_1_S1-… MI_1… T150… T150_vs_T0_… G_T150_…
#>  6 20250429_035_C38312-3_S934995_MI_0_3_S1-C1… MI_0… T0_MI T300_vs_T0_… G_T300_…
#>  7 20250429_022_C38312-25_S935017_MINOCA_150_… MINO… T150… T300_vs_T0_… G_T300_…
#>  8 20250429_043_C38312-14_S935006_MI_300_2_S1… MI_3… T300… T300_vs_T15… G_T300_…
#>  9 20250429_004_C38312-29_S935021_MINOCA_150_… MINO… T150… T300_vs_T15… G_T300_…
#> 10 20250429_034_C38312-15_S935007_MI_300_3_S1… MI_3… T300… interaction… (G_T150…
#> # ℹ 26 more rows
#> 
#> $name
#> [1] "DEA_secondary_dataset.csv"
#> 
```
