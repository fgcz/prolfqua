# Add the contrasts of a two-factor design to a sample annotation

Unites the two factor columns into a `Group` column and appends the
contrasts of
[`generate_contrasts()`](https://fgcz.github.io/prolfqua/reference/generate_contrasts.md)
as the columns `ContrastName` and `Contrast`, in the layout of a
prolfquapp dataset annotation. The shorter side is padded with `NA`
rows: with fewer contrasts than samples the contrast columns, with more
contrasts than samples the sample columns, so the extra rows carry only
a contrast. Nothing is written to disk.

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

## Value

A list with `annot`, the annotation with the `Group`, `ContrastName` and
`Contrast` columns (the two factor columns are dropped), and `name`, a
suggested file name `DEA_<prefix>_<dataset_id>.csv`.

## See also

Other contrasts:
[`generate_contrasts()`](https://fgcz.github.io/prolfqua/reference/generate_contrasts.md),
[`generate_contrasts_for_factor()`](https://fgcz.github.io/prolfqua/reference/generate_contrasts_for_factor.md),
[`interaction_contrasts()`](https://fgcz.github.io/prolfqua/reference/interaction_contrasts.md),
[`level_specific_contrasts()`](https://fgcz.github.io/prolfqua/reference/level_specific_contrasts.md),
[`main_effect_contrasts()`](https://fgcz.github.io/prolfqua/reference/main_effect_contrasts.md)

## Examples

``` r
# 2 genotypes x 3 time points x 6 replicates
design <- expand.grid(
  rep = 1:6, factor_B = c("T0", "T150", "T300"), factor_A = c("WT", "KO"),
  stringsAsFactors = FALSE
)
annotation <- data.frame(
  Name = paste(design$factor_A, design$factor_B, design$rep, sep = "_"),
  factor_A = design$factor_A,
  factor_B = design$factor_B
)
annotation_add_contrasts(annotation, "factor_A", "factor_B", "primary")
#> $annot
#>         Name   Group                         ContrastName
#> 1    WT_T0_1   WT_T0                             WT_vs_KO
#> 2    WT_T0_2   WT_T0                       WT_vs_KO_at_T0
#> 3    WT_T0_3   WT_T0                     WT_vs_KO_at_T150
#> 4    WT_T0_4   WT_T0                     WT_vs_KO_at_T300
#> 5    WT_T0_5   WT_T0   interaction_WT_vs_KO_at_T150_vs_T0
#> 6    WT_T0_6   WT_T0   interaction_WT_vs_KO_at_T300_vs_T0
#> 7  WT_T150_1 WT_T150 interaction_WT_vs_KO_at_T300_vs_T150
#> 8  WT_T150_2 WT_T150                                 <NA>
#> 9  WT_T150_3 WT_T150                                 <NA>
#> 10 WT_T150_4 WT_T150                                 <NA>
#> 11 WT_T150_5 WT_T150                                 <NA>
#> 12 WT_T150_6 WT_T150                                 <NA>
#> 13 WT_T300_1 WT_T300                                 <NA>
#> 14 WT_T300_2 WT_T300                                 <NA>
#> 15 WT_T300_3 WT_T300                                 <NA>
#> 16 WT_T300_4 WT_T300                                 <NA>
#> 17 WT_T300_5 WT_T300                                 <NA>
#> 18 WT_T300_6 WT_T300                                 <NA>
#> 19   KO_T0_1   KO_T0                                 <NA>
#> 20   KO_T0_2   KO_T0                                 <NA>
#> 21   KO_T0_3   KO_T0                                 <NA>
#> 22   KO_T0_4   KO_T0                                 <NA>
#> 23   KO_T0_5   KO_T0                                 <NA>
#> 24   KO_T0_6   KO_T0                                 <NA>
#> 25 KO_T150_1 KO_T150                                 <NA>
#> 26 KO_T150_2 KO_T150                                 <NA>
#> 27 KO_T150_3 KO_T150                                 <NA>
#> 28 KO_T150_4 KO_T150                                 <NA>
#> 29 KO_T150_5 KO_T150                                 <NA>
#> 30 KO_T150_6 KO_T150                                 <NA>
#> 31 KO_T300_1 KO_T300                                 <NA>
#> 32 KO_T300_2 KO_T300                                 <NA>
#> 33 KO_T300_3 KO_T300                                 <NA>
#> 34 KO_T300_4 KO_T300                                 <NA>
#> 35 KO_T300_5 KO_T300                                 <NA>
#> 36 KO_T300_6 KO_T300                                 <NA>
#>                                                                     Contrast
#> 1  (G_WT_T0 + G_WT_T150 + G_WT_T300)/3 - (G_KO_T0 + G_KO_T150 + G_KO_T300)/3
#> 2                                                          G_WT_T0 - G_KO_T0
#> 3                                                      G_WT_T150 - G_KO_T150
#> 4                                                      G_WT_T300 - G_KO_T300
#> 5                              (G_WT_T150 - G_KO_T150) - (G_WT_T0 - G_KO_T0)
#> 6                              (G_WT_T300 - G_KO_T300) - (G_WT_T0 - G_KO_T0)
#> 7                          (G_WT_T300 - G_KO_T300) - (G_WT_T150 - G_KO_T150)
#> 8                                                                       <NA>
#> 9                                                                       <NA>
#> 10                                                                      <NA>
#> 11                                                                      <NA>
#> 12                                                                      <NA>
#> 13                                                                      <NA>
#> 14                                                                      <NA>
#> 15                                                                      <NA>
#> 16                                                                      <NA>
#> 17                                                                      <NA>
#> 18                                                                      <NA>
#> 19                                                                      <NA>
#> 20                                                                      <NA>
#> 21                                                                      <NA>
#> 22                                                                      <NA>
#> 23                                                                      <NA>
#> 24                                                                      <NA>
#> 25                                                                      <NA>
#> 26                                                                      <NA>
#> 27                                                                      <NA>
#> 28                                                                      <NA>
#> 29                                                                      <NA>
#> 30                                                                      <NA>
#> 31                                                                      <NA>
#> 32                                                                      <NA>
#> 33                                                                      <NA>
#> 34                                                                      <NA>
#> 35                                                                      <NA>
#> 36                                                                      <NA>
#> 
#> $name
#> [1] "DEA_primary_dataset.csv"
#> 
annotation_add_contrasts(annotation, "factor_B", "factor_A", "secondary")
#> $annot
#>         Name   Group                         ContrastName
#> 1    WT_T0_1   T0_WT                           T150_vs_T0
#> 2    WT_T0_2   T0_WT                           T300_vs_T0
#> 3    WT_T0_3   T0_WT                         T300_vs_T150
#> 4    WT_T0_4   T0_WT                     T150_vs_T0_at_KO
#> 5    WT_T0_5   T0_WT                     T150_vs_T0_at_WT
#> 6    WT_T0_6   T0_WT                     T300_vs_T0_at_KO
#> 7  WT_T150_1 T150_WT                     T300_vs_T0_at_WT
#> 8  WT_T150_2 T150_WT                   T300_vs_T150_at_KO
#> 9  WT_T150_3 T150_WT                   T300_vs_T150_at_WT
#> 10 WT_T150_4 T150_WT   interaction_T150_vs_T0_at_WT_vs_KO
#> 11 WT_T150_5 T150_WT   interaction_T300_vs_T0_at_WT_vs_KO
#> 12 WT_T150_6 T150_WT interaction_T300_vs_T150_at_WT_vs_KO
#> 13 WT_T300_1 T300_WT                                 <NA>
#> 14 WT_T300_2 T300_WT                                 <NA>
#> 15 WT_T300_3 T300_WT                                 <NA>
#> 16 WT_T300_4 T300_WT                                 <NA>
#> 17 WT_T300_5 T300_WT                                 <NA>
#> 18 WT_T300_6 T300_WT                                 <NA>
#> 19   KO_T0_1   T0_KO                                 <NA>
#> 20   KO_T0_2   T0_KO                                 <NA>
#> 21   KO_T0_3   T0_KO                                 <NA>
#> 22   KO_T0_4   T0_KO                                 <NA>
#> 23   KO_T0_5   T0_KO                                 <NA>
#> 24   KO_T0_6   T0_KO                                 <NA>
#> 25 KO_T150_1 T150_KO                                 <NA>
#> 26 KO_T150_2 T150_KO                                 <NA>
#> 27 KO_T150_3 T150_KO                                 <NA>
#> 28 KO_T150_4 T150_KO                                 <NA>
#> 29 KO_T150_5 T150_KO                                 <NA>
#> 30 KO_T150_6 T150_KO                                 <NA>
#> 31 KO_T300_1 T300_KO                                 <NA>
#> 32 KO_T300_2 T300_KO                                 <NA>
#> 33 KO_T300_3 T300_KO                                 <NA>
#> 34 KO_T300_4 T300_KO                                 <NA>
#> 35 KO_T300_5 T300_KO                                 <NA>
#> 36 KO_T300_6 T300_KO                                 <NA>
#>                                                 Contrast
#> 1      (G_T150_KO + G_T150_WT)/2 - (G_T0_KO + G_T0_WT)/2
#> 2      (G_T300_KO + G_T300_WT)/2 - (G_T0_KO + G_T0_WT)/2
#> 3  (G_T300_KO + G_T300_WT)/2 - (G_T150_KO + G_T150_WT)/2
#> 4                                    G_T150_KO - G_T0_KO
#> 5                                    G_T150_WT - G_T0_WT
#> 6                                    G_T300_KO - G_T0_KO
#> 7                                    G_T300_WT - G_T0_WT
#> 8                                  G_T300_KO - G_T150_KO
#> 9                                  G_T300_WT - G_T150_WT
#> 10         (G_T150_WT - G_T0_WT) - (G_T150_KO - G_T0_KO)
#> 11         (G_T300_WT - G_T0_WT) - (G_T300_KO - G_T0_KO)
#> 12     (G_T300_WT - G_T150_WT) - (G_T300_KO - G_T150_KO)
#> 13                                                  <NA>
#> 14                                                  <NA>
#> 15                                                  <NA>
#> 16                                                  <NA>
#> 17                                                  <NA>
#> 18                                                  <NA>
#> 19                                                  <NA>
#> 20                                                  <NA>
#> 21                                                  <NA>
#> 22                                                  <NA>
#> 23                                                  <NA>
#> 24                                                  <NA>
#> 25                                                  <NA>
#> 26                                                  <NA>
#> 27                                                  <NA>
#> 28                                                  <NA>
#> 29                                                  <NA>
#> 30                                                  <NA>
#> 31                                                  <NA>
#> 32                                                  <NA>
#> 33                                                  <NA>
#> 34                                                  <NA>
#> 35                                                  <NA>
#> 36                                                  <NA>
#> 
#> $name
#> [1] "DEA_secondary_dataset.csv"
#> 
```
