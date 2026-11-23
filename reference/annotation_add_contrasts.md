# Add the contrasts of a two-factor design to a sample annotation

Unites the two factor columns into a `Group` column and appends the
contrasts of
[`generate_contrasts()`](https://fgcz.github.io/prolfqua/reference/generate_contrasts.md)
as the columns `ContrastName` and `Contrast`, padded with `NA` to the
number of samples, in the layout of a prolfquapp dataset annotation.
Nothing is written to disk.

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
#>                                           Name   Group
#> WT_vs_KO                               WT_T0_1   WT_T0
#> WT_vs_KO_at_T0                         WT_T0_2   WT_T0
#> WT_vs_KO_at_T150                       WT_T0_3   WT_T0
#> WT_vs_KO_at_T300                       WT_T0_4   WT_T0
#> interaction_WT_vs_KO_at_T150_vs_T0     WT_T0_5   WT_T0
#> interaction_WT_vs_KO_at_T300_vs_T0     WT_T0_6   WT_T0
#> interaction_WT_vs_KO_at_T300_vs_T150 WT_T150_1 WT_T150
#> ...8                                 WT_T150_2 WT_T150
#> ...9                                 WT_T150_3 WT_T150
#> ...10                                WT_T150_4 WT_T150
#> ...11                                WT_T150_5 WT_T150
#> ...12                                WT_T150_6 WT_T150
#> ...13                                WT_T300_1 WT_T300
#> ...14                                WT_T300_2 WT_T300
#> ...15                                WT_T300_3 WT_T300
#> ...16                                WT_T300_4 WT_T300
#> ...17                                WT_T300_5 WT_T300
#> ...18                                WT_T300_6 WT_T300
#> ...19                                  KO_T0_1   KO_T0
#> ...20                                  KO_T0_2   KO_T0
#> ...21                                  KO_T0_3   KO_T0
#> ...22                                  KO_T0_4   KO_T0
#> ...23                                  KO_T0_5   KO_T0
#> ...24                                  KO_T0_6   KO_T0
#> ...25                                KO_T150_1 KO_T150
#> ...26                                KO_T150_2 KO_T150
#> ...27                                KO_T150_3 KO_T150
#> ...28                                KO_T150_4 KO_T150
#> ...29                                KO_T150_5 KO_T150
#> ...30                                KO_T150_6 KO_T150
#> ...31                                KO_T300_1 KO_T300
#> ...32                                KO_T300_2 KO_T300
#> ...33                                KO_T300_3 KO_T300
#> ...34                                KO_T300_4 KO_T300
#> ...35                                KO_T300_5 KO_T300
#> ...36                                KO_T300_6 KO_T300
#>                                                              ContrastName
#> WT_vs_KO                                                         WT_vs_KO
#> WT_vs_KO_at_T0                                             WT_vs_KO_at_T0
#> WT_vs_KO_at_T150                                         WT_vs_KO_at_T150
#> WT_vs_KO_at_T300                                         WT_vs_KO_at_T300
#> interaction_WT_vs_KO_at_T150_vs_T0     interaction_WT_vs_KO_at_T150_vs_T0
#> interaction_WT_vs_KO_at_T300_vs_T0     interaction_WT_vs_KO_at_T300_vs_T0
#> interaction_WT_vs_KO_at_T300_vs_T150 interaction_WT_vs_KO_at_T300_vs_T150
#> ...8                                                                 <NA>
#> ...9                                                                 <NA>
#> ...10                                                                <NA>
#> ...11                                                                <NA>
#> ...12                                                                <NA>
#> ...13                                                                <NA>
#> ...14                                                                <NA>
#> ...15                                                                <NA>
#> ...16                                                                <NA>
#> ...17                                                                <NA>
#> ...18                                                                <NA>
#> ...19                                                                <NA>
#> ...20                                                                <NA>
#> ...21                                                                <NA>
#> ...22                                                                <NA>
#> ...23                                                                <NA>
#> ...24                                                                <NA>
#> ...25                                                                <NA>
#> ...26                                                                <NA>
#> ...27                                                                <NA>
#> ...28                                                                <NA>
#> ...29                                                                <NA>
#> ...30                                                                <NA>
#> ...31                                                                <NA>
#> ...32                                                                <NA>
#> ...33                                                                <NA>
#> ...34                                                                <NA>
#> ...35                                                                <NA>
#> ...36                                                                <NA>
#>                                                                                                       Contrast
#> WT_vs_KO                             (G_WT_T0 + G_WT_T150 + G_WT_T300)/3 - (G_KO_T0 + G_KO_T150 + G_KO_T300)/3
#> WT_vs_KO_at_T0                                                                               G_WT_T0 - G_KO_T0
#> WT_vs_KO_at_T150                                                                         G_WT_T150 - G_KO_T150
#> WT_vs_KO_at_T300                                                                         G_WT_T300 - G_KO_T300
#> interaction_WT_vs_KO_at_T150_vs_T0                               (G_WT_T150 - G_KO_T150) - (G_WT_T0 - G_KO_T0)
#> interaction_WT_vs_KO_at_T300_vs_T0                               (G_WT_T300 - G_KO_T300) - (G_WT_T0 - G_KO_T0)
#> interaction_WT_vs_KO_at_T300_vs_T150                         (G_WT_T300 - G_KO_T300) - (G_WT_T150 - G_KO_T150)
#> ...8                                                                                                      <NA>
#> ...9                                                                                                      <NA>
#> ...10                                                                                                     <NA>
#> ...11                                                                                                     <NA>
#> ...12                                                                                                     <NA>
#> ...13                                                                                                     <NA>
#> ...14                                                                                                     <NA>
#> ...15                                                                                                     <NA>
#> ...16                                                                                                     <NA>
#> ...17                                                                                                     <NA>
#> ...18                                                                                                     <NA>
#> ...19                                                                                                     <NA>
#> ...20                                                                                                     <NA>
#> ...21                                                                                                     <NA>
#> ...22                                                                                                     <NA>
#> ...23                                                                                                     <NA>
#> ...24                                                                                                     <NA>
#> ...25                                                                                                     <NA>
#> ...26                                                                                                     <NA>
#> ...27                                                                                                     <NA>
#> ...28                                                                                                     <NA>
#> ...29                                                                                                     <NA>
#> ...30                                                                                                     <NA>
#> ...31                                                                                                     <NA>
#> ...32                                                                                                     <NA>
#> ...33                                                                                                     <NA>
#> ...34                                                                                                     <NA>
#> ...35                                                                                                     <NA>
#> ...36                                                                                                     <NA>
#> 
#> $name
#> [1] "DEA_primary_dataset.csv"
#> 
annotation_add_contrasts(annotation, "factor_B", "factor_A", "secondary")
#> $annot
#>                                           Name   Group
#> T150_vs_T0                             WT_T0_1   T0_WT
#> T300_vs_T0                             WT_T0_2   T0_WT
#> T300_vs_T150                           WT_T0_3   T0_WT
#> T150_vs_T0_at_KO                       WT_T0_4   T0_WT
#> T150_vs_T0_at_WT                       WT_T0_5   T0_WT
#> T300_vs_T0_at_KO                       WT_T0_6   T0_WT
#> T300_vs_T0_at_WT                     WT_T150_1 T150_WT
#> T300_vs_T150_at_KO                   WT_T150_2 T150_WT
#> T300_vs_T150_at_WT                   WT_T150_3 T150_WT
#> interaction_T150_vs_T0_at_WT_vs_KO   WT_T150_4 T150_WT
#> interaction_T300_vs_T0_at_WT_vs_KO   WT_T150_5 T150_WT
#> interaction_T300_vs_T150_at_WT_vs_KO WT_T150_6 T150_WT
#> ...13                                WT_T300_1 T300_WT
#> ...14                                WT_T300_2 T300_WT
#> ...15                                WT_T300_3 T300_WT
#> ...16                                WT_T300_4 T300_WT
#> ...17                                WT_T300_5 T300_WT
#> ...18                                WT_T300_6 T300_WT
#> ...19                                  KO_T0_1   T0_KO
#> ...20                                  KO_T0_2   T0_KO
#> ...21                                  KO_T0_3   T0_KO
#> ...22                                  KO_T0_4   T0_KO
#> ...23                                  KO_T0_5   T0_KO
#> ...24                                  KO_T0_6   T0_KO
#> ...25                                KO_T150_1 T150_KO
#> ...26                                KO_T150_2 T150_KO
#> ...27                                KO_T150_3 T150_KO
#> ...28                                KO_T150_4 T150_KO
#> ...29                                KO_T150_5 T150_KO
#> ...30                                KO_T150_6 T150_KO
#> ...31                                KO_T300_1 T300_KO
#> ...32                                KO_T300_2 T300_KO
#> ...33                                KO_T300_3 T300_KO
#> ...34                                KO_T300_4 T300_KO
#> ...35                                KO_T300_5 T300_KO
#> ...36                                KO_T300_6 T300_KO
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
#> ...13                                                                <NA>
#> ...14                                                                <NA>
#> ...15                                                                <NA>
#> ...16                                                                <NA>
#> ...17                                                                <NA>
#> ...18                                                                <NA>
#> ...19                                                                <NA>
#> ...20                                                                <NA>
#> ...21                                                                <NA>
#> ...22                                                                <NA>
#> ...23                                                                <NA>
#> ...24                                                                <NA>
#> ...25                                                                <NA>
#> ...26                                                                <NA>
#> ...27                                                                <NA>
#> ...28                                                                <NA>
#> ...29                                                                <NA>
#> ...30                                                                <NA>
#> ...31                                                                <NA>
#> ...32                                                                <NA>
#> ...33                                                                <NA>
#> ...34                                                                <NA>
#> ...35                                                                <NA>
#> ...36                                                                <NA>
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
#> ...13                                                                                 <NA>
#> ...14                                                                                 <NA>
#> ...15                                                                                 <NA>
#> ...16                                                                                 <NA>
#> ...17                                                                                 <NA>
#> ...18                                                                                 <NA>
#> ...19                                                                                 <NA>
#> ...20                                                                                 <NA>
#> ...21                                                                                 <NA>
#> ...22                                                                                 <NA>
#> ...23                                                                                 <NA>
#> ...24                                                                                 <NA>
#> ...25                                                                                 <NA>
#> ...26                                                                                 <NA>
#> ...27                                                                                 <NA>
#> ...28                                                                                 <NA>
#> ...29                                                                                 <NA>
#> ...30                                                                                 <NA>
#> ...31                                                                                 <NA>
#> ...32                                                                                 <NA>
#> ...33                                                                                 <NA>
#> ...34                                                                                 <NA>
#> ...35                                                                                 <NA>
#> ...36                                                                                 <NA>
#> 
#> $name
#> [1] "DEA_secondary_dataset.csv"
#> 
```
