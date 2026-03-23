# Contrast Facades with Parallel Designs

## Purpose

[`build_contrast_analysis()`](https://wolski.github.io/prolfqua/reference/build_contrast_analysis.md)
provides a common front-end for several contrast backends:

- `lm`
- `limma`
- `lmer`
- `ropeca`
- `lm_missing`
- `lm_impute`
- `deqms`
- `firth`

All of them expose the same basic interface: `get_contrasts()`,
`get_Plotter()`, and `to_wide()`. All examples below return
protein-level contrasts. The important difference is the required input
level:

- `lm`, `limma`, `lm_missing`, `lm_impute`, and `deqms` require
  aggregated protein-level data
- `lmer` and `ropeca` require lower-level measurements such as peptides
  nested within proteins, but still report protein-level contrasts
- `firth` can be used with either aggregated protein-level data or
  nested peptide-level data and still reports protein-level contrasts

This vignette starts from one simulated peptide-level experiment,
aggregates it to protein level, and then demonstrates both families of
facades separately.

## Simulate one experiment

``` r
istar <- sim_lfq_data_peptide_config(Nprot = 80, seed = 42)
istar$config <- old2new(istar$config)

lfq_peptide <- LFQData$new(istar$data, istar$config)
lfq_peptide <- lfq_peptide$get_Transformer()$log2()$lfq

aggregator <- LFQDataAggregator$new(lfq_peptide, "protein")
lfq_protein <- aggregator$medpolish()

lfq_peptide$config$hierarchy_keys()
```

    ## [1] "protein_Id" "peptide_Id"

``` r
lfq_protein$config$hierarchy_keys()
```

    ## [1] "protein_Id"

``` r
lfq_protein$config$nr_children
```

    ## [1] "nr_children_protein_Id"

The aggregation step produces protein-level intensities while keeping
count metadata in the `LFQData` object. The DEqMS facade uses
`lfq_protein$config$nr_children` directly, so no extra count table has
to be passed around.

## Define contrasts

``` r
contrasts <- c(
  "A_vs_Ctrl" = "group_A - group_Ctrl",
  "B_vs_Ctrl" = "group_B - group_Ctrl"
)
```

Two contrasts let us see how each backend handles multiple comparisons
and how FDR correction propagates across contrasts.

## Protein-input facades

The following facades require aggregated input, which in practice means
`lfqdata$subject_Id()` must match `lfqdata$config$hierarchy_keys()`.
`firth` is included here on purpose because it can be fitted directly on
aggregated protein input.

``` r
fa_lm <- build_contrast_analysis(
  lfq_protein,
  "~ group_",
  contrasts,
  method = "lm"
)

fa_limma <- build_contrast_analysis(
  lfq_protein,
  "~ group_",
  contrasts,
  method = "limma"
)

fa_lm_missing <- build_contrast_analysis(
  lfq_protein,
  "~ group_",
  contrasts,
  method = "lm_missing"
)

fa_deqms <- build_contrast_analysis(
  lfq_protein,
  "~ group_",
  contrasts,
  method = "deqms"
)

fa_lm_impute <- build_contrast_analysis(
  lfq_protein,
  "~ group_",
  contrasts,
  method = "lm_impute"
)

fa_firth_protein <- build_contrast_analysis(
  lfq_protein,
  "~ group_",
  contrasts,
  method = "firth"
)
```

Because all protein-input facades share the same interface and report
protein-level contrasts, their outputs can be combined directly.

``` r
results_protein <- bind_rows(
  fa_lm$get_contrasts(),
  fa_limma$get_contrasts(),
  fa_lm_missing$get_contrasts(),
  fa_lm_impute$get_contrasts(),
  fa_deqms$get_contrasts(),
  fa_firth_protein$get_contrasts()
) |>
  dplyr::select(dplyr::any_of(c(
    "facade", "modelName", "protein_Id", "contrast", "avgAbd", "diff", "FDR",
    "statistic", "std.error", "df", "p.value", "conf.low", "conf.high",
    "sigma"
  ))) |>
  dplyr::mutate(
    significant = FDR < 0.1 & abs(diff) > 0.5
  )

results_protein |>
  dplyr::count(facade, name = "n_results")
```

    ## # A tibble: 6 × 2
    ##   facade     n_results
    ##   <chr>          <int>
    ## 1 deqms            157
    ## 2 firth            160
    ## 3 limma            155
    ## 4 lm               157
    ## 5 lm_impute        160
    ## 6 lm_missing       160

For facades that combine several underlying result types, such as
`lm_missing`, the `modelName` column still tells you where individual
rows came from.

``` r
results_protein |>
  dplyr::count(facade, contrast, modelName, name = "n_results")
```

    ## # A tibble: 14 × 4
    ##    facade     contrast  modelName          n_results
    ##    <chr>      <chr>     <chr>                  <int>
    ##  1 deqms      A_vs_Ctrl WaldTest_DEqMS            78
    ##  2 deqms      B_vs_Ctrl WaldTest_DEqMS            79
    ##  3 firth      A_vs_Ctrl WaldTestFirth             80
    ##  4 firth      B_vs_Ctrl WaldTestFirth             80
    ##  5 limma      A_vs_Ctrl limma                     78
    ##  6 limma      B_vs_Ctrl limma                     77
    ##  7 lm         A_vs_Ctrl WaldTest_moderated        78
    ##  8 lm         B_vs_Ctrl WaldTest_moderated        79
    ##  9 lm_impute  A_vs_Ctrl WaldTest_moderated        80
    ## 10 lm_impute  B_vs_Ctrl WaldTest_moderated        80
    ## 11 lm_missing A_vs_Ctrl WaldTest_moderated        78
    ## 12 lm_missing A_vs_Ctrl groupAverage               2
    ## 13 lm_missing B_vs_Ctrl WaldTest_moderated        79
    ## 14 lm_missing B_vs_Ctrl groupAverage               1

## Protein-level volcano comparison

``` r
ggplot(results_protein, aes(x = diff, y = -log10(p.value), color = significant)) +
  geom_point(alpha = 0.6, size = 1.2) +
  facet_grid(contrast ~ facade, scales = "free_y") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey70")) +
  labs(
    x = "log2 fold change",
    y = "-log10(p.value)",
    color = "FDR < 0.1\nand |diff| > 0.5"
  ) +
  theme_minimal(base_size = 12)
```

![Volcano plots for the protein-level facades. Rows are contrasts,
columns are
backends.](ContrastFacades_files/figure-html/volcano_protein-1.png)

Volcano plots for the protein-level facades. Rows are contrasts, columns
are backends.

## Looking at the strongest protein-level hits

``` r
results_protein |>
  dplyr::group_by(facade, contrast) |>
  dplyr::slice_min(order_by = p.value, n = 5, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(facade, contrast, modelName, protein_Id, diff, p.value, FDR)
```

    ## # A tibble: 60 × 7
    ##    facade contrast  modelName      protein_Id    diff  p.value      FDR
    ##    <chr>  <chr>     <chr>          <chr>        <dbl>    <dbl>    <dbl>
    ##  1 deqms  A_vs_Ctrl WaldTest_DEqMS Zci7Jw~7064 -0.718 5.63e-11 4.39e- 9
    ##  2 deqms  A_vs_Ctrl WaldTest_DEqMS 4Y4DYT~0927  0.583 3.43e-10 1.10e- 8
    ##  3 deqms  A_vs_Ctrl WaldTest_DEqMS 6TevMr~7550  0.765 4.24e-10 1.10e- 8
    ##  4 deqms  A_vs_Ctrl WaldTest_DEqMS 0CubNR~0890  0.674 7.29e-10 1.42e- 8
    ##  5 deqms  A_vs_Ctrl WaldTest_DEqMS KVkccD~1805 -0.531 3.00e- 9 4.68e- 8
    ##  6 deqms  B_vs_Ctrl WaldTest_DEqMS fylZqB~3883  0.717 6.57e-13 5.19e-11
    ##  7 deqms  B_vs_Ctrl WaldTest_DEqMS XxJoJB~7286  0.587 2.16e-11 8.51e-10
    ##  8 deqms  B_vs_Ctrl WaldTest_DEqMS f0Cvvj~6658  0.345 2.89e- 9 7.62e- 8
    ##  9 deqms  B_vs_Ctrl WaldTest_DEqMS TR3Ksv~1492 -0.413 1.07e- 8 2.12e- 7
    ## 10 deqms  B_vs_Ctrl WaldTest_DEqMS 4Y4DYT~0927  0.424 2.65e- 8 4.18e- 7
    ## # ℹ 50 more rows

## Proteins that could not be estimated

Every facade has a `get_missing()` method that returns the protein ×
contrast pairs present in the input data but absent from
`get_contrasts()`. This makes it easy to see which proteins each method
fails on and to compare coverage.

``` r
missing_all <- dplyr::bind_rows(
  fa_lm$get_missing() |> dplyr::mutate(facade = "lm"),
  fa_limma$get_missing() |> dplyr::mutate(facade = "limma"),
  fa_lm_missing$get_missing() |> dplyr::mutate(facade = "lm_missing"),
  fa_lm_impute$get_missing() |> dplyr::mutate(facade = "lm_impute"),
  fa_deqms$get_missing() |> dplyr::mutate(facade = "deqms"),
  fa_firth_protein$get_missing() |> dplyr::mutate(facade = "firth")
)

missing_all |>
  dplyr::count(facade, contrast, name = "n_missing") |>
  knitr::kable(caption = "Number of missing protein × contrast pairs per facade")
```

| facade | contrast  | n_missing |
|:-------|:----------|----------:|
| deqms  | A_vs_Ctrl |         2 |
| deqms  | B_vs_Ctrl |         1 |
| limma  | A_vs_Ctrl |         2 |
| limma  | B_vs_Ctrl |         3 |
| lm     | A_vs_Ctrl |         2 |
| lm     | B_vs_Ctrl |         1 |

Number of missing protein × contrast pairs per facade

### Per-sample intensities of the missing proteins

``` r
missing_proteins <- unique(missing_all$protein_Id)

if (length(missing_proteins) > 0) {
  lfq_protein$data |>
    dplyr::filter(protein_Id %in% missing_proteins) |>
    dplyr::select(protein_Id, sampleName,
                  !!rlang::sym(lfq_protein$config$get_response())) |>
    tidyr::pivot_wider(names_from = sampleName,
                       values_from = !!rlang::sym(lfq_protein$config$get_response())) |>
    knitr::kable(digits = 2, caption = "Per-sample intensities of proteins that could not be estimated")
}
```

| protein_Id  | B_V1 | B_V4 | Ctrl_V3 | Ctrl_V4 | B_V2 | B_V3 | Ctrl_V2 | A_V4 | Ctrl_V1 |
|:------------|-----:|-----:|--------:|--------:|-----:|-----:|--------:|-----:|--------:|
| 8mS8sK~0150 | 3.85 | 3.76 |    3.37 |    3.55 |   NA |   NA |      NA |   NA |      NA |
| DTCi0N~0734 |   NA | 4.28 |    4.07 |    4.21 | 4.37 | 4.35 |    4.06 |   NA |      NA |
| OrL0ux~1369 |   NA |   NA |      NA |    3.98 |   NA |   NA |    4.05 | 3.78 |    4.12 |

Per-sample intensities of proteins that could not be estimated

The missing cells (NA) explain why these proteins cannot be estimated —
they lack observations in one or more groups. The `lm_missing` facade
fills in these gaps via group-mean imputation, while `lm_impute` re-fits
the lm model after imputing individual values with the LOD and borrowing
covariance from successful fits. Both should have fewer missing proteins
than plain `lm`.

### Estimates for the missing proteins from `lm_missing` and `lm_impute`

For proteins that plain `lm` could not estimate, the imputation-based
facades can still produce contrast results. The table below shows these
rescued estimates side by side.

``` r
lm_missing_proteins <- fa_lm$get_missing()$protein_Id |> unique()

if (length(lm_missing_proteins) > 0) {
  rescued <- results_protein |>
    dplyr::filter(
      protein_Id %in% lm_missing_proteins,
      facade %in% c("lm_missing", "lm_impute")
    ) |>
    dplyr::arrange(protein_Id, contrast, facade)

  rescued |>
    knitr::kable(
      digits = 3,
      caption = "Contrast estimates from lm_missing and lm_impute for proteins that plain lm could not estimate"
    )
} else {
  cat("All proteins were estimable by plain lm — no rescued estimates to show.")
}
```

| facade     | modelName          | protein_Id  | contrast  | avgAbd |   diff |   FDR | statistic | std.error |    df | p.value | conf.low | conf.high | sigma | significant |
|:-----------|:-------------------|:------------|:----------|-------:|-------:|------:|----------:|----------:|------:|--------:|---------:|----------:|------:|:------------|
| lm_impute  | WaldTest_moderated | 8mS8sK~0150 | A_vs_Ctrl |  3.776 |  0.000 | 1.000 |     0.000 |     0.013 | 4.746 |   1.000 |   -0.154 |     0.154 | 0.059 | FALSE       |
| lm_missing | groupAverage       | 8mS8sK~0150 | A_vs_Ctrl |  3.697 |  0.000 | 1.000 |     0.000 |     0.102 | 2.000 |   1.000 |   -0.437 |     0.437 | 0.102 | FALSE       |
| lm_impute  | WaldTest_moderated | 8mS8sK~0150 | B_vs_Ctrl |  3.776 |  0.001 | 0.922 |     0.119 |     0.013 | 4.746 |   0.910 |   -0.153 |     0.156 | 0.059 | FALSE       |
| lm_missing | WaldTest_moderated | 8mS8sK~0150 | B_vs_Ctrl |  3.632 |  0.339 | 0.009 |     4.446 |     0.102 | 5.746 |   0.005 |    0.150 |     0.527 | 0.076 | FALSE       |
| lm_impute  | WaldTest_moderated | DTCi0N~0734 | A_vs_Ctrl |  3.786 | -0.022 | 0.193 |    -1.776 |     0.013 | 6.746 |   0.121 |   -0.165 |     0.122 | 0.060 | FALSE       |
| lm_missing | groupAverage       | DTCi0N~0734 | A_vs_Ctrl |  3.902 | -0.253 | 0.032 |    -4.417 |     0.057 | 4.000 |   0.012 |   -0.412 |    -0.094 | 0.070 | FALSE       |
| lm_impute  | WaldTest_moderated | DTCi0N~0734 | B_vs_Ctrl |  3.803 |  0.013 | 0.386 |     1.032 |     0.013 | 6.746 |   0.338 |   -0.131 |     0.156 | 0.060 | FALSE       |
| lm_missing | WaldTest_moderated | DTCi0N~0734 | B_vs_Ctrl |  4.224 |  0.222 | 0.006 |     4.194 |     0.057 | 7.746 |   0.003 |    0.072 |     0.372 | 0.065 | FALSE       |
| lm_impute  | WaldTest_moderated | OrL0ux~1369 | A_vs_Ctrl |  3.784 | -0.018 | 0.293 |    -1.481 |     0.013 | 4.746 |   0.202 |   -0.172 |     0.137 | 0.059 | FALSE       |
| lm_missing | WaldTest_moderated | OrL0ux~1369 | A_vs_Ctrl |  3.913 | -0.276 | 0.023 |    -3.752 |     0.084 | 5.757 |   0.010 |   -0.433 |    -0.118 | 0.064 | FALSE       |
| lm_impute  | WaldTest_moderated | OrL0ux~1369 | B_vs_Ctrl |  3.784 | -0.018 | 0.247 |    -1.460 |     0.013 | 4.746 |   0.207 |   -0.172 |     0.137 | 0.059 | FALSE       |
| lm_missing | groupAverage       | OrL0ux~1369 | B_vs_Ctrl |  3.879 | -0.207 | 0.094 |    -3.473 |     0.060 | 2.000 |   0.074 |   -0.463 |     0.049 | 0.073 | FALSE       |

Contrast estimates from lm_missing and lm_impute for proteins that plain
lm could not estimate

## Peptide-input facades

The mixed-effects `lmer` facade and `ropeca` require lower-level
measurements below the analysis subject. The `firth` facade can also
operate directly on peptide-level `LFQData`. `firth` is shown a second
time here on purpose, because it can also be fitted on peptide input.
All three still return protein-level contrasts.

``` r
fa_lmer <- build_contrast_analysis(
  lfq_peptide,
  "~ group_ + (1 | peptide_Id) + (1 | sampleName)",
  contrasts,
  method = "lmer"
)

fa_ropeca <- build_contrast_analysis(
  lfq_peptide,
  "~ group_",
  contrasts,
  method = "ropeca"
)

fa_firth_peptide <- build_contrast_analysis(
  lfq_peptide,
  "~ group_",
  contrasts,
  method = "firth"
)
```

`ropeca` aggregates peptide evidence back to proteins, whereas `lmer`
models the nested peptide structure directly before reporting
protein-level contrasts. Peptide-level `firth` also reports
protein-level contrasts. Proteins with exactly one peptide are fitted
without an added peptide term, while proteins with multiple peptides are
fitted with the lowest hierarchy key added internally.

``` r
results_peptide <- bind_rows(
  fa_lmer$get_contrasts(),
  fa_ropeca$get_contrasts(),
  fa_firth_peptide$get_contrasts()
) |>
  dplyr::select(dplyr::any_of(c(
    "facade", "modelName", "protein_Id", "contrast", "avgAbd", "diff", "FDR",
    "statistic", "std.error", "df", "p.value", "conf.low", "conf.high",
    "sigma"
  ))) |>
  dplyr::mutate(
    significant = FDR < 0.1 & abs(diff) > 0.5
  )

results_peptide |>
  dplyr::count(facade, name = "n_results")
```

    ## # A tibble: 3 × 2
    ##   facade n_results
    ##   <chr>      <int>
    ## 1 firth        160
    ## 2 lmer         102
    ## 3 ropeca       157

``` r
ggplot(results_peptide, aes(x = diff, y = -log10(p.value), color = significant)) +
  geom_point(alpha = 0.6, size = 1.2) +
  facet_grid(contrast ~ facade, scales = "free_y") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey70")) +
  labs(
    x = "log2 fold change",
    y = "-log10(p.value)",
    color = "FDR < 0.1\nand |diff| > 0.5"
  ) +
  theme_minimal(base_size = 12)
```

![Volcano plots for the peptide-level facades. Rows are contrasts,
columns are
backends.](ContrastFacades_files/figure-html/volcano_peptide-1.png)

Volcano plots for the peptide-level facades. Rows are contrasts, columns
are backends.

## Remarks

The facades make it easy to benchmark alternative contrast backends
without rewriting the analysis pipeline:

- the protein-level facades now enforce aggregation before modelling
- the peptide-level facades now explicitly require lower-level hierarchy
  below the analysis subject
- the shared facade API still makes it straightforward to compare
  methods once the data level is chosen consistently

The results are comparable at the API level, but the comparison is only
meaningful when methods using the same biological unit are plotted
together.
