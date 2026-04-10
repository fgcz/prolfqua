# Contrast Facades with Parallel Designs

## Purpose

[`build_contrast_analysis()`](https://wolski.github.io/prolfqua/reference/build_contrast_analysis.md)
provides a common front-end for several contrast backends:

- `lm`
- `limma`
- `limma_impute`
- `limpa` (requires `AggregateLimpa` for DPC-based aggregation with
  standard errors)
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

- `lm`, `limma`, `limma_impute`, `lm_missing`, `lm_impute`, and `deqms`
  require aggregated protein-level data
- `limpa` requires protein-level data from `AggregateLimpa` (which
  provides standard errors and observation counts for vooma precision
  weighting)
- `lmer` and `ropeca` require lower-level measurements such as peptides
  nested within proteins, but still report protein-level contrasts
- `firth` can be used with either aggregated protein-level data or
  nested peptide-level data and still reports protein-level contrasts

This vignette starts from one simulated peptide-level experiment,
aggregates it to protein level, and then demonstrates both families of
facades separately.

## Simulate one experiment

``` r
options(prolfqua.vectorize = TRUE)

istar <- sim_lfq_data_peptide_config(Nprot = 80, seed = 42)


lfq_peptide <- LFQData$new(istar$data, istar$config)
lfq_peptide <- lfq_peptide$get_Transformer()$log2()$lfq

lfq_protein <- lfq_peptide$get_Aggregator()$aggregate()

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
`lfqdata$subject_id()` must match `lfqdata$config$hierarchy_keys()`.
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

fa_limma_impute <- build_contrast_analysis(
  lfq_protein,
  "~ group_",
  contrasts,
  method = "limma_impute"
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

#### limpa facade (DPC-based aggregation)

The `limpa` facade is different from the other protein-level facades: it
requires its own aggregation step via `AggregateLimpa`, which uses
limpa’s Detection Probability Curve (DPC) to aggregate peptides to
proteins while producing per-protein, per-sample standard errors. These
SEs feed into a bivariate vooma variance model for precision weighting.

``` r
lfq_limpa <- AggregateLimpa$new(lfq_peptide, "protein")$aggregate()

fa_limpa <- build_contrast_analysis(
  lfq_limpa,
  "~ group_",
  contrasts,
  method = "limpa"
)
```

Because all protein-input facades share the same interface and report
protein-level contrasts, their outputs can be combined directly.

``` r
# Proteins missing in the baseline lm facade (used to flag rescued proteins)
lm_missing_ids <- fa_lm$get_missing() |>
  dplyr::select(protein_Id, contrast) |>
  dplyr::mutate(rescued = TRUE)

results_limpa <- if (requireNamespace("limpa", quietly = TRUE)) {
  fa_limpa$get_contrasts()
} else {
  data.frame()
}

results_protein <- bind_rows(
  fa_lm$get_contrasts(),
  fa_limma$get_contrasts(),
  fa_limma_impute$get_contrasts(),
  fa_lm_missing$get_contrasts(),
  fa_lm_impute$get_contrasts(),
  fa_deqms$get_contrasts(),
  fa_firth_protein$get_contrasts(),
  results_limpa
) |>
  dplyr::select(dplyr::any_of(c(
    "facade", "modelName", "protein_Id", "contrast", "avgAbd", "diff", "FDR",
    "statistic", "std.error", "df", "p.value", "conf.low", "conf.high",
    "sigma"
  ))) |>
  dplyr::left_join(lm_missing_ids, by = c("protein_Id", "contrast")) |>
  dplyr::mutate(
    rescued = dplyr::coalesce(rescued, FALSE),
    significant = FDR < 0.1 & abs(diff) > 0.5
  )

results_protein |>
  dplyr::count(facade, name = "n_results")
```

    ## # A tibble: 8 × 2
    ##   facade       n_results
    ##   <chr>            <int>
    ## 1 deqms              157
    ## 2 firth              157
    ## 3 limma              155
    ## 4 limma_impute       160
    ## 5 limpa              160
    ## 6 lm                 157
    ## 7 lm_impute          160
    ## 8 lm_missing         160

For facades that combine several underlying result types, such as
`lm_missing`, the `modelName` column still tells you where individual
rows came from.

``` r
results_protein |>
  dplyr::count(facade, contrast, modelName, name = "n_results")
```

    ## # A tibble: 18 × 4
    ##    facade       contrast  modelName          n_results
    ##    <chr>        <chr>     <chr>                  <int>
    ##  1 deqms        A_vs_Ctrl WaldTest_DEqMS            78
    ##  2 deqms        B_vs_Ctrl WaldTest_DEqMS            79
    ##  3 firth        A_vs_Ctrl WaldTestFirth             78
    ##  4 firth        B_vs_Ctrl WaldTestFirth             79
    ##  5 limma        A_vs_Ctrl limma                     78
    ##  6 limma        B_vs_Ctrl limma                     77
    ##  7 limma_impute A_vs_Ctrl limma                     80
    ##  8 limma_impute B_vs_Ctrl limma                     80
    ##  9 limpa        A_vs_Ctrl limpa                     80
    ## 10 limpa        B_vs_Ctrl limpa                     80
    ## 11 lm           A_vs_Ctrl WaldTest_moderated        78
    ## 12 lm           B_vs_Ctrl WaldTest_moderated        79
    ## 13 lm_impute    A_vs_Ctrl WaldTest_moderated        80
    ## 14 lm_impute    B_vs_Ctrl WaldTest_moderated        80
    ## 15 lm_missing   A_vs_Ctrl WaldTest_moderated        78
    ## 16 lm_missing   A_vs_Ctrl groupAverage               2
    ## 17 lm_missing   B_vs_Ctrl WaldTest_moderated        79
    ## 18 lm_missing   B_vs_Ctrl groupAverage               1

## Protein-level volcano comparison

### Standard facades

``` r
standard_facades <- c("lm", "limma", "deqms")
results_standard <- results_protein |>
  dplyr::filter(facade %in% standard_facades)

ggplot(results_standard, aes(x = diff, y = -log10(p.value), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  facet_grid(contrast ~ facade, scales = "free_y") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey70")) +
  labs(x = "log2 fold change", y = "-log10(p.value)", color = "FDR < 0.1\nand |diff| > 0.5") +
  theme_minimal(base_size = 12)
```

![Volcano plots for the standard protein-level facades (lm, limma,
rlm-based).](ContrastFacades_files/figure-html/volcano_standard-1.png)

Volcano plots for the standard protein-level facades (lm, limma,
rlm-based).

### Imputation and missingness facades

Rescued proteins (missing in plain `lm`) are shown as large red diamonds
so they stand out clearly. This group includes the LOD imputation
facades (`lm_missing`, `lm_impute`, `limma_impute`) and the firth
logistic regression facade which models missingness directly.

``` r
impute_facades <- c("lm_missing", "lm_impute", "limma_impute", "firth")
results_impute <- results_protein |>
  dplyr::filter(facade %in% impute_facades) |>
  dplyr::mutate(facade = factor(facade, levels = impute_facades))

ggplot(results_impute, aes(x = diff, y = -log10(p.value))) +
  geom_point(data = dplyr::filter(results_impute, !rescued),
             aes(color = significant), alpha = 0.5, size = 1.2) +
  geom_point(data = dplyr::filter(results_impute, rescued),
             color = "red", shape = 18, size = 5, alpha = 1) +
  facet_grid(contrast ~ facade, scales = "free_y") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey70")) +
  labs(x = "log2 fold change", y = "-log10(p.value)", color = "FDR < 0.1\nand |diff| > 0.5") +
  theme_minimal(base_size = 12)
```

![Volcano plots for the imputation and missingness facades. Large red
diamonds mark proteins rescued (missing in plain
lm).](ContrastFacades_files/figure-html/volcano_impute-1.png)

Volcano plots for the imputation and missingness facades. Large red
diamonds mark proteins rescued (missing in plain lm).

## Looking at the strongest protein-level hits

``` r
results_protein |>
  dplyr::group_by(facade, contrast) |>
  dplyr::slice_min(order_by = p.value, n = 5, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(facade, contrast, modelName, protein_Id, diff, p.value, FDR)
```

    ## # A tibble: 80 × 7
    ##    facade contrast  modelName      protein_Id    diff  p.value      FDR
    ##    <chr>  <chr>     <chr>          <chr>        <dbl>    <dbl>    <dbl>
    ##  1 deqms  A_vs_Ctrl WaldTest_DEqMS Zci7Jw~7064 -0.722 3.79e-12 2.96e-10
    ##  2 deqms  A_vs_Ctrl WaldTest_DEqMS 6TevMr~7550  0.765 3.46e-11 1.35e- 9
    ##  3 deqms  A_vs_Ctrl WaldTest_DEqMS 4Y4DYT~0927  0.583 5.59e-11 1.45e- 9
    ##  4 deqms  A_vs_Ctrl WaldTest_DEqMS 0CubNR~0890  0.674 9.73e-11 1.90e- 9
    ##  5 deqms  A_vs_Ctrl WaldTest_DEqMS KVkccD~1805 -0.524 1.25e- 9 1.95e- 8
    ##  6 deqms  B_vs_Ctrl WaldTest_DEqMS fylZqB~3883  0.717 1.05e-13 8.28e-12
    ##  7 deqms  B_vs_Ctrl WaldTest_DEqMS XxJoJB~7286  0.587 4.63e-12 1.83e-10
    ##  8 deqms  B_vs_Ctrl WaldTest_DEqMS f0Cvvj~6658  0.345 5.39e-10 1.42e- 8
    ##  9 deqms  B_vs_Ctrl WaldTest_DEqMS TR3Ksv~1492 -0.413 3.73e- 9 7.37e- 8
    ## 10 deqms  B_vs_Ctrl WaldTest_DEqMS 4Y4DYT~0927  0.424 8.16e- 9 1.29e- 7
    ## # ℹ 70 more rows

## Proteins that could not be estimated

Every facade has a `get_missing()` method that returns the protein ×
contrast pairs present in the input data but absent from
`get_contrasts()`. This makes it easy to see which proteins each method
fails on and to compare coverage.

``` r
missing_limpa <- if (requireNamespace("limpa", quietly = TRUE)) {
  fa_limpa$get_missing() |> dplyr::mutate(facade = "limpa")
} else {
  data.frame()
}

missing_all <- dplyr::bind_rows(
  fa_lm$get_missing() |> dplyr::mutate(facade = "lm"),
  fa_limma$get_missing() |> dplyr::mutate(facade = "limma"),
  fa_limma_impute$get_missing() |> dplyr::mutate(facade = "limma_impute"),
  fa_lm_missing$get_missing() |> dplyr::mutate(facade = "lm_missing"),
  fa_lm_impute$get_missing() |> dplyr::mutate(facade = "lm_impute"),
  fa_deqms$get_missing() |> dplyr::mutate(facade = "deqms"),
  fa_firth_protein$get_missing() |> dplyr::mutate(facade = "firth"),
  missing_limpa
)

missing_all |>
  dplyr::count(facade, contrast, name = "n_missing") |>
  knitr::kable(caption = "Number of missing protein × contrast pairs per facade")
```

| facade | contrast  | n_missing |
|:-------|:----------|----------:|
| deqms  | A_vs_Ctrl |         2 |
| deqms  | B_vs_Ctrl |         1 |
| firth  | A_vs_Ctrl |         2 |
| firth  | B_vs_Ctrl |         1 |
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
fills in these gaps via group-mean imputation, while `lm_impute` and
`limma_impute` re-fit after imputing individual values with the LOD and
borrowing covariance from successful fits. These should have fewer
missing proteins than plain `lm` or `limma`.

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
      facade %in% c("lm_missing", "lm_impute", "limma_impute", "limpa")
    ) |>
    dplyr::arrange(protein_Id, contrast, facade)

  rescued |>
    knitr::kable(
      digits = 3,
      caption = "Contrast estimates from lm_missing, lm_impute, and limma_impute for proteins that plain lm could not estimate"
    )
} else {
  cat("All proteins were estimable by plain lm — no rescued estimates to show.")
}
```

| facade       | modelName          | protein_Id  | contrast  | avgAbd |   diff |   FDR | statistic | std.error |     df | p.value | conf.low | conf.high | sigma | rescued | significant |
|:-------------|:-------------------|:------------|:----------|-------:|-------:|------:|----------:|----------:|-------:|--------:|---------:|----------:|------:|:--------|:------------|
| limma_impute | limma              | 8mS8sK~0150 | A_vs_Ctrl |  3.776 |  0.000 | 1.000 |     0.000 |     0.063 |  4.468 |   1.000 |   -0.167 |     0.167 | 0.089 | TRUE    | FALSE       |
| limpa        | limpa              | 8mS8sK~0150 | A_vs_Ctrl |  2.798 | -0.615 | 0.226 |    -1.435 |     0.429 | 30.965 |   0.161 |   -1.489 |     0.259 | 0.965 | TRUE    | FALSE       |
| lm_impute    | WaldTest_moderated | 8mS8sK~0150 | A_vs_Ctrl |  3.776 |  0.000 | 1.000 |     0.000 |     0.065 |  4.310 |   1.000 |   -0.234 |     0.234 | 0.087 | TRUE    | FALSE       |
| lm_missing   | groupAverage       | 8mS8sK~0150 | A_vs_Ctrl |     NA |     NA |    NA |        NA |     0.102 |  2.000 |      NA |       NA |        NA | 0.102 | TRUE    | NA          |
| limma_impute | limma              | 8mS8sK~0150 | B_vs_Ctrl |  3.784 |  0.018 | 0.803 |     0.279 |     0.063 |  4.468 |   0.792 |   -0.150 |     0.185 | 0.089 | FALSE   | FALSE       |
| limpa        | limpa              | 8mS8sK~0150 | B_vs_Ctrl |  3.245 |  0.279 | 0.578 |     0.662 |     0.422 | 30.965 |   0.513 |   -0.581 |     1.140 | 0.965 | FALSE   | FALSE       |
| lm_impute    | WaldTest_moderated | 8mS8sK~0150 | B_vs_Ctrl |  3.784 |  0.018 | 0.798 |     0.286 |     0.065 |  4.310 |   0.788 |   -0.217 |     0.252 | 0.087 | FALSE   | FALSE       |
| lm_missing   | WaldTest_moderated | 8mS8sK~0150 | B_vs_Ctrl |  3.632 |  0.339 | 0.020 |     3.690 |     0.102 |  5.377 |   0.012 |    0.108 |     0.570 | 0.092 | FALSE   | FALSE       |
| limma_impute | limma              | DTCi0N~0734 | A_vs_Ctrl |  3.902 | -0.253 | 0.017 |    -3.996 |     0.063 |  6.468 |   0.006 |   -0.405 |    -0.101 | 0.090 | TRUE    | FALSE       |
| limpa        | limpa              | DTCi0N~0734 | A_vs_Ctrl |  3.550 | -0.982 | 0.020 |    -2.714 |     0.362 | 30.965 |   0.011 |   -1.719 |    -0.244 | 0.991 | TRUE    | TRUE        |
| lm_impute    | WaldTest_moderated | DTCi0N~0734 | A_vs_Ctrl |  3.902 | -0.253 | 0.017 |    -4.057 |     0.065 |  6.310 |   0.006 |   -0.467 |    -0.040 | 0.088 | TRUE    | FALSE       |
| lm_missing   | groupAverage       | DTCi0N~0734 | A_vs_Ctrl |     NA |     NA |    NA |        NA |     0.057 |  4.000 |      NA |       NA |        NA | 0.070 | TRUE    | NA          |
| limma_impute | limma              | DTCi0N~0734 | B_vs_Ctrl |  4.112 |  0.166 | 0.051 |     2.626 |     0.063 |  6.468 |   0.037 |    0.014 |     0.319 | 0.090 | FALSE   | FALSE       |
| limpa        | limpa              | DTCi0N~0734 | B_vs_Ctrl |  4.145 |  0.208 | 0.619 |     0.581 |     0.358 | 30.965 |   0.565 |   -0.522 |     0.939 | 0.991 | FALSE   | FALSE       |
| lm_impute    | WaldTest_moderated | DTCi0N~0734 | B_vs_Ctrl |  4.112 |  0.166 | 0.049 |     2.665 |     0.065 |  6.310 |   0.035 |   -0.047 |     0.380 | 0.088 | FALSE   | FALSE       |
| lm_missing   | WaldTest_moderated | DTCi0N~0734 | B_vs_Ctrl |  4.224 |  0.222 | 0.016 |     3.499 |     0.057 |  7.377 |   0.009 |    0.040 |     0.403 | 0.078 | FALSE   | FALSE       |
| limma_impute | limma              | OrL0ux~1369 | A_vs_Ctrl |  3.879 | -0.207 | 0.050 |    -3.293 |     0.063 |  4.468 |   0.026 |   -0.374 |    -0.039 | 0.089 | FALSE   | FALSE       |
| limpa        | limpa              | OrL0ux~1369 | A_vs_Ctrl |  3.497 | -0.881 | 0.025 |    -2.630 |     0.335 | 30.965 |   0.013 |   -1.563 |    -0.198 | 0.960 | FALSE   | TRUE        |
| lm_impute    | WaldTest_moderated | OrL0ux~1369 | A_vs_Ctrl |  3.879 | -0.207 | 0.049 |    -3.369 |     0.065 |  4.310 |   0.025 |   -0.441 |     0.028 | 0.087 | FALSE   | FALSE       |
| lm_missing   | WaldTest_moderated | OrL0ux~1369 | A_vs_Ctrl |  3.913 | -0.276 | 0.055 |    -2.951 |     0.084 |  5.340 |   0.029 |   -0.480 |    -0.072 | 0.081 | FALSE   | FALSE       |
| limma_impute | limma              | OrL0ux~1369 | B_vs_Ctrl |  3.879 | -0.207 | 0.038 |    -3.293 |     0.063 |  4.468 |   0.026 |   -0.374 |    -0.039 | 0.089 | TRUE    | FALSE       |
| limpa        | limpa              | OrL0ux~1369 | B_vs_Ctrl |  3.297 | -1.281 | 0.006 |    -3.200 |     0.400 | 30.965 |   0.003 |   -2.097 |    -0.464 | 0.960 | TRUE    | TRUE        |
| lm_impute    | WaldTest_moderated | OrL0ux~1369 | B_vs_Ctrl |  3.879 | -0.207 | 0.036 |    -3.369 |     0.065 |  4.310 |   0.025 |   -0.441 |     0.028 | 0.087 | TRUE    | FALSE       |
| lm_missing   | groupAverage       | OrL0ux~1369 | B_vs_Ctrl |     NA |     NA |    NA |        NA |     0.060 |  2.000 |      NA |       NA |        NA | 0.073 | TRUE    | NA          |

Contrast estimates from lm_missing, lm_impute, and limma_impute for
proteins that plain lm could not estimate

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
- the `limpa` facade uses its own DPC-based aggregation
  (`AggregateLimpa`) which produces standard errors propagated into
  vooma precision weights
- the peptide-level facades now explicitly require lower-level hierarchy
  below the analysis subject
- the shared facade API still makes it straightforward to compare
  methods once the data level is chosen consistently

The results are comparable at the API level, but the comparison is only
meaningful when methods using the same biological unit are plotted
together.
