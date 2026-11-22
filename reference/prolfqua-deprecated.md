# Deprecated functions and arguments

These functions and arguments were renamed or retired in prolfqua 1.9.5.
The old names still work and warn; following the policy in
[prolfqua-package](https://fgcz.github.io/prolfqua/reference/prolfqua-package.md),
they call [`.Defunct()`](https://rdrr.io/r/base/Defunct.html) in a later
release.

## Usage

``` r
LFQDataToSummarizedExperiment(...)

LR_test(...)

lfq_power_t_test_quantiles_V2(...)

R6_extract_values(...)

list_to_AnalysisConfiguration(...)

list_to_ContrastConfiguration(...)

summary_ROPECA_median_p_scaled(...)

missigness_histogram(...)
```

## Arguments

- ...:

  passed to the replacement function

## Value

the value of the replacement function

## Details

Functions:

- `LFQDataToSummarizedExperiment()`: use
  [`lfqdata_to_summarized_experiment()`](https://fgcz.github.io/prolfqua/reference/lfqdata_to_summarized_experiment.md)

- `LR_test()`: use
  [`lr_test()`](https://fgcz.github.io/prolfqua/reference/lr_test.md)

- `lfq_power_t_test_quantiles_V2()`: use
  [`lfq_power_t_test_quantiles()`](https://fgcz.github.io/prolfqua/reference/lfq_power_t_test_quantiles.md)

- `R6_extract_values()`: use
  [`r6_extract_values()`](https://fgcz.github.io/prolfqua/reference/r6_extract_values.md)

- `list_to_AnalysisConfiguration()`: use
  [`list_to_analysis_configuration()`](https://fgcz.github.io/prolfqua/reference/list_to_analysis_configuration.md)

- `list_to_ContrastConfiguration()`: use
  [`list_to_contrast_configuration()`](https://fgcz.github.io/prolfqua/reference/list_to_contrast_configuration.md)

- `summary_ROPECA_median_p_scaled()`: use
  [`summary_ropeca_median_p_scaled()`](https://fgcz.github.io/prolfqua/reference/summary_ropeca_median_p_scaled.md)

- `missigness_histogram()` and `LFQDataPlotter$missigness_histogram()`:
  use
  [`missingness_histogram()`](https://fgcz.github.io/prolfqua/reference/missingness_histogram.md)
  and `LFQDataPlotter$missingness_histogram()`

Simulator arguments of
[`sim_lfq_data()`](https://fgcz.github.io/prolfqua/reference/sim_lfq_data.md),
[`sim_lfq_data_peptide_config()`](https://fgcz.github.io/prolfqua/reference/sim_lfq_data_peptide_config.md),
[`sim_lfq_data_protein_config()`](https://fgcz.github.io/prolfqua/reference/sim_lfq_data_protein_config.md),
[`sim_lfq_data_2factor_config()`](https://fgcz.github.io/prolfqua/reference/sim_lfq_data_2factor_config.md),
[`sim_build_models_lm()`](https://fgcz.github.io/prolfqua/reference/sim_build_models_lm.md),
[`sim_build_models_lmer()`](https://fgcz.github.io/prolfqua/reference/sim_build_models_lmer.md)
and
[`sim_build_models_logistf()`](https://fgcz.github.io/prolfqua/reference/sim_build_models_logistf.md):

- `Nprot`: use `n_proteins`

- `N`: use `n_per_group`

- `PEPTIDE`: use `peptide`

- `TWO`: use `two_factor`

Arguments of
[`volcano_plotly()`](https://fgcz.github.io/prolfqua/reference/volcano_plotly.md)
and
[`multigroup_volcano()`](https://fgcz.github.io/prolfqua/reference/multigroup_volcano.md):

- `.data`: use `data`

Argument of
[`merge_contrasts_results()`](https://fgcz.github.io/prolfqua/reference/merge_contrasts_results.md):

- `model_name`: never used; the merged model is named after its inputs,
  `<prefer>_<add>`
