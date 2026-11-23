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
  `LFQData$get_Stats()$power_t_test_quantiles()`

- `R6_extract_values()`: use
  [`r6_extract_values()`](https://fgcz.github.io/prolfqua/reference/r6_extract_values.md)

- `list_to_AnalysisConfiguration()`: use
  [`list_to_analysis_configuration()`](https://fgcz.github.io/prolfqua/reference/list_to_analysis_configuration.md)

- `list_to_ContrastConfiguration()`: use
  [`list_to_contrast_configuration()`](https://fgcz.github.io/prolfqua/reference/list_to_contrast_configuration.md)

- `summary_ROPECA_median_p_scaled()`: use
  `build_contrast_analysis(method = "ropeca_nested")`

- `missigness_histogram()` and `LFQDataPlotter$missigness_histogram()`:
  use `LFQData$get_Plotter()$missingness_histogram()`

Helpers unexported in prolfqua 1.9.5. They were documented only on
internal help pages; each has a public route, and the deprecation
warning names it. Their help pages are kept for this release.

- `LFQData` methods:
  [`complete_cases()`](https://fgcz.github.io/prolfqua/reference/complete_cases.md),
  [`hierarchy_counts()`](https://fgcz.github.io/prolfqua/reference/hierarchy_counts.md),
  [`summarize_hierarchy()`](https://fgcz.github.io/prolfqua/reference/summarize_hierarchy.md),
  [`sample_subset()`](https://fgcz.github.io/prolfqua/reference/sample_subset.md)
  (`$get_sample()`),
  [`tidy_to_wide_config()`](https://fgcz.github.io/prolfqua/reference/tidy_to_wide_config.md)
  (`$data_wide()`),
  [`table_factors()`](https://fgcz.github.io/prolfqua/reference/table_factors.md)
  and
  [`table_factors_size()`](https://fgcz.github.io/prolfqua/reference/table_factors_size.md)
  (`$factors()`)

- `LFQData$get_Plotter()`:
  [`plot_heatmap()`](https://fgcz.github.io/prolfqua/reference/plot_heatmap.md),
  [`plot_heatmap_cor()`](https://fgcz.github.io/prolfqua/reference/plot_heatmap_cor.md),
  [`plot_hierarchies_boxplot()`](https://fgcz.github.io/prolfqua/reference/plot_hierarchies_boxplot.md),
  [`plot_hierarchies_boxplot_df()`](https://fgcz.github.io/prolfqua/reference/plot_hierarchies_boxplot_df.md),
  [`plot_intensity_distribution_violin()`](https://fgcz.github.io/prolfqua/reference/plot_intensity_distribution_violin.md),
  [`plot_intensity_distribution_density()`](https://fgcz.github.io/prolfqua/reference/plot_intensity_distribution_violin.md),
  [`plot_na_heatmap()`](https://fgcz.github.io/prolfqua/reference/plot_na_heatmap.md),
  [`plot_pca()`](https://fgcz.github.io/prolfqua/reference/plot_pca.md),
  [`plot_raster()`](https://fgcz.github.io/prolfqua/reference/plot_raster.md),
  [`plot_sample_correlation()`](https://fgcz.github.io/prolfqua/reference/plot_sample_correlation.md),
  [`pairs_smooth()`](https://fgcz.github.io/prolfqua/reference/pairs_smooth.md),
  [`panel_cor()`](https://fgcz.github.io/prolfqua/reference/panel_cor.md),
  [`missingness_histogram()`](https://fgcz.github.io/prolfqua/reference/missingness_histogram.md),
  [`upset_missing_stats()`](https://fgcz.github.io/prolfqua/reference/upset_missing_stats.md)

- `LFQData$get_Stats()`:
  [`summarize_stats()`](https://fgcz.github.io/prolfqua/reference/summarize_stats.md),
  [`summarize_stats_quantiles()`](https://fgcz.github.io/prolfqua/reference/summarize_stats.md),
  [`lfq_power_t_test_proteins()`](https://fgcz.github.io/prolfqua/reference/lfq_power_t_test_proteins.md),
  [`lfq_power_t_test_quantiles()`](https://fgcz.github.io/prolfqua/reference/lfq_power_t_test_quantiles.md),
  [`plot_stat_density()`](https://fgcz.github.io/prolfqua/reference/plot_stat_density.md),
  [`plot_stat_density_median()`](https://fgcz.github.io/prolfqua/reference/plot_stat_density_median.md),
  [`plot_stat_violin()`](https://fgcz.github.io/prolfqua/reference/plot_stat_violin.md),
  [`plot_stat_violin_median()`](https://fgcz.github.io/prolfqua/reference/plot_stat_violin_median.md),
  [`plot_stdv_vs_mean()`](https://fgcz.github.io/prolfqua/reference/plot_stdv_vs_mean.md),
  [`poolvar()`](https://fgcz.github.io/prolfqua/reference/pooled_var.md),
  [`compute_pooled()`](https://fgcz.github.io/prolfqua/reference/pooled_var.md)

- `LFQData$get_Summariser()`:
  [`missingness_per_condition()`](https://fgcz.github.io/prolfqua/reference/missingness_per_condition.md),
  [`missingness_per_condition_cumsum()`](https://fgcz.github.io/prolfqua/reference/missingness_per_condition_cumsum.md)

- `LFQData$get_Transformer()`:
  [`robust_scale()`](https://fgcz.github.io/prolfqua/reference/robust_scale.md),
  [`get_robscales()`](https://fgcz.github.io/prolfqua/reference/get_robscales.md),
  [`scale_with_subset()`](https://fgcz.github.io/prolfqua/reference/scale_with_subset.md),
  [`apply_to_response_matrix()`](https://fgcz.github.io/prolfqua/reference/apply_to_response_matrix.md),
  [`transform_work_intensity()`](https://fgcz.github.io/prolfqua/reference/transform_work_intensity.md),
  [`response_matrix_as_tibble()`](https://fgcz.github.io/prolfqua/reference/response_matrix_as_tibble.md)

- `LFQData$get_Aggregator()`:
  [`estimate_intensity()`](https://fgcz.github.io/prolfqua/reference/estimate_intensity.md),
  [`aggregate_intensity_top_n()`](https://fgcz.github.io/prolfqua/reference/aggregate_intensity_top_n.md),
  [`rank_peptide_by_intensity()`](https://fgcz.github.io/prolfqua/reference/rank_peptide_by_intensity.md),
  [`plot_estimate()`](https://fgcz.github.io/prolfqua/reference/plot_estimate.md),
  [`plot_hierarchies_line()`](https://fgcz.github.io/prolfqua/reference/plot_hierarchies_line.md),
  [`make_reduced_hierarchy_config()`](https://fgcz.github.io/prolfqua/reference/make_reduced_hierarchy_config.md)

- [`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md)
  and the [Model](https://fgcz.github.io/prolfqua/reference/Model.md)
  and
  [Contrasts](https://fgcz.github.io/prolfqua/reference/Contrasts.md)
  classes:
  [`compute_contrast()`](https://fgcz.github.io/prolfqua/reference/compute_contrast.md),
  [`compute_lmer_contrast()`](https://fgcz.github.io/prolfqua/reference/compute_lmer_contrast.md),
  [`contrasts_linfct()`](https://fgcz.github.io/prolfqua/reference/contrasts_linfct.md),
  [`linfct_from_model()`](https://fgcz.github.io/prolfqua/reference/linfct_from_model.md),
  [`linfct_matrix_contrasts()`](https://fgcz.github.io/prolfqua/reference/linfct_matrix_contrasts.md),
  [`get_complete_model_fit()`](https://fgcz.github.io/prolfqua/reference/get_complete_model_fit.md),
  [`get_anova_df()`](https://fgcz.github.io/prolfqua/reference/get_anova_df.md),
  [`model_analyse()`](https://fgcz.github.io/prolfqua/reference/model_analyse.md),
  [`model_summary()`](https://fgcz.github.io/prolfqua/reference/model_summary.md),
  [`build_model_logistf()`](https://fgcz.github.io/prolfqua/reference/build_model_logistf.md),
  [`encode_bin_resp()`](https://fgcz.github.io/prolfqua/reference/encode_bin_resp.md),
  [`make_interaction_column()`](https://fgcz.github.io/prolfqua/reference/make_interaction_column.md),
  [`moderated_p_limma()`](https://fgcz.github.io/prolfqua/reference/moderated_p_limma.md),
  [`moderated_p_deqms()`](https://fgcz.github.io/prolfqua/reference/moderated_p_deqms.md),
  [`moderated_p_deqms_long()`](https://fgcz.github.io/prolfqua/reference/moderated_p_deqms_long.md),
  [`summary_ropeca_median_p_scaled()`](https://fgcz.github.io/prolfqua/reference/summary_ropeca_median_p_scaled.md),
  [`get_p_values_pbeta()`](https://fgcz.github.io/prolfqua/reference/get_p_values_pbeta.md),
  [`multigroup_volcano()`](https://fgcz.github.io/prolfqua/reference/multigroup_volcano.md)

- simulation helpers for tests, use
  [`sim_lfq_data_protein_config()`](https://fgcz.github.io/prolfqua/reference/sim_lfq_data_protein_config.md)
  or
  [`sim_lfq_data_peptide_config()`](https://fgcz.github.io/prolfqua/reference/sim_lfq_data_peptide_config.md)
  with
  [`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md):
  [`sim_build_models_lm()`](https://fgcz.github.io/prolfqua/reference/sim_build_models_lm.md),
  [`sim_build_models_lmer()`](https://fgcz.github.io/prolfqua/reference/sim_build_models_lmer.md),
  [`sim_build_models_logistf()`](https://fgcz.github.io/prolfqua/reference/sim_build_models_logistf.md),
  [`sim_make_model_lm()`](https://fgcz.github.io/prolfqua/reference/sim_make_model_lm.md),
  [`sim_make_model_lmer()`](https://fgcz.github.io/prolfqua/reference/sim_make_model_lmer.md)

- moved to prolfquapp:
  [`script_copy_helper_vec()`](https://fgcz.github.io/prolfqua/reference/script_copy_helper_vec.md),
  [`separate_hierarchy()`](https://fgcz.github.io/prolfqua/reference/separate_hierarchy.md)

- [`get_contrast()`](https://fgcz.github.io/prolfqua/reference/get_contrast.md):
  removed; only prolfquatakeoutmissing used it, which has its own copy

- [`matrix_to_tibble()`](https://fgcz.github.io/prolfqua/reference/matrix_to_tibble.md):
  use `tibble::as_tibble(x, rownames = )`

Simulator arguments of
[`sim_lfq_data()`](https://fgcz.github.io/prolfqua/reference/sim_lfq_data.md),
[`sim_lfq_data_peptide_config()`](https://fgcz.github.io/prolfqua/reference/sim_lfq_data_peptide_config.md),
[`sim_lfq_data_protein_config()`](https://fgcz.github.io/prolfqua/reference/sim_lfq_data_protein_config.md)
and
[`sim_lfq_data_2factor_config()`](https://fgcz.github.io/prolfqua/reference/sim_lfq_data_2factor_config.md):

- `Nprot`: use `n_proteins`

- `N`: use `n_per_group`

- `PEPTIDE`: use `peptide`

- `TWO`: use `two_factor`

Argument of
[`volcano_plotly()`](https://fgcz.github.io/prolfqua/reference/volcano_plotly.md):

- `.data`: use `data`

Argument of
[`merge_contrasts_results()`](https://fgcz.github.io/prolfqua/reference/merge_contrasts_results.md):

- `model_name`: never used; the merged model is named after its inputs,
  `<prefer>_<add>`
