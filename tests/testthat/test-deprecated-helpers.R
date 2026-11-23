## Helpers unexported in prolfqua 1.9.5 (rOpenSci self-review WP19.3): exported wrappers warn and forward unchanged

deprecated_helpers <- c(
  "complete_cases",
  "hierarchy_counts",
  "summarize_hierarchy",
  "sample_subset",
  "tidy_to_wide_config",
  "table_factors",
  "plot_heatmap",
  "plot_heatmap_cor",
  "plot_hierarchies_boxplot",
  "plot_hierarchies_boxplot_df",
  "plot_intensity_distribution_violin",
  "plot_intensity_distribution_density",
  "plot_na_heatmap",
  "plot_pca",
  "plot_raster",
  "plot_sample_correlation",
  "pairs_smooth",
  "panel_cor",
  "missingness_histogram",
  "upset_missing_stats",
  "summarize_stats",
  "summarize_stats_quantiles",
  "lfq_power_t_test_proteins",
  "lfq_power_t_test_quantiles",
  "plot_stat_density",
  "plot_stat_density_median",
  "plot_stat_violin",
  "plot_stat_violin_median",
  "plot_stdv_vs_mean",
  "poolvar",
  "compute_pooled",
  "table_factors_size",
  "missingness_per_condition",
  "missingness_per_condition_cumsum",
  "robust_scale",
  "get_robscales",
  "scale_with_subset",
  "apply_to_response_matrix",
  "transform_work_intensity",
  "response_matrix_as_tibble",
  "estimate_intensity",
  "aggregate_intensity_top_n",
  "rank_peptide_by_intensity",
  "plot_estimate",
  "plot_hierarchies_line",
  "make_reduced_hierarchy_config",
  "compute_contrast",
  "compute_lmer_contrast",
  "contrasts_linfct",
  "linfct_from_model",
  "linfct_matrix_contrasts",
  "get_complete_model_fit",
  "get_anova_df",
  "model_analyse",
  "model_summary",
  "build_model_logistf",
  "encode_bin_resp",
  "make_interaction_column",
  "moderated_p_limma",
  "moderated_p_deqms",
  "moderated_p_deqms_long",
  "summary_ropeca_median_p_scaled",
  "get_p_values_pbeta",
  "multigroup_volcano",
  "sim_build_models_lm",
  "sim_build_models_lmer",
  "sim_build_models_logistf",
  "sim_make_model_lm",
  "sim_make_model_lmer",
  "script_copy_helper_vec",
  "separate_hierarchy",
  "get_contrast",
  "matrix_to_tibble"
)

test_that("every unexported helper is still exported for one release and deprecates first", {
  exports <- getNamespaceExports("prolfqua")
  expect_length(deprecated_helpers, 73)
  expect_setequal(intersect(deprecated_helpers, exports), deprecated_helpers)
  expect_true("ContrastsFacadeBase" %in% exports)
  # covr rewrites function bodies, so the first call is its counter, not .Deprecated()
  skip_on_covr()
  first_call <- vapply(
    deprecated_helpers,
    function(f) {
      deparse(body(get(f, envir = asNamespace("prolfqua")))[[2]][[1]])
    },
    character(1)
  )
  expect_true(all(first_call == ".Deprecated"), info = toString(names(first_call)[first_call != ".Deprecated"]))
})

test_that("deprecated helpers forward the call unchanged and name the public route", {
  sim <- sim_lfq_data_peptide_config(n_proteins = 5)
  lfq <- LFQData$new(sim$data, sim$config)
  lfq_log <- lfq$get_Transformer()$log2()$lfq

  expect_warning(old <- summarize_stats(lfq), "LFQData\\$get_Stats\\(\\)\\$stats\\(\\)", class = "deprecatedWarning")
  expect_identical(old, .summarize_stats(lfq))

  # a default that refers to another argument (hierarchy = hierarchy_keys) is evaluated by the internal function
  expect_warning(
    old <- summarize_hierarchy(lfq$data_long(), lfq$hierarchy_keys()),
    "LFQData\\$summarize_hierarchy\\(\\)",
    class = "deprecatedWarning"
  )
  expect_identical(old, .summarize_hierarchy(lfq$data_long(), lfq$hierarchy_keys()))

  expect_warning(old <- tidy_to_wide_config(lfq, as.matrix = TRUE), class = "deprecatedWarning")
  expect_identical(old, .tidy_to_wide_config(lfq, as.matrix = TRUE))

  expect_warning(old <- suppressMessages(estimate_intensity(lfq_log)), class = "deprecatedWarning")
  expect_identical(old$data, suppressMessages(.estimate_intensity(lfq_log))$data)
})

test_that("helpers that leave prolfqua say where they went", {
  expect_warning(
    suppressMessages(script_copy_helper_vec("CITATION", workdir = withr::local_tempdir())),
    "prolfquapp::script_copy_helper_vec",
    class = "deprecatedWarning"
  )
  data <- tibble::tibble(protein_Id = "p1", group_A = 2, group_B = 1)
  expect_warning(
    suppressMessages(get_contrast(data, "protein_Id", c(AvsB = "group_A - group_B"))),
    "only prolfquatakeoutmissing used it",
    class = "deprecatedWarning"
  )
  expect_warning(
    matrix_to_tibble(matrix(1, 1, 1, dimnames = list("a", "x"))),
    "tibble::as_tibble",
    class = "deprecatedWarning"
  )
})
