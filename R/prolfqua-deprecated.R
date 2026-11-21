#' Deprecated functions and arguments
#'
#' These functions and arguments were renamed in prolfqua 1.9.5. The old names still work and warn; following the
#' policy in [prolfqua-package], they call `.Defunct()` in a later release.
#'
#' Functions:
#' * `LFQDataToSummarizedExperiment()`: use [lfqdata_to_summarized_experiment()]
#' * `LR_test()`: use [lr_test()]
#' * `lfq_power_t_test_quantiles_V2()`: use [lfq_power_t_test_quantiles()]
#' * `R6_extract_values()`: use [r6_extract_values()]
#' * `list_to_AnalysisConfiguration()`: use [list_to_analysis_configuration()]
#' * `list_to_ContrastConfiguration()`: use [list_to_contrast_configuration()]
#' * `summary_ROPECA_median_p_scaled()`: use [summary_ropeca_median_p_scaled()]
#' * `missigness_histogram()` and `LFQDataPlotter$missigness_histogram()`: use [missingness_histogram()] and
#'   `LFQDataPlotter$missingness_histogram()`
#'
#' Simulator arguments of [sim_lfq_data()], [sim_lfq_data_peptide_config()], [sim_lfq_data_protein_config()],
#' [sim_lfq_data_2factor_config()], [sim_build_models_lm()], [sim_build_models_lmer()] and
#' [sim_build_models_logistf()]:
#' * `Nprot`: use `n_proteins`
#' * `N`: use `n_per_group`
#' * `PEPTIDE`: use `peptide`
#' * `TWO`: use `two_factor`
#'
#' Arguments of [volcano_plotly()] and [multigroup_volcano()]:
#' * `.data`: use `data`
#'
#' @param ... passed to the replacement function
#' @return the value of the replacement function
#' @name prolfqua-deprecated
#' @keywords internal
NULL

#' @rdname prolfqua-deprecated
#' @export
LFQDataToSummarizedExperiment <- function(...) {
  .Deprecated("lfqdata_to_summarized_experiment")
  lfqdata_to_summarized_experiment(...)
}

#' @rdname prolfqua-deprecated
#' @export
LR_test <- function(...) {
  .Deprecated("lr_test")
  lr_test(...)
}

#' @rdname prolfqua-deprecated
#' @export
lfq_power_t_test_quantiles_V2 <- function(...) {
  .Deprecated("lfq_power_t_test_quantiles")
  lfq_power_t_test_quantiles(...)
}

#' @rdname prolfqua-deprecated
#' @export
R6_extract_values <- function(...) {
  .Deprecated("r6_extract_values")
  r6_extract_values(...)
}

#' @rdname prolfqua-deprecated
#' @export
list_to_AnalysisConfiguration <- function(...) {
  .Deprecated("list_to_analysis_configuration")
  list_to_analysis_configuration(...)
}

#' @rdname prolfqua-deprecated
#' @export
list_to_ContrastConfiguration <- function(...) {
  .Deprecated("list_to_contrast_configuration")
  list_to_contrast_configuration(...)
}

#' @rdname prolfqua-deprecated
#' @export
summary_ROPECA_median_p_scaled <- function(...) {
  .Deprecated("summary_ropeca_median_p_scaled")
  summary_ropeca_median_p_scaled(...)
}

#' @rdname prolfqua-deprecated
#' @export
missigness_histogram <- function(...) {
  .Deprecated("missingness_histogram")
  missingness_histogram(...)
}

#' Warn that an argument is deprecated and return its value for the replacement
#' @param old deprecated argument name
#' @param new replacement argument name
#' @param value value given to the deprecated argument
#' @noRd
.deprecate_arg <- function(old, new, value) {
  .Deprecated(msg = sprintf("Argument `%s` is deprecated; use `%s` instead.", old, new))
  value
}
