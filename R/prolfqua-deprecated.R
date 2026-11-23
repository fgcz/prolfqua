#' Deprecated functions and arguments
#'
#' These functions and arguments were renamed or retired in prolfqua 1.9.5. The old names still work and warn;
#' following the policy in [prolfqua-package], they call `.Defunct()` in a later release.
#'
#' Functions:
#' * `LFQDataToSummarizedExperiment()`: use [lfqdata_to_summarized_experiment()]
#' * `LR_test()`: use [lr_test()]
#' * `lfq_power_t_test_quantiles_V2()`: use `LFQData$get_Stats()$power_t_test_quantiles()`
#' * `R6_extract_values()`: use [r6_extract_values()]
#' * `list_to_AnalysisConfiguration()`: use [list_to_analysis_configuration()]
#' * `list_to_ContrastConfiguration()`: use [list_to_contrast_configuration()]
#' * `summary_ROPECA_median_p_scaled()`: use `build_contrast_analysis(method = "ropeca_nested")`
#' * `missigness_histogram()` and `LFQDataPlotter$missigness_histogram()`: use
#'   `LFQData$get_Plotter()$missingness_histogram()`
#'
#' Helpers unexported in prolfqua 1.9.5. They were documented only on internal help pages; each has a public
#' route, and the deprecation warning names it. Their help pages are kept for this release.
#' * `LFQData` methods: `complete_cases()`, `hierarchy_counts()`, `summarize_hierarchy()`, `sample_subset()`
#'   (`$get_sample()`), `tidy_to_wide_config()` (`$data_wide()`), `table_factors()` and `table_factors_size()`
#'   (`$factors()`)
#' * `LFQData$get_Plotter()`: `plot_heatmap()`, `plot_heatmap_cor()`, `plot_hierarchies_boxplot()`,
#'   `plot_hierarchies_boxplot_df()`, `plot_intensity_distribution_violin()`, `plot_intensity_distribution_density()`,
#'   `plot_na_heatmap()`, `plot_pca()`, `plot_raster()`, `plot_sample_correlation()`, `pairs_smooth()`,
#'   `panel_cor()`, `missingness_histogram()`, `upset_missing_stats()`
#' * `LFQData$get_Stats()`: `summarize_stats()`, `summarize_stats_quantiles()`, `lfq_power_t_test_proteins()`,
#'   `lfq_power_t_test_quantiles()`, `plot_stat_density()`, `plot_stat_density_median()`, `plot_stat_violin()`,
#'   `plot_stat_violin_median()`, `plot_stdv_vs_mean()`, `poolvar()`, `compute_pooled()`
#' * `LFQData$get_Summariser()`: `missingness_per_condition()`, `missingness_per_condition_cumsum()`
#' * `LFQData$get_Transformer()`: `robust_scale()`, `get_robscales()`, `scale_with_subset()`,
#'   `apply_to_response_matrix()`, `transform_work_intensity()`, `response_matrix_as_tibble()`
#' * `LFQData$get_Aggregator()`: `estimate_intensity()`, `aggregate_intensity_top_n()`,
#'   `rank_peptide_by_intensity()`, `plot_estimate()`, `plot_hierarchies_line()`, `make_reduced_hierarchy_config()`
#' * [build_contrast_analysis()] and the [Model] and [Contrasts] classes: `compute_contrast()`,
#'   `compute_lmer_contrast()`, `contrasts_linfct()`, `linfct_from_model()`, `linfct_matrix_contrasts()`,
#'   `get_complete_model_fit()`, `get_anova_df()`, `model_analyse()`, `model_summary()`, `build_model_logistf()`,
#'   `encode_bin_resp()`, `make_interaction_column()`, `moderated_p_limma()`, `moderated_p_deqms()`,
#'   `moderated_p_deqms_long()`, `summary_ropeca_median_p_scaled()`, `get_p_values_pbeta()`, `multigroup_volcano()`
#' * simulation helpers for tests, use [sim_lfq_data_protein_config()] or [sim_lfq_data_peptide_config()] with
#'   [build_contrast_analysis()]: `sim_build_models_lm()`, `sim_build_models_lmer()`, `sim_build_models_logistf()`,
#'   `sim_make_model_lm()`, `sim_make_model_lmer()`
#' * moved to prolfquapp: `script_copy_helper_vec()`, `separate_hierarchy()`
#' * `get_contrast()`: removed; only prolfquatakeoutmissing used it, which has its own copy
#' * `matrix_to_tibble()`: use `tibble::as_tibble(x, rownames = )`
#'
#' Simulator arguments of [sim_lfq_data()], [sim_lfq_data_peptide_config()], [sim_lfq_data_protein_config()] and
#' [sim_lfq_data_2factor_config()]:
#' * `Nprot`: use `n_proteins`
#' * `N`: use `n_per_group`
#' * `PEPTIDE`: use `peptide`
#' * `TWO`: use `two_factor`
#'
#' Argument of [volcano_plotly()]:
#' * `.data`: use `data`
#'
#' Argument of [merge_contrasts_results()]:
#' * `model_name`: never used; the merged model is named after its inputs, `<prefer>_<add>`
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
  .Deprecated("LFQData$get_Stats()$power_t_test_quantiles()", package = "prolfqua")
  .lfq_power_t_test_quantiles(...)
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
  .Deprecated("build_contrast_analysis(method = \"ropeca_nested\")", package = "prolfqua")
  .summary_ropeca_median_p_scaled(...)
}

#' @rdname prolfqua-deprecated
#' @export
missigness_histogram <- function(...) {
  .Deprecated("LFQData$get_Plotter()$missingness_histogram()", package = "prolfqua")
  .missingness_histogram(...)
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

# Helpers unexported in prolfqua 1.9.5 (rOpenSci self-review WP19.3) ----

#' Complete cases
#'
#' The tidy table does not need to contain missing data.
#' This function re-establishes the missing observations in a sample.
#'
#' @param lfqdata [LFQData] object with the data and its configuration
#' @return a tibble of the long-format data with one row for every combination of sample and hierarchy element;
#'   added rows have `NA` intensities
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$complete_cases()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
complete_cases <- function(lfqdata) {
  .Deprecated("LFQData$complete_cases()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .complete_cases
  eval(call, parent.frame())
}

#' Count distinct elements for each level of hierarchy
#'
#' E.g. number of proteins, peptides, precursors in the dataset. Counts over the
#' whole dataset, whether observed or not.
#'
#' @param pdata data.frame
#' @param hierarchy_keys character vector - all hierarchy column names
#' @return a one-row tibble with one column per hierarchy key holding the number of distinct elements
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$hierarchy_counts()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
hierarchy_counts <- function(pdata, hierarchy_keys) {
  .Deprecated("LFQData$hierarchy_counts()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .hierarchy_counts
  eval(call, parent.frame())
}

#' Summarize hierarchy counts
#'
#' E.g compute number of peptides for each protein: for each entry of
#' `hierarchy` (and of `factors`), the number of distinct entries of the other
#' hierarchy levels, in columns named `<level>_n`. Counts experiment-wide,
#' whether observed or not.
#'
#' @param pdata data.frame
#' @param hierarchy_keys character vector - all hierarchy column names
#' @param hierarchy character vector - hierarchy level to group by (default hierarchy_keys)
#' @param factors character vector - factor columns to include
#' @return a tibble with one row per combination of the `factors` and `hierarchy` columns and, for each other
#'   hierarchy level, a column `<level>_n` with the number of its distinct entries
#'
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$summarize_hierarchy()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
summarize_hierarchy <- function(pdata, hierarchy_keys, hierarchy = hierarchy_keys, factors = character()) {
  .Deprecated("LFQData$summarize_hierarchy()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .summarize_hierarchy
  eval(call, parent.frame())
}

#' Sample subset of proteins/peptides/precursors
#' @param size size of sample
#' @param pdata tidy table
#' @param hierarchy_keys_depth character vector - hierarchy columns at current depth
#' @return a tibble with all rows of `pdata` for `size` randomly drawn hierarchy elements
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_sample()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
sample_subset <- function(size, pdata, hierarchy_keys_depth) {
  .Deprecated("LFQData$get_sample()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .sample_subset
  eval(call, parent.frame())
}

#' Transform long to wide format
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param as.matrix if TRUE return matrix, otherwise data.frame
#' @param file_name if TRUE use file_name as column labels, otherwise sample_name
#' @param sep separator for row IDs when as.matrix = TRUE
#' @param value column name to pivot (default: lfqdata$response())
#' @param row_keys columns identifying a row of the wide table (default: the
#'   hierarchy keys). They must identify one value per sample.
#' @return list with data, rowdata, and annotation (colData)
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$data_wide()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
tidy_to_wide_config <- function(
  lfqdata,
  as.matrix = FALSE,
  file_name = FALSE,
  sep = "~lfq~",
  value = lfqdata$response(),
  row_keys = lfqdata$hierarchy_keys()
) {
  .Deprecated("LFQData$data_wide()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .tidy_to_wide_config
  eval(call, parent.frame())
}

#' Table of distinct factors (sample annotation)
#' @param pdata data.frame
#' @param file_name character - file name column
#' @param sample_name character - sample name column
#' @param factor_keys character vector - factor column names
#' @return a tibble with the distinct combinations of the file name, sample name and factor columns, sorted by
#'   sample name
#'
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$factors()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
table_factors <- function(pdata, file_name, sample_name, factor_keys) {
  .Deprecated("LFQData$factors()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .table_factors
  eval(call, parent.frame())
}

#' Plot a row z-score heatmap with annotations
#'
#' @param matrix numeric matrix - wide-format intensity data
#' @param annotation data.frame - sample annotation
#' @param factor_keys character vector - factor column names for annotation
#' @param sample_name character - sample name column
#' @param na_fraction maximum fraction of missing values per row. Rows that
#'   cannot be meaningfully z-scored are always removed.
#' @param show_rownames if TRUE shows row names, default FALSE
#' @param max_rownames_chars maximum displayed row label length
#' @param max_sample_label_chars maximum displayed sample label length. Labels
#'   keep their suffix because sample prefixes are often shared.
#' @param ... passed to [ComplexHeatmap::Heatmap()]
#' @return invisibly, a ComplexHeatmap `Heatmap` object of the row z-scored data; `NULL` with a warning if
#'   `matrix` has no rows or no row can be z-scored
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Plotter()$heatmap()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
plot_heatmap <- function(
  matrix,
  annotation,
  factor_keys,
  sample_name,
  na_fraction = 0.4,
  show_rownames = FALSE,
  max_rownames_chars = 60,
  max_sample_label_chars = 20,
  ...
) {
  .Deprecated("LFQData$get_Plotter()$heatmap()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .plot_heatmap
  eval(call, parent.frame())
}

#' Plot a sample correlation heatmap with annotations
#'
#' @param matrix numeric matrix - wide-format intensity data
#' @param annotation data.frame - sample annotation
#' @param factor_keys character vector - factor column names for annotation
#' @param sample_name character - sample name column
#' @param R2 logical - plot R-squared instead of correlation
#' @param color colours, spread evenly from the lowest to the highest correlation (default white to red)
#' @param max_sample_label_chars maximum displayed sample label length. Labels
#'   keep their suffix because sample prefixes are often shared.
#' @param method correlation method: `"pearson"` (default), `"spearman"` or `"kendall"` for [stats::cor()] on
#'   pairwise complete observations, or a function of the matrix returning a sample-by-sample correlation matrix
#' @param ... passed to [ComplexHeatmap::Heatmap()]
#' @return invisibly, a ComplexHeatmap `Heatmap` object of the clustered sample correlations (or R^2 values)
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Plotter()$heatmap_cor()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
plot_heatmap_cor <- function(
  matrix,
  annotation,
  factor_keys,
  sample_name,
  R2 = FALSE,
  color = c("white", "red"),
  max_sample_label_chars = 20,
  method = "pearson",
  ...
) {
  .Deprecated("LFQData$get_Plotter()$heatmap_cor()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .plot_heatmap_cor
  eval(call, parent.frame())
}

#' Boxplot of feature intensities by factor levels
#'
#' @param pdata data.frame
#' @param title name to show
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param facet_grid_on on which variable to run facet_grid
#' @param beeswarm use beeswarm default FALSE
#' @param show_mean show mean values
#' @param pb progress bar
#' @return a ggplot object with boxplots of the intensities in `pdata` per combination of factor levels, optionally
#'   with the individual points and the means
#'
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Plotter()$write_boxplots()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
plot_hierarchies_boxplot <- function(
  pdata,
  title,
  lfqdata,
  facet_grid_on = NULL,
  beeswarm = TRUE,
  show_mean = TRUE,
  pb
) {
  .Deprecated("LFQData$get_Plotter()$write_boxplots()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .plot_hierarchies_boxplot
  eval(call, parent.frame())
}

#' Generate feature-level boxplots for all proteins
#' @param pdata data.frame
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param hierarchy e.g. protein_Id default relevant_hierarchy_keys
#' @param facet_grid_on default NULL
#' @return a tibble with one row per element of `hierarchy`, a column (named after the `hierarchy` columns joined
#'   by `+`) identifying it, and a list column `boxplot` of ggplot objects made by [plot_hierarchies_boxplot()]
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Plotter()$boxplots()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
plot_hierarchies_boxplot_df <- function(
  pdata,
  lfqdata,
  hierarchy = lfqdata$relevant_hierarchy_keys(),
  facet_grid_on = NULL
) {
  .Deprecated("LFQData$get_Plotter()$boxplots()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .plot_hierarchies_boxplot_df
  eval(call, parent.frame())
}

#' Visualize intensity distributions
#' @param pdata data.frame
#' @param sample_name character - sample column name
#' @param response character - intensity column name
#' @param is_transformed logical - is intensity log-transformed?
#' @return `plot_intensity_distribution_violin()`: a ggplot object with one violin of the intensities per sample
#'   and the median marked.
#'
#' @import ggplot2
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Plotter()$intensity_distribution_violin()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
plot_intensity_distribution_violin <- function(pdata, sample_name, response, is_transformed = FALSE) {
  .Deprecated("LFQData$get_Plotter()$intensity_distribution_violin()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .plot_intensity_distribution_violin
  eval(call, parent.frame())
}

#' visualize intensity distributions
#' @param pdata data.frame
#' @param sample_name character - sample column name
#' @param response character - intensity column name
#' @param is_transformed logical - is intensity log-transformed?
#' @param legend show legend. If `NA`, hide automatically when the number of
#'   samples is larger than `max_legend_samples`.
#' @param max_legend_samples maximum number of samples for automatic legend
#'   display.
#' @return `plot_intensity_distribution_density()`: a ggplot object with one density line of the intensities per
#'   sample.
#' @rdname plot_intensity_distribution_violin
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Plotter()$intensity_distribution_density()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
plot_intensity_distribution_density <- function(
  pdata,
  sample_name,
  response,
  is_transformed = FALSE,
  legend = NA,
  max_legend_samples = 16
) {
  .Deprecated("LFQData$get_Plotter()$intensity_distribution_density()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .plot_intensity_distribution_density
  eval(call, parent.frame())
}

#' Plot a heatmap of missing values
#' @param matrix numeric matrix - wide-format intensity data
#' @param annotation data.frame - sample annotation
#' @param factor_keys character vector - factor column names for annotation
#' @param sample_name character - sample name column
#' @param limitrows max rows to display
#' @param distance distance method for clustering
#' @return a ComplexHeatmap `Heatmap` object marking the missing values of the rows that have any; `NULL` if fewer
#'   than two such rows
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Plotter()$na_heatmap()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
plot_na_heatmap <- function(matrix, annotation, factor_keys, sample_name, limitrows = 10000, distance = "binary") {
  .Deprecated("LFQData$get_Plotter()$na_heatmap()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .plot_na_heatmap
  eval(call, parent.frame())
}

#' Plot a principal component analysis of the samples
#' @param matrix numeric matrix - wide-format intensity data
#' @param annotation data.frame - sample annotation
#' @param sample_name character - sample name column
#' @param factor_keys character vector - factor column names (first for color, second for shape)
#' @param PC which principal components to plot
#' @param add_txt show sample labels
#' @param nudge label nudge distance
#' @return a ggplot object of the sample scores on the principal components `PC`, coloured by the first and shaped
#'   by the second of `factor_keys`
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Plotter()$pca()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
plot_pca <- function(matrix, annotation, sample_name, factor_keys, PC = c(1, 2), add_txt = FALSE, nudge = 0.1) {
  .Deprecated("LFQData$get_Plotter()$pca()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .plot_pca
  eval(call, parent.frame())
}

#' Plot a heatmap without clustering (use to show NAs)
#' @param matrix numeric matrix - wide-format intensity data
#' @param annotation data.frame - sample annotation
#' @param factor_keys character vector - factor column names for annotation
#' @param sample_name character - sample name column
#' @param arrange either mean or var
#' @param not_na if true than arrange by nr of NA's first and then by arrange
#' @param show_rownames logical, show row names in heatmap
#' @param max_rownames_chars maximum displayed row label length
#' @param max_sample_label_chars maximum displayed sample label length. Labels
#'   keep their suffix because sample prefixes are often shared.
#' @param ... additional arguments passed to [ComplexHeatmap::Heatmap()]
#' @return invisibly, a ComplexHeatmap `Heatmap` object of the abundances with rows ordered by `arrange` and
#'   unclustered columns; `NULL` with a warning if `matrix` has fewer than two rows
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Plotter()$raster()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
plot_raster <- function(
  matrix,
  annotation,
  factor_keys,
  sample_name,
  arrange = c("mean", "var"),
  not_na = FALSE,
  show_rownames = FALSE,
  max_rownames_chars = 60,
  max_sample_label_chars = 20,
  ...
) {
  .Deprecated("LFQData$get_Plotter()$raster()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .plot_raster
  eval(call, parent.frame())
}

#' Visualize correlation among samples
#' @param matrix numeric matrix - wide-format intensity data (samples as columns)
#' @param method correlation method, as for [plot_heatmap_cor()]
#' @param ... passed to [corrplot::corrplot.mixed()], e.g. `tl.cex` and `number.cex` for text size (default 1)
#' @return invisibly, the list returned by [corrplot::corrplot.mixed()] with the correlation matrix `corr` and the
#'   cell positions `corrPos`; called for its side effect of drawing the plot
#'
#' @rdname plot_sample_correlation
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Plotter()$sample_correlation()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
plot_sample_correlation <- function(matrix, method = "pearson", ...) {
  .Deprecated("LFQData$get_Plotter()$sample_correlation()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .plot_sample_correlation
  eval(call, parent.frame())
}

#' Pairs plot with smoothed scatter panels
#' @param data data matrix or data.frame as normally passed to pairs
#' @param legend  add legend to plots
#' @param ... params usually passed to pairs
#' @return `NULL`, invisibly; called for its side effect of drawing the pairs plot
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Plotter()$pairs_smooth()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
pairs_smooth <- function(data, legend = FALSE, ...) {
  .Deprecated("LFQData$get_Plotter()$pairs_smooth()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .pairs_smooth
  eval(call, parent.frame())
}

#' Correlation panel for pairs plots (default lower panel of `pairs_smooth()`)
#' @param x numeric data
#' @param y numeric data
#' @param ... not used
#' @param digits number of digits to display
#' @return `NULL`, invisibly; called for its side effect of writing the correlation, R^2 and p-value into the
#'   current panel
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Plotter()$pairs_smooth()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
panel_cor <- function(x, y, digits = 2, ...) {
  .Deprecated("LFQData$get_Plotter()$pairs_smooth()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .panel_cor
  eval(call, parent.frame())
}

#' Histogram summarizing missingness
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param showempty show empty values
#' @param factors factor columns to use
#' @param alpha transparency
#' @return a ggplot object with density plots of the mean intensity per number of missing values, faceted by
#'   factor levels
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Plotter()$missingness_histogram()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
missingness_histogram <- function(lfqdata, showempty = FALSE, factors = lfqdata$relevant_factor_keys(), alpha = 0.1) {
  .Deprecated("LFQData$get_Plotter()$missingness_histogram()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .missingness_histogram
  eval(call, parent.frame())
}

#' Prepare data for an UpSetR plot of all samples
#'
#' @param lfqdata [LFQData] object with the data and its configuration
#' @return a list with `data`, a data.frame with one row per hierarchy element and one 0/1 presence column per
#'   sample, and `nsets`, the number of sample columns
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Plotter()$upset_missing()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
upset_missing_stats <- function(lfqdata) {
  .Deprecated("LFQData$get_Plotter()$upset_missing()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .upset_missing_stats
  eval(call, parent.frame())
}

#' Compute mean, sd, and CV for all Peptides, or proteins, for all interactions and all samples.
#'
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param factor_key character vector - factor columns to group by (default: relevant_factor_keys)
#' @return `summarize_stats()`: a tibble with one row per hierarchy element and level of
#'   `factor_key`, with the columns `nrReplicates`, `nrMeasured`, `nrNAs`, `sd`, `var`, `meanAbundance`,
#'   `medianAbundance`, `CV` (untransformed data only) and `interaction`.
#' @rdname summarize_stats
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Stats()$stats()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
summarize_stats <- function(lfqdata, factor_key = lfqdata$relevant_factor_keys()) {
  .Deprecated("LFQData$get_Stats()$stats()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .summarize_stats
  eval(call, parent.frame())
}

#' summarize stats output (compute quantiles)
#' @param stats_res result of running `summarize_stats`
#' @param factor_keys_depth character vector - factor columns at current depth
#' @param stats summarize either sd or CV
#' @param probs for which quantiles 10, 20 etc.
#' @return `summarize_stats_quantiles()`: a list with `long`, a tibble with one row per factor level and quantile
#'   (columns `probs` and `quantiles`), and `wide`, the same quantiles with one column per factor level.
#' @rdname summarize_stats
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Stats()$stats_quantiles()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
summarize_stats_quantiles <- function(
  stats_res,
  factor_keys_depth,
  stats = c("sd", "CV"),
  probs = c(0.1, 0.25, 0.5, 0.75, 0.9)
) {
  .Deprecated("LFQData$get_Stats()$stats_quantiles()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .summarize_stats_quantiles
  eval(call, parent.frame())
}

#' Compute theoretical sample sizes from factor level standard deviations
#' @param stats_res data.frame `summarize_stats` output
#' @param delta effect size you are interested in
#' @param power of test
#' @param sig.level P-Value
#' @param min.n smallest n to determine
#' @return a tibble with one row per complete row of `stats_res` and value of `delta`, with the added columns
#'   `delta`, `N_exact` (samples per group needed) and `N` (rounded up)
#'
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Stats()$power_t_test()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
lfq_power_t_test_proteins <- function(stats_res, delta = c(0.59, 1, 2), power = 0.8, sig.level = 0.05, min.n = 1.5) {
  .Deprecated("LFQData$get_Stats()$power_t_test()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .lfq_power_t_test_proteins
  eval(call, parent.frame())
}

#' Estimate sample sizes from standard deviation quantiles
#' @param quantile_sd output of `summarize_stats_quantiles`
#' @param delta effect size you are interested in
#' @param power of test
#' @param sig.level P-Value
#' @param min.n smallest n to determine
#' @return a tibble with the rows of `quantile_sd` repeated for each `delta`, with the added columns `sdtrimmed`
#'   (the sd floored at the smallest sd resolvable with `min.n` samples), `N_exact`, `N` (rounded up) and `delta`
#'
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Stats()$power_t_test_quantiles()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
lfq_power_t_test_quantiles <- function(quantile_sd, delta = c(0.59, 1, 2), power = 0.8, sig.level = 0.05, min.n = 1.5) {
  .Deprecated("LFQData$get_Stats()$power_t_test_quantiles()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .lfq_power_t_test_quantiles
  eval(call, parent.frame())
}

#' Plot the density or ECDF of sd, mean or CV
#' @param pdata data.frame with statistics
#' @param factor_key character - factor column name for colouring
#' @param stat sd, mean or CV
#' @param ggstat either density or ecdf
#' @return a ggplot object with one density or ECDF line of `stat` per level of `factor_key`
#'
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Stats()$density()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
plot_stat_density <- function(pdata, factor_key, stat = c("CV", "meanAbundance", "sd"), ggstat = c("density", "ecdf")) {
  .Deprecated("LFQData$get_Stats()$density()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .plot_stat_density
  eval(call, parent.frame())
}

#' Plot the density or ECDF of sd, mean or CV for intensities below and above the median
#' @param pdata data.frame with statistics
#' @param factor_key character - factor column name for faceting
#' @param stat sd, mean or CV
#' @param ggstat either density or ecdf
#' @return a ggplot object with density or ECDF lines of `stat` for features above and below the median mean
#'   intensity, faceted by `factor_key`
#'
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Stats()$density_median()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
plot_stat_density_median <- function(
  pdata,
  factor_key,
  stat = c("CV", "meanAbundance", "sd"),
  ggstat = c("density", "ecdf")
) {
  .Deprecated("LFQData$get_Stats()$density_median()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .plot_stat_density_median
  eval(call, parent.frame())
}

#' Violin plot of sd, CV or mean
#'
#' @param pdata data.frame with statistics
#' @param factor_keys_depth character vector - factor columns for grouping
#' @param stat either CV, mean or sd
#' @return a ggplot object with one violin of `stat` per combination of factor levels, with the median marked
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Stats()$violin()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
plot_stat_violin <- function(pdata, factor_keys_depth, stat = c("CV", "meanAbundance", "sd")) {
  .Deprecated("LFQData$get_Stats()$violin()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .plot_stat_violin
  eval(call, parent.frame())
}

#' Violin plot of sd, CV or mean for intensities below and above the median
#' @param pdata data.frame with statistics
#' @param factor_key character - factor column name for x-axis
#' @param stat either CV, mean or sd
#' @return a ggplot object with violins of `stat` per level of `factor_key` and the median and interquartile range
#'   marked, faceted into features below and above the median mean intensity
#'
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Stats()$violin_median()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
plot_stat_violin_median <- function(pdata, factor_key, stat = c("CV", "meanAbundance", "sd")) {
  .Deprecated("LFQData$get_Stats()$violin_median()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .plot_stat_violin_median
  eval(call, parent.frame())
}

#' Plot standard deviation against mean to assess variance stability
#' @param pdata data.frame with statistics
#' @param factor_keys_depth character vector - factor columns for faceting
#' @param size how many points to sample (since scatter plot to slow for all)
#' @return a ggplot scatter plot of `sd` against `meanAbundance` with a loess fit, one facet per combination of
#'   factor levels
#'
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Stats()$stdv_vs_mean()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
plot_stdv_vs_mean <- function(pdata, factor_keys_depth, size = 2000) {
  .Deprecated("LFQData$get_Stats()$stdv_vs_mean()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .plot_stdv_vs_mean
  eval(call, parent.frame())
}

#' pooled variance
#' @param res1 per-group statistics from [summarize_stats()]
#' @param config AnalysisConfiguration of the data
#' @rdname pooled_var
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Stats()$power_fdr_pilot()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
poolvar <- function(res1, config) {
  .Deprecated("LFQData$get_Stats()$power_fdr_pilot()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .poolvar
  eval(call, parent.frame())
}

#' Compute pooled variance
#'
#' following the documentation here:
#' https://online.stat.psu.edu/stat500/lesson/7/7.3/7.3.1/7.3.1.1
#'
#' @param x data.frame
#' @return data.frame
#' @rdname pooled_var
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Stats()$power_fdr_pilot()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
compute_pooled <- function(x) {
  .Deprecated("LFQData$get_Stats()$power_fdr_pilot()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .compute_pooled
  eval(call, parent.frame())
}

#' Table of distinct factors with group sizes
#' @param pdata data.frame
#' @param file_name character - file name column
#' @param sample_name character - sample name column
#' @param factor_keys character vector - all factor column names
#' @param factor_keys_depth character vector - factor columns at current depth
#' @return a tibble with one row per combination of the `factor_keys_depth` levels and a column `n` with the
#'   number of samples
#'
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$factors()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
table_factors_size <- function(pdata, file_name, sample_name, factor_keys, factor_keys_depth) {
  .Deprecated("LFQData$factors()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .table_factors_size
  eval(call, parent.frame())
}

#' Summarize missing in condition as barplot
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param factors factor columns to use
#' @return a list with `data`, a wide tibble with one row per condition and one column per number of missing
#'   values holding the feature count, and `figure`, a ggplot bar chart of these counts
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Summariser()$missingness_per_group()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
missingness_per_condition <- function(lfqdata, factors = lfqdata$relevant_factor_keys()) {
  .Deprecated("LFQData$get_Summariser()$missingness_per_group()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .missingness_per_condition
  eval(call, parent.frame())
}

#' Cumulative counts of missing values per condition
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param factors factor columns to use
#' @return a list with `data`, a wide tibble with one row per condition and one column per number of missing
#'   values holding the cumulative feature count, and `figure`, a ggplot bar chart of these counts
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Summariser()$missingness_per_group_cumsum()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
missingness_per_condition_cumsum <- function(lfqdata, factors = lfqdata$relevant_factor_keys()) {
  .Deprecated("LFQData$get_Summariser()$missingness_per_group_cumsum()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .missingness_per_condition_cumsum
  eval(call, parent.frame())
}

#' Robustly scale a matrix by medians and MADs
#' @param data numeric matrix
#' @param dim margin to scale: 2 scales each column (sample), 1 each row
#' @param preserve_mean if TRUE, add the mean of the medians back
#' @return a matrix like `data`, centred by the median of each column (`dim = 2`) or row (`dim = 1`) and divided
#'   by its MAD relative to the mean MAD; only centred, with a warning, if a MAD is zero; with
#'   `preserve_mean = TRUE` the mean of the medians is added back
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Transformer()$robscale()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
robust_scale <- function(data, dim = 2, preserve_mean = FALSE) {
  .Deprecated("LFQData$get_Transformer()$robscale()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .robust_scale
  eval(call, parent.frame())
}

#' Compute the median and MAD of each sample
#' @param lfqdata [LFQData] object with the data and its configuration
#' @return a list with the numeric vectors `medians` and `mads`, with one element per sample
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Transformer()$get_scales()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
get_robscales <- function(lfqdata) {
  .Deprecated("LFQData$get_Transformer()$get_scales()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .lfqdata_robscales
  eval(call, parent.frame())
}

#' Scale data using a subset of the data
#'
#' this should reduce the overall variance.
#'
#' @param lfqdata LFQData object with full dataset
#' @param lfqsubset LFQData object with subset for computing scales
#' @param preserve_mean default FALSE - sets mean to zero
#' @param colname name of the scaled intensity column; default NULL uses \code{<response>_subset_scaled}
#' @return list with data, scales, and colname
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Transformer()$robscale_subset()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
scale_with_subset <- function(lfqdata, lfqsubset, preserve_mean = FALSE, colname = NULL) {
  .Deprecated("LFQData$get_Transformer()$robscale_subset()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .scale_with_subset
  eval(call, parent.frame())
}

#' Apply function requiring a matrix to tidy table
#'
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param .func function taking and returning a matrix
#' @param .funcname name of function (used for creating new column)
#' @return list with `data` (data.frame) and `colname` (new column name)
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Transformer()$intensity_matrix()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
apply_to_response_matrix <- function(lfqdata, .func, .funcname = NULL) {
  .Deprecated("LFQData$get_Transformer()$intensity_matrix()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .apply_to_response_matrix
  eval(call, parent.frame())
}

#' Transform intensity column by applying a function
#' @param pdata data.frame
#' @param response character, name of the response column to transform
#' @param .func function to transform intensities e.g. log2
#' @param .funcname name of function (used for creating new column name)
#' @return list with `data` (data.frame), `colname` (new column name) and `zeros` (logical, rows whose zero
#'   intensity became `NA`)
#' @details For `log`, `log2` and `log10`, zero intensities mean "not detected" and become `NA`; LFQData records
#'   them in its provenance. Negative intensities give `NaN` with a warning, which fitting and aggregation then
#'   reject.
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Transformer()$intensity_array()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
transform_work_intensity <- function(pdata, response, .func, .funcname = NULL) {
  .Deprecated("LFQData$get_Transformer()$intensity_array()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .transform_work_intensity
  eval(call, parent.frame())
}

#' Takes matrix of responses and converts into tibble
#'
#' @param pdata matrix with rownames encoding hierarchy keys
#' @param value name of column to store values in
#' @param config AnalysisConfiguration (needed for column name mapping)
#' @param data optional data.frame to join back to
#' @param sep separator used to unite the hierarchy keys
#' @return a long-format tibble with the hierarchy key, sample name and `value` columns; if `data`
#'   is given, `data` joined with these values
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Transformer()$intensity_matrix()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
response_matrix_as_tibble <- function(pdata, value, config, data = NULL, sep = "~lfq~") {
  .Deprecated("LFQData$get_Transformer()$intensity_matrix()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .response_matrix_as_tibble
  eval(call, parent.frame())
}

#' Aggregates e.g. protein abundances from peptide abundances
#'
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param method "medpolish" (Tukey's median polish, column `medpolish`) or "rlm" (robust regression with
#'   `MASS::rlm`, column `lmrob`)
#' @param control convergence controls merged over the defaults. medpolish: `maxiter` (default 10) and `eps`
#'   (default 0.01, the change of the sum of absolute residuals that stops the polish, on the log2 scale of the
#'   input); rlm: `maxit` (default 20) and `acc` (default 1e-4, relative change of the IRLS residuals).
#' @return returns list with data (data.frame), config (AnalysisConfiguration) and failures (tibble with the
#'   subject `id` and `warning` text of every warning raised while estimating one subject; the warnings are not
#'   re-raised)
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Aggregator("medpolish")$aggregate()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
estimate_intensity <- function(lfqdata, method = c("medpolish", "rlm"), control = list()) {
  .Deprecated("LFQData$get_Aggregator(\"medpolish\")$aggregate()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .estimate_intensity
  eval(call, parent.frame())
}

#' Aggregates top N intensities
#'
#' run \link{rank_peptide_by_intensity} first
#' @param ranked_data data.frame with ranked peptides
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param .func function to use for aggregation
#' @param N default 3 top intensities.
#' @return list with data and new reduced configuration (config)
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Aggregator("topN")$aggregate()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
aggregate_intensity_top_n <- function(ranked_data, lfqdata, .func, N = 3) {
  .Deprecated("LFQData$get_Aggregator(\"topN\")$aggregate()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .aggregate_intensity_top_n
  eval(call, parent.frame())
}

#' Rank precursors or peptides by intensity
#' @param pdata data.frame
#' @param response character - intensity column name
#' @param hierarchy_keys character vector - all hierarchy column names
#' @return data.frame
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Aggregator("topN")$aggregate()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
rank_peptide_by_intensity <- function(pdata, response, hierarchy_keys) {
  .Deprecated("LFQData$get_Aggregator(\"topN\")$aggregate()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .rank_peptide_by_intensity
  eval(call, parent.frame())
}

#' Plot feature data and result of aggregation
#'
#' @param lfqdata LFQData object (original data)
#' @param lfqdata_agg LFQData object (aggregated data)
#' @param show.legend logical, show legend in plot
#' @return a tibble with one row per aggregated element and the columns `hierarchy_id`, `data` (the feature
#'   intensities), `other` (the aggregated intensities) and `plots` (ggplot objects of the feature intensities with
#'   the aggregate drawn in black)
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Aggregator()$plot()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
plot_estimate <- function(lfqdata, lfqdata_agg, show.legend = FALSE) {
  .Deprecated("LFQData$get_Aggregator()$plot()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .plot_estimate
  eval(call, parent.frame())
}

#' Plot peptide intensities of protein as a function of the sample and factor
#'
#' @param res data.frame
#' @param protein_name title of plot
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param show.legend logical, show legend in plot
#' @return a ggplot object with one line per feature across samples, coloured by peptide, faceted by factor
#'   levels and titled `protein_name`
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Aggregator()$plot()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
plot_hierarchies_line <- function(res, protein_name, lfqdata, show.legend = FALSE) {
  .Deprecated("LFQData$get_Aggregator()$plot()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .plot_hierarchies_line
  eval(call, parent.frame())
}

#' Make reduced hierarchy configuration
#' @param config AnalysisConfiguration
#' @param work_intensity work intensity column
#' @param hierarchy new reduced hierarchy
#' @return AnalysisConfiguration with reduced hierarchy
#' @note Deprecated in prolfqua 1.9.5; use \code{LFQData$get_Aggregator()$aggregate()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
make_reduced_hierarchy_config <- function(config, work_intensity, hierarchy) {
  .Deprecated("LFQData$get_Aggregator()$aggregate()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .make_reduced_hierarchy_config
  eval(call, parent.frame())
}

#' Compute contrasts for possibly incomplete models
#'
#' only keeps non NA coefficients.
#'
#' @param m linear model generated using lm
#' @param linfct linear function
#' @param confint confidence interval default 0.95
#' @return a data.frame with one row per row of `linfct` and columns `lhs`, `sigma`, `df`, `estimate`,
#'   `std.error`, `statistic`, `p.value`, `conf.low` and `conf.high`; the estimates are `NA` for rows that use a
#'   coefficient the model could not estimate
#'
#' @note Deprecated in prolfqua 1.9.5; use \code{build_contrast_analysis()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
compute_contrast <- function(m, linfct, confint = 0.95) {
  .Deprecated("build_contrast_analysis()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .compute_estimable_contrast
  eval(call, parent.frame())
}

#' Compute contrasts of a mixed-effects model with lmerTest
#' @param model mixed effects model
#' @param linfct linear function
#' @param ddf method to determine denominator degrees of freedom
#' @return a tibble with one row per row of `linfct` and columns `lhs`, `estimate`, `std.error`, `df`,
#'   `statistic`, `conf.low`, `conf.high`, `p.value` and `sigma`; `NA` (with a warning) if `linfct` has no rows or
#'   the model is rank-deficient
#' @note Deprecated in prolfqua 1.9.5; use \code{build_contrast_analysis()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
compute_lmer_contrast <- function(model, linfct, ddf = c("Satterthwaite", "Kenward-Roger")) {
  .Deprecated("build_contrast_analysis()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .compute_lmer_contrast
  eval(call, parent.frame())
}

#' Compute group averages or contrasts for each model
#'
#' used in p2621, p2109
#'
#' @param models model table with one fit per subject in the column `linear_model`, e.g. the `model_df` of
#'   [model_analyse()]
#' @param linfct linear functions: a matrix with one row per contrast and one column per model coefficient, or
#'   a list with one such matrix per model
#' @param subject_id columns identifying a subject
#' @param contrastfun function computing the contrasts of one fit, e.g. [compute_lmer_contrast()]
#' @return a tibble with one row per subject and row of `linfct`, holding the `subject_id` columns, the columns
#'   returned by `contrastfun`, and `isSingular`, `sigma.model` and `df.residual.model`; subjects whose contrasts
#'   failed are dropped with a warning
#' @note Deprecated in prolfqua 1.9.5; use \code{build_contrast_analysis()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
contrasts_linfct <- function(models, linfct, subject_id = "protein_Id", contrastfun = .compute_lmer_contrast) {
  .Deprecated("build_contrast_analysis()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .contrasts_linfct
  eval(call, parent.frame())
}

#' Get group-level linear functions from a fitted model
#'
#' @details The fitted model supplies its coefficient names and observed
#'   categorical factor combinations. Each combination is expressed as weights
#'   on the fitted coefficients; weights for an individual factor level are
#'   averaged equally over the observed combinations containing that level.
#'   Thus \code{linfct_factors} describes marginal factor-level means and
#'   \code{linfct_interactions} describes cell means. With treatment coding,
#'   a reference-level mean includes the intercept; a non-reference mean also
#'   includes its treatment coefficient. The \emph{Under the hood} article
#'   shows a complete \code{~ A * B} example. These helpers assume categorical
#'   predictors and the usual intercept-based coding; they are not a general
#'   marginal-means implementation for continuous predictors or arbitrary coding.
#' @param as_list if TRUE, return the factor and interaction matrices as a list; if FALSE, as one matrix
#' @param m linear model
#' @return if `as_list = TRUE`, a list with the matrices `linfct_factors` (one row per factor level) and
#'   `linfct_interactions` (one row per observed factor combination), each with one column per model coefficient;
#'   otherwise the two matrices bound by rows into a single matrix
#' @note Deprecated in prolfqua 1.9.5; use \code{Contrasts$get_linfct()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
linfct_from_model <- function(m, as_list = TRUE) {
  .Deprecated("Contrasts$get_linfct()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .linfct_from_model
  eval(call, parent.frame())
}

#' Convert named group-mean contrasts to coefficient weights
#'
#' @details Rows of \code{linfct} represent group or interaction means and
#'   columns represent fitted coefficients. Contrast expressions combine these
#'   rows, so \code{"group_A - group_Ctrl"} becomes their row difference.
#'   Quote interaction names containing colons with backticks. Expressions may
#'   reference earlier named contrasts. A row naming a group unavailable in a
#'   particular model is dropped with a warning. The returned matrix has one
#'   row per estimable contrast and coefficient columns in model order.
#'
#' @param linfct linear functions as created by linfct_from_model
#' @param contrasts named character vector of contrasts to determine linear functions for
#' @param p.message print messages default FALSE
#' @return a matrix with one row per computable contrast and one column per model coefficient
#'
#' @note Deprecated in prolfqua 1.9.5; use \code{Contrasts$get_linfct()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
linfct_matrix_contrasts <- function(linfct, contrasts, p.message = FALSE) {
  .Deprecated("Contrasts$get_linfct()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .linfct_matrix_contrasts
  eval(call, parent.frame())
}

#' Retrieve complete model fits
#' @param complete_models model table, the `model_df` of [model_analyse()]
#' @return the rows of `complete_models` that have a model fit, the maximal number of estimated coefficients and
#'   more than one residual degree of freedom
#' @note Deprecated in prolfqua 1.9.5; use \code{build_contrast_analysis()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
get_complete_model_fit <- function(complete_models) {
  .Deprecated("build_contrast_analysis()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .get_complete_model_fit
  eval(call, parent.frame())
}

#' Create an ANOVA extractor returning a data frame
#' @param test one case-sensitive test choice: "F", "Chisq", "LRT", "Rao" or "Cp".
#'   The fitted model's ANOVA method must support the chosen test.
#' @return An \code{\link{AnovaExtractor}} object.
#' @note Deprecated in prolfqua 1.9.5; use \code{Model$get_anova()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
get_anova_df <- function(test = "F") {
  .Deprecated("Model$get_anova()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .get_anova_df
  eval(call, parent.frame())
}

#' Fit one model per subject with a modelling strategy
#'
#' used in project p2901
#'
#' @param pepIntensity long-format data with one row per observation
#' @param model_strategy modelling strategy, e.g. from [strategy_lm()]
#' @param subject_id columns identifying a subject; one model is fitted per subject
#' @param model_name name stored with the fitted models
#' @param label short label naming this fitting pass in progress output;
#'   defaults to `model_name`.
#' @param progress progress reporter selector passed to [.make_progress()];
#'   `NULL` (default) keeps the legacy terminal `progress::progress_bar`, a
#'   `function(i, total, label)` or `"message"` emits batch-visible progress.
#'   Defaults to `getOption("prolfqua.progress", NULL)`.
#' @return a list with `model_name` and `model_df`, a tibble with one row per subject holding the nested `data`,
#'   the fitted `linear_model` (the error message for a failed fit), `fit_warnings`, `has_model_fit`,
#'   `isSingular`, `df.residual`, `sigma`, `nr_coef`, `nr_coef_not_NA`, `converged` and `iterations`
#' @note Deprecated in prolfqua 1.9.5; use \code{build_contrast_analysis()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
model_analyse <- function(
  pepIntensity,
  model_strategy,
  subject_id = "protein_Id",
  model_name = "Model",
  label = model_name,
  progress = getOption("prolfqua.progress", NULL)
) {
  .Deprecated("build_contrast_analysis()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .model_analyse
  eval(call, parent.frame())
}

#' Summarize modelling and error reporting
#' @param mod model table see \code{\link{build_model}}
#' @return a list with the tables `exists` (counts of subjects with and without a model fit) and `isSingular`
#'   (counts of singular and non-singular fits)
#' @note Deprecated in prolfqua 1.9.5; use \code{build_contrast_analysis()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
model_summary <- function(mod) {
  .Deprecated("build_contrast_analysis()", package = "prolfqua")
  res <- list()
  res$exists <- table(mod$model_df$has_model_fit)
  res$isSingular <- table(mod$model_df$isSingular)
  res
}

#' Fit Firth logistic regression to binary detection data
#'
#' Fits Firth's logistic regression per subject to the binary detection
#' response of aggregated data, one row per subject and sample.
#' @param data aggregated LFQData with a binary detection response
#' @param formula model formula with the binary response
#' @param control convergence controls, see \code{\link{strategy_logistf}}
#' @return a [ModelFirth] object holding one Firth logistic model per subject
#' @note Deprecated in prolfqua 1.9.5; use \code{build_contrast_analysis(method = "firth")}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
build_model_logistf <- function(data, formula, control = list()) {
  .Deprecated("build_contrast_analysis(method = \"firth\")", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .build_model_logistf
  eval(call, parent.frame())
}

#' Encode intensities as a binary detection response
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param name column name for binary response
#' @return a tibble of the long-format data with an added integer column `name`: 1 if the response is observed,
#'   0 if it is `NA`
#' @note Deprecated in prolfqua 1.9.5; use \code{build_contrast_analysis(method = "firth")}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
encode_bin_resp <- function(lfqdata, name = "bin_resp") {
  .Deprecated("build_contrast_analysis(method = \"firth\")", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .encode_bin_resp
  eval(call, parent.frame())
}

#' Create an interaction column from factors
#' @param data data frame
#' @param columns factor columns to combine
#' @param sep separator between the levels in the combined labels
#' @return `data` with an added factor column `interaction` combining the values of `columns`, each prefixed with
#'   its column name and joined by `sep`
#' @note Deprecated in prolfqua 1.9.5; use \code{build_contrast_analysis()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
make_interaction_column <- function(data, columns, sep = ".") {
  .Deprecated("build_contrast_analysis()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .make_interaction_column
  eval(call, parent.frame())
}

#' Moderate p-values - limma approach
#' @param contrast_df contrast table with the columns `sigma`, the residual degrees of freedom and the estimate
#' @param df name of the residual degrees-of-freedom column
#' @param estimate name of the effect-estimate column
#' @param robust passed to [limma::squeezeVar()]: estimate the prior robustly
#' @param confint confidence level of the moderated intervals
#' @param variance_floor optional positive lower bound for posterior variances
#' @return `contrast_df`, ungrouped, with the added columns `moderated.var.post`, `moderated.var.prior`,
#'   `moderated.df.prior`, `moderated.std.error`, `moderated.statistic`, `moderated.df.total`,
#'   `moderated.p.value`, `moderated.conf.low` and `moderated.conf.high`
#' @note Deprecated in prolfqua 1.9.5; use \code{build_contrast_analysis()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
moderated_p_limma <- function(
  contrast_df,
  df = "df",
  estimate = "diff",
  robust = FALSE,
  confint = 0.95,
  variance_floor = NULL
) {
  .Deprecated("build_contrast_analysis()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .moderated_p_limma
  eval(call, parent.frame())
}

#' DEqMS-style count-dependent variance moderation
#'
#' Applies count-dependent empirical Bayes variance shrinkage to a contrast
#' result table. Proteins quantified from many peptides get less shrinkage;
#' proteins from few peptides get more.
#'
#' @param mm data.frame from one contrast group with columns: sigma, df,
#'   statistic, std.error, and the estimate column
#' @param count_col name of column with peptide/PSM count per protein
#' @param df name of the degrees of freedom column
#' @param estimate name of the fold change column
#' @param loess_span span parameter for LOESS fit (default 0.75)
#' @param confint confidence level for intervals (default 0.95)
#' @return data.frame with added moderated.* columns
#' @note Deprecated in prolfqua 1.9.5; use \code{build_contrast_analysis(method = "deqms")}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
moderated_p_deqms <- function(mm, count_col, df = "df", estimate = "diff", loess_span = 0.75, confint = 0.95) {
  .Deprecated("build_contrast_analysis(method = \"deqms\")", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .moderated_p_deqms
  eval(call, parent.frame())
}

#' DEqMS-style moderation for long contrast table
#'
#' Splits by contrast group and applies \code{\link{moderated_p_deqms}} to
#' each group independently.
#'
#' @param mm contrast result data.frame (long format)
#' @param count_col name of column with peptide/PSM count per protein
#' @param group_by_col column to group by (default "contrast")
#' @param estimate name of the fold change column
#' @param loess_span span parameter for LOESS fit (default 0.75)
#' @return data.frame with moderated columns added
#' @note Deprecated in prolfqua 1.9.5; use \code{build_contrast_analysis(method = "deqms")}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
moderated_p_deqms_long <- function(mm, count_col, group_by_col = "contrast", estimate = "diff", loess_span = 0.75) {
  .Deprecated("build_contrast_analysis(method = \"deqms\")", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .moderated_p_deqms_long
  eval(call, parent.frame())
}

#' Compute protein-level fold changes and p-values using the beta distribution
#'
#' takes p-value of the scaled p-value
#'
#' @param contrasts_data data frame
#' @param contrast name of column with contrast identifier
#' @param subject_id name of column with typically protein Id
#' @param estimate name of column with effect size estimate
#' @param statistic statistic name of column with statistic (typically t-statistics)
#' @param p.value name of column with moderated.p.value
#' @param max.n used to limit the number of peptides in probability computation.
#' @return data.frame with columns
#'
#'
#' @note Deprecated in prolfqua 1.9.5; use \code{build_contrast_analysis(method = "ropeca_nested")}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
summary_ropeca_median_p_scaled <- function(
  contrasts_data,
  contrast = "contrast",
  subject_id = "protein_Id",
  estimate = "diff",
  statistic = "statistic",
  p.value = "moderated.p.value",
  max.n = 10
) {
  .Deprecated("build_contrast_analysis(method = \"ropeca_nested\")", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .summary_ropeca_median_p_scaled
  eval(call, parent.frame())
}

#' Protein p-value from the p-value of the median fold-change peptide
#' @param median.p.value median of the peptide p-values of a protein
#' @param n.obs number of peptides
#' @param max.n limit number of peptides per protein.
#' @return a numeric vector of p-values, the beta-distribution probability that the median of
#'   `min(n.obs, max.n)` uniform p-values is at most `median.p.value`
#' @note Deprecated in prolfqua 1.9.5; use \code{build_contrast_analysis(method = "ropeca_nested")}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
get_p_values_pbeta <- function(median.p.value, n.obs, max.n = 10) {
  .Deprecated("build_contrast_analysis(method = \"ropeca_nested\")", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .get_p_values_pbeta
  eval(call, parent.frame())
}

#' Volcano plot for multiple contrasts
#' @param data data in long format
#' @param effect column containing effect sizes
#' @param significance column containing p-values, q.values etc
#' @param contrast column with contrast
#' @param colour colouring of points
#' @param xintercept fc thresholds
#' @param yintercept significance threshold
#' @param label column containing labels
#' @param size label text size in mm (default 3, about 8.5 pt)
#' @param segment.size controls size of lines
#' @param segment.alpha controls visibility of lines
#' @param scales parameter to ggplot2::facet_wrap
#' @param max_nr_significant_text maximum number of significant labels to display
#' @param .data deprecated; use `data`
#' @return a ggplot object with one volcano panel per contrast; if `label` is given, significant points are labelled
#' @note Deprecated in prolfqua 1.9.5; use \code{ModelInterface$coef_volcano()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
multigroup_volcano <- function(
  data,
  effect = "fc",
  significance = "p.adjust",
  contrast = "condition",
  colour = "colour",
  xintercept = c(-2, 2),
  yintercept = 0.05,
  label = NULL,
  size = 3,
  segment.size = 0.3,
  segment.alpha = 0.3,
  scales = "fixed",
  max_nr_significant_text = 20,
  .data
) {
  .Deprecated("ModelInterface$coef_volcano()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .multigroup_volcano_labelled
  eval(call, parent.frame())
}

#' Build linear models on simulated data for testing
#' @param model design: `"parallel2"`, `"parallel3"`, `"factors"` or `"interaction"`
#' @param n_proteins number of simulated proteins
#' @param with_missing include simulated missing measurements
#' @param weight_missing weight controlling simulated missingness
#' @param Nprot deprecated; use `n_proteins`
#' @return a [Model] object with one linear model per simulated protein, fitted with [build_model()]
#' @note Deprecated in prolfqua 1.9.5; use \code{sim_lfq_data_protein_config() with build_contrast_analysis()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
sim_build_models_lm <- function(
  model = c("parallel2", "parallel3", "factors", "interaction"),
  n_proteins = 10,
  with_missing = TRUE,
  weight_missing = 1,
  Nprot
) {
  .Deprecated("sim_lfq_data_protein_config() with build_contrast_analysis()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .sim_build_models_lm
  eval(call, parent.frame())
}

#' Build mixed-effects models on simulated data for testing
#' @inheritParams sim_build_models_lm
#' @return a [Model] object with one mixed-effects model per simulated protein, fitted with [build_model()]
#' @note Deprecated in prolfqua 1.9.5; use \code{sim_lfq_data_peptide_config() with build_contrast_analysis()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
sim_build_models_lmer <- function(
  model = c("parallel2", "parallel3", "factors", "interaction"),
  n_proteins = 10,
  with_missing = TRUE,
  weight_missing = 1,
  Nprot
) {
  .Deprecated("sim_lfq_data_peptide_config() with build_contrast_analysis()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .sim_build_models_lmer
  eval(call, parent.frame())
}

#' Build Firth logistic models on simulated data for testing
#' @inheritParams sim_build_models_lm
#' @return a [ModelFirth] object fitted with [build_model_logistf()]
#' @note Deprecated in prolfqua 1.9.5; use \code{sim_lfq_data_protein_config() with build_contrast_analysis()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
sim_build_models_logistf <- function(
  model = c("parallel2", "parallel3", "factors", "interaction"),
  n_proteins = 10,
  with_missing = TRUE,
  weight_missing = 1,
  Nprot
) {
  .Deprecated("sim_lfq_data_protein_config() with build_contrast_analysis()", package = "prolfqua")
  if (!missing(Nprot)) {
    n_proteins <- .deprecate_arg("Nprot", "n_proteins", Nprot)
  }
  model <- match.arg(model)
  if (model != "parallel3") {
    istar <- prolfqua::sim_lfq_data_2factor_config(
      n_proteins = n_proteins,
      with_missing = with_missing,
      weight_missing = weight_missing
    )
  } else {
    istar <- prolfqua::sim_lfq_data_protein_config(
      n_proteins = n_proteins,
      with_missing = with_missing,
      weight_missing = weight_missing
    )
  }
  istar$data <- encode_bin_resp(LFQData$new(istar$data, istar$config))
  istar$config$bin_resp <- "bin_resp"
  istar <- prolfqua::LFQData$new(istar$data, istar$config)

  model <- if (model == "factors") {
    "~ Treatment + Background"
  } else if (model == "interaction") {
    "~ Treatment * Background"
  } else if (model == "parallel2") {
    "~ Treatment"
  } else if (model == "parallel3") {
    "~ group_"
  } else {
    NULL
  }
  model_function <- paste0(istar$get_config()$bin_resp, model)
  mod <- build_model_logistf(
    istar,
    model_function
  )
  mod
}

#' Fit a linear model to simulated data for examples
#' @param model design: `"parallel2"`, `"parallel3"`, `"factors"` or `"interaction"`
#' @return an `lm` object fitted to one simulated protein
#' @note Deprecated in prolfqua 1.9.5; use \code{sim_lfq_data_protein_config() with build_contrast_analysis()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
sim_make_model_lm <- function(model = c("parallel2", "parallel3", "factors", "interaction")) {
  .Deprecated("sim_lfq_data_protein_config() with build_contrast_analysis()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .sim_make_model_lm
  eval(call, parent.frame())
}

#' Fit a mixed-effects model to simulated data for examples
#' @param model design: `"parallel2"`, `"parallel3"`, `"factors"` or `"interaction"`
#' @return an `lmerModLmerTest` object, the first successful fit among ten simulated proteins
#' @note Deprecated in prolfqua 1.9.5; use \code{sim_lfq_data_peptide_config() with build_contrast_analysis()}. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
sim_make_model_lmer <- function(model = c("parallel2", "parallel3", "factors", "interaction")) {
  .Deprecated("sim_lfq_data_peptide_config() with build_contrast_analysis()", package = "prolfqua")
  call <- match.call()
  call[[1L]] <- .sim_make_model_lmer
  eval(call, parent.frame())
}

#' Copy script and other files from a package to a working directory
#' @param runscripts paths of the files relative to the package root (or \code{inst/})
#' @param workdir directory to copy the files to
#' @param packagename package to copy the files from
#' @param overwrite overwrite existing files in \code{workdir}, default TRUE
#' @return paths of the copied files
#' @note Deprecated in prolfqua 1.9.5; use \code{prolfquapp::script_copy_helper_vec()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
script_copy_helper_vec <-
  function(runscripts, workdir = getwd(), packagename = "prolfqua", overwrite = TRUE) {
    .Deprecated("prolfquapp::script_copy_helper_vec()", package = "prolfqua")
    res <- NULL
    for (scripts in runscripts) {
      src_script <- find_package_file(packagename, scripts)
      dest_script <- file.path(workdir, basename(scripts))
      inform_prolfqua("copy ", src_script, " to ", dest_script)
      if (is.null(src_script)) {
        warning(sprintf("Could not copy script file %s.", dest_script), call. = FALSE)
      } else if (file.copy(src_script, dest_script, overwrite = overwrite)) {
        res <- c(res, dest_script)
      } else {
        warning(sprintf("Could not copy script file %s to %s.", src_script, dest_script), call. = FALSE)
      }
    }
    inform_prolfqua(sprintf("your working directory now should contain: %d new files:", length(res)))
    res
  }

#' Separates hierarchy columns into starting columns
#'
#'
#' @param data data.frame
#' @param config AnalysisConfiguration
#' @return `data` with each hierarchy key at the current depth split at `config$sep` into its component columns,
#'   which are added after the key column
#' @note Deprecated in prolfqua 1.9.5; use \code{prolfquapp::separate_hierarchy()}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
separate_hierarchy <- function(data, config) {
  .Deprecated("prolfquapp::separate_hierarchy()", package = "prolfqua")
  for (hkey in config$hierarchy_keys_depth()) {
    if (length(config$hierarchy[[hkey]]) == 1 && hkey == config$hierarchy[[hkey]]) {} else {
      data <- data |>
        tidyr::separate_wider_delim(
          dplyr::all_of(hkey),
          delim = config$sep,
          names = config$hierarchy[[hkey]],
          too_few = "align_start",
          too_many = "drop",
          cols_remove = FALSE
        ) |>
        dplyr::relocate(dplyr::all_of(hkey), .before = dplyr::all_of(config$hierarchy[[hkey]][1]))
    }
  }
  data
}

#' Compute fold changes given Contrasts
#' @param data data.frame
#' @param hierarchy_keys hierarchy key columns of the AnalysisConfiguration
#' @param contrasts list of contrasts to compute
#' @return a tibble with one row per hierarchy element and contrast and columns of the `hierarchy_keys`,
#'   `group_1`, `group_2`, `estimate`, `group_1_name`, `group_2_name` and `contrast`
#' @note Deprecated in prolfqua 1.9.5 and will be removed; only prolfquatakeoutmissing used it. See
#'   \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
get_contrast <- function(data, hierarchy_keys, contrasts) {
  .Deprecated(
    msg = paste(
      "'get_contrast' is deprecated and will be removed from prolfqua;",
      "only prolfquatakeoutmissing used it, which now has its own copy."
    )
  )
  res <- vector(mode = "list", length(contrasts))
  names(res) <- names(contrasts)
  for (i in seq_along(contrasts)) {
    cname <- names(contrasts)[i]

    # A contrast must be a difference "LHS - RHS" (e.g. "group_A - group_B" or
    # "(group_A + group_B)/2 - group_Ctrl"). group_1 / group_2 are the evaluated
    # left / right sides; estimate is their difference.
    sx <- .contrast_sides_expr(contrasts[i])
    if (is.null(sx)) {
      abort_prolfqua(
        "get_contrast: contrast '",
        cname,
        "' (",
        contrasts[i],
        ") is not of the required form 'LHS - RHS'. A contrast must be a ",
        "difference of two group expressions, e.g. 'group_A - group_B' or ",
        "'(group_A + group_B)/2 - group_Ctrl'."
      )
    }

    # A contrast referencing a level absent from `data` errors here at mutate
    # ("object '<level>' not found") -- once per contrast, not per row.
    dd <- dplyr::mutate(data, group_1 = !!sx$lhs, group_2 = !!sx$rhs, estimate = !!sx$full)
    # Materialize this contrast under its own name so a *later* nested contrast
    # can reference it as a group, e.g. "T_C_gv_WT - T_C_gv_KO" after
    # T_C_gv_WT / T_C_gv_KO have been defined. Sides may therefore be a level
    # present in `data` or the name of an earlier contrast.
    data[[cname]] <- dd[["estimate"]]
    df <- dplyr::select(dd, dplyr::all_of(hierarchy_keys), "group_1", "group_2", "estimate")
    df$group_1_name <- rlang::as_label(sx$lhs)
    df$group_2_name <- rlang::as_label(sx$rhs)
    df$contrast <- cname
    res[[cname]] <- df
  }
  dplyr::ungroup(dplyr::bind_rows(res))
}

#' Convert a matrix or data.frame to a tibble (taken from tidyquant)
#'
#' @param x a matrix
#' @param preserve_row_names give name to rownames column, if NULL discard rownames
#' @param ... further parameters passed to as_tibble
#' @return The computed result.
#' @note Deprecated in prolfqua 1.9.5; use \code{tibble::as_tibble(x, rownames = )}. See \link{prolfqua-deprecated}.
#' @keywords internal
#' @export
matrix_to_tibble <- function(x, preserve_row_names = "row.names", ...) {
  .Deprecated("tibble::as_tibble(x, rownames = )", package = "prolfqua")
  if (!(is.matrix(x) || is.data.frame(x))) {
    abort_bad_argument("x", "be a matrix or data.frame")
  }
  if (!is.null(preserve_row_names)) {
    row.names <- rownames(x)
    if (!is.null(row.names)) {
      dplyr::bind_cols(
        tibble::tibble(!!preserve_row_names := row.names),
        tibble::as_tibble(x, ...)
      )
    } else {
      warning("No row names to preserve. Object otherwise converted to tibble successfully.", call. = FALSE)
      tibble::as_tibble(x, ...)
    }
  } else {
    tibble::as_tibble(x, ...)
  }
}
