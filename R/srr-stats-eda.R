#' EDA statistical software standards
#'
#' Standards 0.3.0, Exploratory Data Analysis. Assessment per standard in srr/standards-compliance.md.
#'
#' @srrstats {EA1.1} Long-format mass-spectrometry LFQ data (intensities per protein/peptide/precursor and sample)
#'   described in DESCRIPTION and vignettes/prolfqua.Rmd.
#' @srrstats {EA1.2} Questions (QC, variance and CV, sample size, differential expression) stated in the
#'   QualityControlAndSampleSizeEstimation and prolfqua vignettes.
#' @srrstats {EA1.3} Each exported function documents its input as LFQData, AnalysisConfiguration or data.frame in
#'   @param.
#' @srrstats {EA2.0} AnalysisConfiguration names the index columns (hierarchy, file_name, sample_name, factors) used
#'   for all joins and reshapes.
#' @srrstats {EA2.2} Index columns are identified through the AnalysisConfiguration R6 class carried by every LFQData.
#' @srrstats {EA2.2a} AnalysisConfiguration (R/AnalysisConfiguration.R) is the class system for index columns.
#' @srrstats {EA2.3} Joins use explicit `by =` keys taken from the configuration or function arguments, never assumed
#'   column names.
#' @srrstats {EA3.0} build_contrast_analysis() runs one of about 20 methods on all proteins at once; LFQDataStats
#'   reports variance, CV and sample size for all proteins.
#' @srrstats {EA3.1} All facades return the same contrast schema (ContrastsInterface column descriptions);
#'   merge_contrasts_results() and the ComparingMethods vignette compare methods.
#' @srrstats {EA5.4} Axes use ggplot2 pretty breaks; PCA axis labels round the explained variance
#'   (R/tidyMS_plotting.R).
#' @srrstats {EA6.0a} Tests check classes and types of returned objects (expect_s3_class, expect_type, expect_r6_class
#'   throughout tests/testthat).
#' @srrstats {EA6.0c} Tests check column names of returned tables (e.g. test-ContrastsFacades.R, test-LFQDataStats.R).
#' @srrstats {EA6.0e} Single values tested with expect_equal() and explicit tolerance (e.g. test-Contrasts.R).
#' @srrstats {EA5.5} Mass-spectrometry intensities carry no physical unit, so none can be obtained from the
#'   input. Axis titles name the response column, which records the transformation (e.g. `log2_abundance`);
#'   untransformed intensities are drawn on a log10 axis; effect axes name the effect and score columns.
#' @noRd
NULL

#' NA_standards
#'
#' @srrstatsNA {EA2.2b} Index columns are identified through a class system (EA2.2a), not table attributes.
#' @srrstatsNA {EA2.4} Input is one long table plus its configuration; no multi-table input exists.
#' @srrstatsNA {EA2.5} As EA2.4: a single long table, whose index columns follow EA2.0-EA2.3.
#' @srrstatsNA {EA5.6} No JavaScript library is bundled; interactive plots use the plotly package.
#' @noRd
NULL
