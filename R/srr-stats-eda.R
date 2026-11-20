#' EDA statistical software standards
#'
#' Standards 0.3.0, Exploratory Data Analysis. Assessment per standard in srr/standards-compliance.md.
#'
#' @srrstatsTODO {EA1.0} target audience not stated explicitly; block H, TODO/prolfqua/PLAN_srr_block_H_eda.md.
#' @srrstats {EA1.1} Long-format mass-spectrometry LFQ data (intensities per protein/peptide/precursor and sample)
#'   described in DESCRIPTION and vignettes/prolfqua.Rmd.
#' @srrstats {EA1.2} Questions (QC, variance and CV, sample size, differential expression) stated in the
#'   QualityControlAndSampleSizeEstimation and prolfqua vignettes.
#' @srrstats {EA1.3} Each exported function documents its input as LFQData, AnalysisConfiguration or data.frame in
#'   @param.
#' @srrstats {EA2.0} AnalysisConfiguration names the index columns (hierarchy, file_name, sample_name, factors) used
#'   for all joins and reshapes.
#' @srrstatsTODO {EA2.1} LFQData$new(setup = FALSE, the default) skips the check; block H,
#'   TODO/prolfqua/PLAN_srr_block_H_eda.md.
#' @srrstats {EA2.2} Index columns are identified through the AnalysisConfiguration R6 class carried by every LFQData.
#' @srrstats {EA2.2a} AnalysisConfiguration (R/AnalysisConfiguration.R) is the class system for index columns.
#' @srrstats {EA2.3} Joins use explicit `by =` keys taken from the configuration or function arguments, never assumed
#'   column names.
#' @srrstatsTODO {EA2.6} no handling/tests of vectors with extra attributes; block H,
#'   TODO/prolfqua/PLAN_srr_block_H_eda.md.
#' @srrstats {EA3.0} build_contrast_analysis() runs one of about 20 methods on all proteins at once; LFQDataStats
#'   reports variance, CV and sample size for all proteins.
#' @srrstats {EA3.1} All facades return the same contrast schema (ContrastsInterface column descriptions);
#'   merge_contrasts_results() and the ComparingMethods vignette compare methods.
#' @srrstatsTODO {EA4.0} data.frame input returns tibble; untested; block H, TODO/prolfqua/PLAN_srr_block_H_eda.md.
#' @srrstatsTODO {EA4.1} no digits/precision parameter; block H, TODO/prolfqua/PLAN_srr_block_H_eda.md.
#' @srrstatsTODO {EA4.2} print/plot/summary for LFQData and facades; block H, TODO/prolfqua/PLAN_srr_block_H_eda.md.
#' @srrstatsTODO {EA5.0} no accessibility consideration; block H, TODO/prolfqua/PLAN_srr_block_H_eda.md.
#' @srrstatsTODO {EA5.0a} small label text by default; block H, TODO/prolfqua/PLAN_srr_block_H_eda.md.
#' @srrstatsTODO {EA5.0b} colour-blind-unsafe defaults; block H, TODO/prolfqua/PLAN_srr_block_H_eda.md.
#' @srrstatsTODO {EA5.1} overrides ignore accessibility; block H, TODO/prolfqua/PLAN_srr_block_H_eda.md.
#' @srrstatsTODO {EA5.2} most output relies on default tibble formatting; block H,
#'   TODO/prolfqua/PLAN_srr_block_H_eda.md.
#' @srrstatsTODO {EA5.3} summary statistics do not state storage mode explicitly; block H,
#'   TODO/prolfqua/PLAN_srr_block_H_eda.md.
#' @srrstats {EA5.4} Axes use ggplot2 pretty breaks; PCA axis labels round the explained variance
#'   (R/tidyMS_plotting.R).
#' @srrstatsTODO {EA5.5} no units on axes; block H, TODO/prolfqua/PLAN_srr_block_H_eda.md.
#' @srrstatsTODO {EA6.0} Closed with EA6.0b and EA6.0d; block H, TODO/prolfqua/PLAN_srr_block_H_eda.md.
#' @srrstats {EA6.0a} Tests check classes and types of returned objects (expect_s3_class, expect_type, expect_r6_class
#'   throughout tests/testthat).
#' @srrstatsTODO {EA6.0b} many return values untested for dimensions; block H, TODO/prolfqua/PLAN_srr_block_H_eda.md.
#' @srrstats {EA6.0c} Tests check column names of returned tables (e.g. test-ContrastsFacades.R, test-LFQDataStats.R).
#' @srrstatsTODO {EA6.0d} tests of column classes/types; block H, TODO/prolfqua/PLAN_srr_block_H_eda.md.
#' @srrstats {EA6.0e} Single values tested with expect_equal() and explicit tolerance (e.g. test-Contrasts.R).
#' @srrstatsTODO {EA6.1} no vdiffr snapshots; most plots only class-checked; block H,
#'   TODO/prolfqua/PLAN_srr_block_H_eda.md.
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
