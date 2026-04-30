#' @importFrom forcats fct_relevel
#' @importFrom dplyr across all_of anti_join arrange bind_cols bind_rows case_when count desc distinct filter group_by
#' @importFrom dplyr inner_join left_join mutate nest_by
#' @importFrom dplyr rename select starts_with summarize ungroup
#'
#' @importFrom ggplot2 aes element_text facet_grid ggplot ggtitle geom_boxplot
#' @importFrom ggplot2 geom_line geom_violin guides
#' @importFrom ggplot2 stat_summary scale_x_continuous scale_y_continuous theme
#'
#' @importFrom ggrepel geom_text_repel
#' @importFrom graphics par rect hist
#' @importFrom graphics abline pairs plot text
#' @importFrom grDevices colorRampPalette dev.off pdf png rainbow
#' @importFrom pheatmap pheatmap
#' @importFrom plotly ggplotly
#' @importFrom purrr map map2 map2_dbl map_lgl map_chr map_dbl reduce map_if map_dfc map_int map_df
#' @importFrom rlang := sym syms .data
#' @importFrom stats as.formula cor
#' @importFrom stats .lm.fit glm model.matrix residuals
#' @importFrom tidyr nest nesting separate_rows unite unnest separate
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom tibble add_column as_tibble column_to_rownames tibble
#' @importFrom stats anova coef coefficients confint cor.test delete.response
#'   df.residual ecdf formula lm mad median medpolish
#' @importFrom stats na.omit p.adjust pbeta power.t.test prcomp predict pt qt
#'   quantile sd setNames sigma terms update vcov
#' @importFrom stats approx df fisher.test ftable lm.fit pchisq pf pnorm qf rgeom rlnorm rnorm uniroot
#' @importFrom utils combn data head read.csv tail unzip
#' @importFrom limma lmFit contrasts.fit eBayes topTable
NULL


## @importFrom vsn justvsn

# Suppress R CMD check NOTEs for variables used in NSE (dplyr/ggplot2)
utils::globalVariables(c(
  ".",
  ".env",
  "Freq",
  "Group",
  "PC",
  "Replicate",
  "group",
  "isSingular",
  "linear_model",
  "meanAbundance",
  "nr_peptides",
  "percent_variance_explained",
  "reference_mean",
  "reference_median",
  "subject",
  "y"
))

#' Internal Functions by category
#' @family aggregation
#' @family benchmarking
#' @family concrete_configuration
#' @family configuration
#' @family deprecated
#' @family MaxQuant
#' @family modelling
#' @family FragPipe
#' @family plotting
#' @family preprocessing
#' @family stats
#' @family summary
#' @family transitionCorrelation
#' @family utilities
#' @family vignetteHelpers
#' @family workflows
#' @name INTERNAL_FUNCTIONS_BY_FAMILY
#'
#' @return A documentation topic.
NULL

#' Package Data
#' @family data
#' @name PACKAGE_DATA
#' @return A documentation topic.
NULL


#' SAINT express output
#'
#' SAINT express output produced by running the function
#' runSaint
#'
#' @family data
#' @keywords internal
#' @docType data
#' @format a linear model
data_SAINTe_output <- NULL

#' example dataset for testing
#' @family data
#' @docType data
#' @keywords internal
#' @format a data frame
x5463yzwer453bbb <- NULL


#' Benchmark data Example
#' @format A data frame
#' @family data
#' @docType data
#' @keywords internal
data_benchmarkExample <- NULL

#' example data for check
#' @family data
#' @docType data
#' @keywords internal
data_checksummarizationrobust87 <- NULL

#' example data for check
#' @family data
#' @docType data
#' @keywords internal
data_checksummarizerobust <- NULL

#' example data for check
#' @family data
#' @docType data
#' @keywords internal
#'
data_checksummarizerobust69 <- NULL

#' example data for check
#' @family data
#' @docType data
#' @keywords internal
#'
data_correlatedPeptideList <- NULL

#' example data for check
#' @family data
#' @docType data
#' @keywords internal
#'
data_IonstarProtein_subsetNorm <- NULL

#' example data for check
#' @family data
#' @docType data
#' @keywords internal
data_ionstar <- NULL

#' example data for check of scores produced based on confusion matrix
#' @family data
#' @docType data
#' @keywords internal
data_test_confusion_matrix_scores <- NULL

#' example data for check
#' @family data
#' @docType data
#' @keywords internal
data_skylinePRMSample_A <- NULL

#' example data for check
#' @family data
#' @docType data
#' @keywords internal
data_skylineSRM_HL_A <- NULL

#' example data for check
#' @family data
#' @docType data
#' @keywords internal
data_spectronautDIA250_A <- NULL
