#' @importFrom conflicted conflict_prefer
#' @importFrom forcats fct_relevel
#' @importFrom dplyr all_of anti_join arrange bind_cols bind_rows case_when count desc distinct filter group_by group_by_at
#' @importFrom dplyr inner_join left_join mutate nest_by
#' @importFrom dplyr one_of rename select select_at summarize_at summarize ungroup vars
#'
#' @importFrom ggplot2 aes_string element_text facet_grid ggplot ggtitle geom_boxplot
#' @importFrom ggplot2 geom_line geom_violin guides
#' @importFrom ggplot2 stat_summary scale_x_continuous scale_y_continuous theme
#'
#' @importFrom ggrepel geom_text_repel
#' @importFrom GGally ggpairs
#' @importFrom glue glue
#' @importFrom graphics par rect hist
#' @importFrom graphics abline pairs plot text
#' @importFrom grDevices colorRampPalette dev.off pdf png rainbow
#' @importFrom htmlwidgets saveWidget
#' @importFrom pheatmap pheatmap
#' @importFrom plotly ggplotly
#' @importFrom purrr map map2 map2_dbl map_lgl map_chr map_dbl reduce map_if map_dfc map_int map_df
#' @importFrom rlang := UQ sym syms .data
#' @importFrom stats as.formula cor
#' @importFrom stats .lm.fit glm model.matrix residuals
#' @importFrom tidyr gather nest nesting separate_rows spread unite unnest unnest_legacy separate
#' @importFrom tidyr pivot_wider
#' @importFrom tibble add_column as_tibble column_to_rownames tibble
#' @importFrom stats anova coef coefficients confint cor.test df.residual ecdf formula lm mad median medpolish
#' @importFrom stats na.omit p.adjust pbeta power.t.test prcomp predict pt qt quantile sd setNames sigma terms update vcov
#' @importFrom stringr str_trim
#' @import tidyverse
#' @importFrom utils combn data head read.csv tail unzip
#' @importFrom yaml write_yaml
NULL


## @importFrom vsn justvsn
## @importFrom prozor readPeptideFasta

#' Internal Functions by category
#' @family configuration
#' @family concrete_configuration
#' @family summary
#' @family aggregation
#' @family plotting
#' @family utilities
#' @family modelling
#' @family workflows
#' @family vignetteHelpers
#' @family deprecated
#' @family transitionCorrelation
#' @family preprocessing
#' @family benchmarking
#' @family MaxQuant
#' @family MSFragger
#' @family stats
#' @name INTERNAL_FUNCTIONS_BY_FAMILY
#'
NULL

#' Package Data
#' @family data
#' @name PACKAGE_DATA
NULL

# data_basicModel_p1807 <- prolfqua::basicModel_p1807
# usethis::use_data(data_basicModel_p1807, compress = TRUE, overwrite = TRUE)
# file.remove("data/basicModel_p1807.rda")

#' Model for debugging
#' @family data
#' @keywords internal
#' @docType data
#' @format a linear model
"data_basicModel_p1807"

# data_benchmarkExample <- prolfqua::benchmarkDataExample
# usethis::use_data(data_benchmarkExample, compress = TRUE, overwrite = TRUE)
# file.remove("data/benchmarkDataExample.rda")

#' Benchmark data Example
#' @format A data frame
#' @family data
#' @docType data
#' @keywords internal
"data_benchmarkExample"

# data_checksummarizationrobust87 <- prolfqua::checksummarizationrobust87
# usethis::use_data(data_checksummarizationrobust87, compress = TRUE, overwrite = TRUE)
# file.remove("data/checksummarizationrobust87.rda")

#' example data for check
#' @family data
#' @docType data
#' @keywords internal
'data_checksummarizationrobust87'

# data_checksummarizerobust <- prolfqua::checksummarizerobust
# usethis::use_data(data_checksummarizerobust, compress = TRUE, overwrite = TRUE)
# file.remove("data/checksummarizerobust.rda")

#' example data for check
#' @family data
#' @docType data
#' @keywords internal
'data_checksummarizerobust'

# data_checksummarizerobust69 <- prolfqua::checksummarizerobust69
# usethis::use_data(data_checksummarizerobust69, compress = TRUE, overwrite = TRUE)
# file.remove("data/checksummarizerobust69.rda")

#' example data for check
#' @family data
#' @docType data
#' @keywords internal
#'
'data_checksummarizerobust69'

#data_correlatedPeptideList <- prolfqua::correlatedPeptideList
#usethis::use_data(data_correlatedPeptideList, compress = TRUE, overwrite = TRUE)
#file.remove("data/correlatedPeptideList.rda")

#' example data for check
#' @family data
#' @docType data
#' @keywords internal
#'
'data_correlatedPeptideList'

# data_IonstarProtein_subsetNorm <- prolfqua::dataIonstarProtein_subsetNorm
# usethis::use_data(data_IonstarProtein_subsetNorm, compress = TRUE, overwrite = TRUE)
# file.remove("data/dataIonstarProtein_subsetNorm.rda")

#' example data for check
#' @family data
#' @keywords internal
#'
'data_IonstarProtein_subsetNorm'


#' example data for check
#' @family data
#' @keywords internal
'data_exampleForRopeca'

# data_factor_levelContrasts <- factor_levelContrasts
# usethis::use_data(data_factor_levelContrasts, compress = TRUE, overwrite = TRUE)
# file.remove("data/factor_levelContrasts.rda")

#' example data for check
#' @family data
#' @keywords internal
'data_factor_levelContrasts'

# data_interactionModel_p1807 <- interactionModel_p1807
# usethis::use_data(data_interactionModel_p1807, compress = TRUE, overwrite = TRUE)
# file.remove("data/interactionModel_p1807.rda")

#' example data for check
#' @family data
#' @keywords internal
'data_interactionModel_p1807'

#data_ionstar <- ionstar
#usethis::use_data(data_ionstar, compress = TRUE, overwrite = TRUE)
#file.remove("data/ionstar.rda")

#' example data for check
#' @family data
#' @keywords internal
'data_ionstar'

#data_lm_models_to_test <- lm_models_to_test
#usethis::use_data(data_lm_models_to_test, compress = TRUE, overwrite = TRUE)
#file.remove("data/lm_models_to_test.rda")

#' example data for check
#' @family data
#' @keywords internal
'data_lm_models_to_test'

# data_modellingResult_A <- modellingResult_A
# usethis::use_data(data_modellingResult_A, compress = TRUE, overwrite = TRUE)
# file.remove("data/modellingResult_A.rda")

#' example data for check
#' @family data
#' @keywords internal
'data_modellingResult_A'

#data_models_interaction <- models_interaction
#usethis::use_data(data_models_interaction,compress = TRUE, overwrite = TRUE)
#file.remove("data/models_interaction.rda")

#' example data for check
#' @family data
#' @keywords internal
'data_models_interaction'

# data_multigroupFC <- multigroupFCDATA
# usethis::use_data(data_multigroupFC,compress = TRUE, overwrite = TRUE)
# file.remove("data/multigroupFCDATA.rda")

#' example data for check
#' @family data
#' @keywords internal
'data_multigroupFC'


# data_skylinePRMSample_A <- skylinePRMSampleData_A
# usethis::use_data(data_skylinePRMSample_A,compress = TRUE, overwrite = TRUE)
# file.remove("data/skylinePRMSampleData_A.rda")

#' example data for check
#' @family data
#' @keywords internal
'data_skylinePRMSample_A'


# data_skylineSRM_HL_A <- skylineSRM_HL_A
# usethis::use_data(data_skylineSRM_HL_A,compress = TRUE, overwrite = TRUE)
# file.remove("data/skylineSRM_HL_A.rda")

#' example data for check
#' @family data
#' @keywords internal
'data_skylineSRM_HL_A'


#data_spectronautDIA250_A <- spectronautDIAData250_A
#usethis::use_data(data_spectronautDIA250_A,compress = TRUE, overwrite = TRUE)
#file.remove("data/spectronautDIAData250_A.rda")

#' example data for check
#' @family data
#' @keywords internal
'data_spectronautDIA250_A'


#' example of 2 Factor data
#' @family data
#' @format LFQData
'data_Yeast2Factor'

