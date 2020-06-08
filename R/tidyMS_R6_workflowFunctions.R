#' correlation preprocessing
#'
#' @export
#' @family workflow functions
#' @examples
#' library(LFQService)
#' rm(list=ls())
#' config <- LFQServiceData::spectronautDIAData250_config$clone(deep=TRUE)
#' config$parameter$min_nr_of_notNA  <- 20
#' data <- LFQServiceData::spectronautDIAData250_analysis
#' res <- workflow_correlation_preprocessing(data,config)
#' names(res)
#'
workflow_correlation_preprocessing_protein_intensities <- function(data, config, minCorrelation = 0.7){
  stat_input <- hierarchy_counts(data, config)

  data_NA_QVal <- filter_byQValue(data, config)
  stat_qval <- hierarchy_counts(data_NA_QVal, config)

  # remove transitions with large numbers of NA's
  data_NA_QVal <- rankPrecursorsByNAs(data_NA_QVal, config)
  data_NA_QVal <- data_NA_QVal %>% dplyr::filter(srm_NrNotNAs > config$parameter$min_nr_of_notNA)
  stat_min_nr_of_notNA <- hierarchy_counts(data_NA_QVal, config)

  # remove single hit wonders
  data_NA_QVal <- nr_B_in_A(data_NA_QVal,config$table$hierarchyKeys()[1], config$table$hierarchyKeys()[2])
  c_name <- make_name(config$table$hierarchyKeys()[1], config$table$hierarchyKeys()[2])
  data_NA_QVal <- dplyr::filter(data_NA_QVal, !!sym(c_name) >= config$parameter$min_peptides_protein )

  stat_min_peptides_protein  <- hierarchy_counts(data_NA_QVal, config)

  # filter decorrelated.
  data_NA_QVal <- transform_work_intensity(data_NA_QVal, config, log2)
  data_NA_QVal <- markDecorrelated(data_NA_QVal, config, minCorrelation = minCorrelation)
  keepCorrelated <- dplyr::filter(data_NA_QVal, srm_decorelated == FALSE)

  stat_correlated  <- hierarchy_counts(keepCorrelated, config)

  # TODO check if you are not aggregating log transformed intensities
  # rank precursors by intensity
  keepCorrelated <- rankPrecursorsByIntensity(keepCorrelated, config)
  qvalFiltImputed <- impute_correlationBased(keepCorrelated, config)
  mean_na <- function(x){mean(x, na.rm = TRUE)}
  proteinIntensities <- aggregateTopNIntensities(qvalFiltImputed, config, func = mean_na,N=3)

  # collect stats
  stats <- list(stat_input = stat_input,
                stat_qval = stat_qval,
                stat_min_nr_of_notNA = stat_min_nr_of_notNA,
                stat_min_peptides_protein = stat_min_peptides_protein,
                stat_correlated = stat_correlated
  )
  x <- dplyr::bind_rows(stats)
  stats <- add_column(x, processing = names(stats),.before = 1)


  return(list(data = proteinIntensities$data, stats = stats, newconfig = proteinIntensities$newconfig))
}
#' Deprectated
#' @export
workflow_correlation_preprocessing <-function(data, config, minCorrelation = 0.7){
  warning("this function name is deprecated, use workflow_correlation_preprocessing_protein_intensities instead.")
  workflow_correlation_preprocessing_protein_intensities(data, config, minCorrelation)
}

#' apply correlation filtering and impute missing values
#' @export
#' @examples
#' rm(list=ls())
#' library(tidyverse)
#' library(LFQService)
#' config <- LFQServiceData::spectronautDIAData250_config$clone(deep=TRUE)
#' config$parameter$min_nr_of_notNA  <- 20
#' data <- LFQServiceData::spectronautDIAData250_analysis
#' res <- workflow_corr_filter_impute(data,config)
workflow_corr_filter_impute <- function(data,config, minCorrelation =0.6){
  stat_input <- hierarchy_counts(data, config)

  data_NA_QVal <- filter_byQValue(data, config)
  stat_qval <- hierarchy_counts(data_NA_QVal, config)

  # remove transitions with large numbers of NA's
  data_NA_QVal <- rankPrecursorsByNAs(data_NA_QVal, config)
  data_NA_QVal <- data_NA_QVal %>% dplyr::filter(srm_NrNotNAs > config$parameter$min_nr_of_notNA)
  stat_min_nr_of_notNA <- hierarchy_counts(data_NA_QVal, config)

  # remove single hit wonders
  data_NA_QVal <- nr_B_in_A(data_NA_QVal,config$table$hierarchyKeys()[1], config$table$hierarchyKeys()[2])
  c_name <- make_name(config$table$hierarchyKeys()[1], config$table$hierarchyKeys()[2])
  data_NA_QVal <- dplyr::filter(data_NA_QVal, !!sym(c_name) >= config$parameter$min_peptides_protein )

  stat_min_peptides_protein  <- hierarchy_counts(data_NA_QVal, config)

  # filter decorrelated.
  data_NA_QVal <- transform_work_intensity(data_NA_QVal, config, log2)
  data_NA_QVal <- markDecorrelated(data_NA_QVal, config, minCorrelation = minCorrelation)
  keepCorrelated <- dplyr::filter(data_NA_QVal, srm_decorelated == FALSE)

  stat_correlated  <- hierarchy_counts(keepCorrelated, config)
  keepCorrelated <- rankPrecursorsByIntensity(keepCorrelated, config)
  qvalFiltImputed <- impute_correlationBased(keepCorrelated, config)


  stats <- list(stat_input = stat_input,
                stat_qval = stat_qval,
                stat_min_nr_of_notNA = stat_min_nr_of_notNA,
                stat_min_peptides_protein = stat_min_peptides_protein,
                stat_correlated = stat_correlated
  )
  x <- dplyr::bind_rows(stats)
  stats <- add_column(x, processing = names(stats),.before = 1)
  return(qvalFiltImputed)
}

#' filter QVAlues and NA's and factor information
#'
#' @export
#' @family workflow functions
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' rm(list=ls())
#' config <- LFQServiceData::spectronautDIAData250_config$clone(deep=TRUE)
#' data <- LFQServiceData::spectronautDIAData250_analysis
#' hierarchy_counts(data, config)
#' tmp <-workflow_DIA_NA_preprocessing(data, config)
#' hierarchy_counts(tmp$data, config)
#' tmp <-workflow_DIA_NA_preprocessing(data, config, percent=70)
#' hierarchy_counts(tmp$data, config)
#' stopifnot(FALSE==(is.grouped_df(tmp$data)))
workflow_DIA_NA_preprocessing <- function(data,
                                          config,
                                          percent = 60,
                                          hierarchy_level = 2,
                                          factor_level = 1,
                                          min_peptides_protein = config$parameter$min_peptides_protein)
{
  stat_input <- hierarchy_counts(data, config)

  data_NA_QVal <- filter_byQValue(data, config)
  stat_qval <- hierarchy_counts(data_NA_QVal, config)

  resNACondition <- filter_factor_levels_by_missing(data_NA_QVal,
                                                    config,
                                                    percent = percent)

  stat_naFilter <- hierarchy_counts(resNACondition, config)
  protID <- summarize_hierarchy(resNACondition,config) %>%
    dplyr::filter(!!sym(paste0(config$table$hierarchyKeys()[hierarchy_level],"_n"))
                  >= min_peptides_protein)

  data_NA_QVal_condition <- protID %>%
    dplyr::select(config$table$hierarchyKeys()[1]) %>%
    dplyr::inner_join(resNACondition)

  # Complete cases
  data_NA_QVal_condition <- complete_cases( data_NA_QVal_condition , config)
  stat_peptidFitler <- hierarchy_counts(data_NA_QVal_condition, config)
  stats = list(stat_input=stat_input,
               stat_qval = stat_qval,
               stat_naFilter = stat_naFilter,
               stat_peptidFitler = stat_peptidFitler
  )
  return(list(data=data_NA_QVal_condition,stats=stats))
}



#' Get Protein Intensities after QValue NA filtering and medianpolish
#'
#' @export
#' @family workflow functions
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' rm(list=ls())
#' config <- LFQServiceData::spectronautDIAData250_config$clone(deep=TRUE)
#' data <- LFQServiceData::spectronautDIAData250_analysis
#' tmp <-workflow_DIA_Q_NA_filtered_medpolish_protein_intensities(data, config, hierarchy_level=2)
#' nrow(tmp$data)
#' config <- LFQServiceData::spectronautDIAData250_config$clone(deep=TRUE)
#' res <-workflow_DIA_Q_NA_filtered_medpolish_protein_intensities(data, config, hierarchy_level=1)
#' nrow(res$data)
#' hierarchy_counts(res$data, res$config)
workflow_DIA_Q_NA_filtered_medpolish_protein_intensities <- function(data,
                                                           config,
                                                           percent = 60,
                                                           hierarchy_level=1,
                                                           factor_level=1){

  data_NA_QVal_condition <- workflow_DIA_NA_preprocessing(data, config=config,
                                                          percent=percent,
                                                          hierarchy_level = 2,
                                                          factor_level = factor_level
  )$data
  resDataLog <- LFQService::transform_work_intensity(data_NA_QVal_condition , config, log2)
  resDataLog <- applyToIntensityMatrix(resDataLog, config, robust_scale)
  figs3 <- intensity_summary_by_hkeys(resDataLog, config, medpolishPly)
  figs3 <- figs3("unnest") # retrieve default result
  return(figs3)
}

#' Deprecated
#' @export
workflow_Q_NA_filtered_Hierarchy <- function(data,
                                             config,
                                             percent = 60,
                                             hierarchy_level=1,
                                             factor_level=1)
{
  warning("this function name is deprecated, use workflow_DIA_Q_NA_filtered_medpolish_protein_intensities instead.")
  workflow_DIA_Q_NA_filtered_medpolish_protein_intensities(data,
                                                 config,
                                                 percent = 60,
                                                 hierarchy_level=1,
                                                 factor_level=1)
}

