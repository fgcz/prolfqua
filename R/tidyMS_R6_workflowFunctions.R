#' correlation preprocessing
#'
#' @export
#' @family workflow functions
#' @examples
#' library(LFQService)
#' rm(list=ls())
#' config <- spectronautDIAData250_config$clone(deep=T)
#' config$parameter$min_nr_of_notNA  <- 20
#' data <- spectronautDIAData250_analysis
#' res <- workflow_correlation_preprocessing(data,config)
#' names(res)
#'
workflow_correlation_preprocessing <- function(data, config, minCorrelation = 0.7){
  stat_input <- hierarchyCounts(data, config)

  data_NA <- removeLarge_Q_Values(data, config)
  data_NA <- summariseQValues(data_NA, config)

  data_NA_QVal <- data_NA %>% filter_at( "srm_QValueMin" , all_vars(. < config$parameter$qValThreshold )   )

  stat_qval <- hierarchyCounts(data_NA_QVal, config)

  # remove transitions with large numbers of NA's
  data_NA_QVal <- rankPrecursorsByNAs(data_NA_QVal, config)
  data_NA_QVal <- data_NA_QVal %>% dplyr::filter(srm_NrNotNAs > config$parameter$min_nr_of_notNA)
  stat_min_nr_of_notNA <- hierarchyCounts(data_NA_QVal, config)

  # remove single hit wonders
  data_NA_QVal <- nr_B_in_A(data_NA_QVal,config$table$hierarchyKeys()[1], config$table$hierarchyKeys()[2])
  c_name <- make_name(config$table$hierarchyKeys()[1], config$table$hierarchyKeys()[2])
  data_NA_QVal <- dplyr::filter(data_NA_QVal, !!sym(c_name) >= config$parameter$min_peptides_protein )

  stat_min_peptides_protein  <- hierarchyCounts(data_NA_QVal, config)

  # filter decorrelated.
  data_NA_QVal <- transform_work_intensity(data_NA_QVal, config, log2)
  data_NA_QVal <- markDecorrelated(data_NA_QVal, config, minCorrelation = minCorrelation)
  keepCorrelated <- dplyr::filter(data_NA_QVal, srm_decorelated == FALSE)

  stat_correlated  <- hierarchyCounts(keepCorrelated, config)

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
  x <- bind_rows(stats)
  stats <- add_column(x, processing = names(stats),.before = 1)


  return(list(data = proteinIntensities$data, stats = stats, newconfig = proteinIntensities$newconfig))
}

#' apply correlation filtering and impuation
#' @export
#' @examples
#' rm(list=ls())
#' library(tidyverse)
#' library(LFQService)
#' config <- spectronautDIAData250_config$clone(deep=T)
#' config$parameter$min_nr_of_notNA  <- 20
#' data <- spectronautDIAData250_analysis
#' res <- workflow_corr_filter_impute(data,config)
workflow_corr_filter_impute <- function(data,config, minCorrelation =0.6){
  stat_input <- hierarchyCounts(data, config)

  data_NA <- removeLarge_Q_Values(data, config)
  data_NA <- summariseQValues(data_NA, config)

  data_NA_QVal <- data_NA %>% filter_at( "srm_QValueMin" , all_vars(. < config$parameter$qValThreshold )   )

  stat_qval <- hierarchyCounts(data_NA_QVal, config)

  # remove transitions with large numbers of NA's
  data_NA_QVal <- rankPrecursorsByNAs(data_NA_QVal, config)
  data_NA_QVal <- data_NA_QVal %>% dplyr::filter(srm_NrNotNAs > config$parameter$min_nr_of_notNA)
  stat_min_nr_of_notNA <- hierarchyCounts(data_NA_QVal, config)

  # remove single hit wonders
  data_NA_QVal <- nr_B_in_A(data_NA_QVal,config$table$hierarchyKeys()[1], config$table$hierarchyKeys()[2])
  c_name <- make_name(config$table$hierarchyKeys()[1], config$table$hierarchyKeys()[2])
  data_NA_QVal <- dplyr::filter(data_NA_QVal, !!sym(c_name) >= config$parameter$min_peptides_protein )

  stat_min_peptides_protein  <- hierarchyCounts(data_NA_QVal, config)

  # filter decorrelated.
  data_NA_QVal <- transform_work_intensity(data_NA_QVal, config, log2)
  data_NA_QVal <- markDecorrelated(data_NA_QVal, config, minCorrelation = minCorrelation)
  keepCorrelated <- dplyr::filter(data_NA_QVal, srm_decorelated == FALSE)

  stat_correlated  <- hierarchyCounts(keepCorrelated, config)
  keepCorrelated <- rankPrecursorsByIntensity(keepCorrelated, config)
  qvalFiltImputed <- impute_correlationBased(keepCorrelated, config)


  stats <- list(stat_input = stat_input,
                stat_qval = stat_qval,
                stat_min_nr_of_notNA = stat_min_nr_of_notNA,
                stat_min_peptides_protein = stat_min_peptides_protein,
                stat_correlated = stat_correlated
  )
  x <- bind_rows(stats)
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
#' config <- spectronautDIAData250_config$clone(deep=T)
#' data <- spectronautDIAData250_analysis
#' hierarchyCounts(data, config)
#' tmp <-workflow_NA_preprocessing(data, config)
#' tmp <-workflow_NA_preprocessing(data, config, percent=70)
#' config$get
#' hierarchyCounts(tmp, config)
#' stopifnot(FALSE==(is.grouped_df(tmp)))
workflow_NA_preprocessing <- function(data, config, percent = 60, factor_level = 1){
  stat_input <- hierarchyCounts(data, config)
  data_NA <- removeLarge_Q_Values(data, config)
  data_NA <- summariseQValues(data_NA, config)

  data_NA_QVal <- data_NA %>% filter_at( "srm_QValueMin" , all_vars(. < config$parameter$qValThreshold )   )

  stat_qval <- hierarchyCounts(data_NA_QVal, config)

  resNACondition <- proteins_WithXPeptidesInCondition(data_NA_QVal,
                                                      config,
                                                      percent = percent,
                                                      factor_level = factor_level)
  data_NA_QVal_condition <- inner_join(resNACondition, data_NA_QVal )

  # Complete cases
  data_NA_QVal_condition <- completeCases( data_NA_QVal_condition , config)
  return(data_NA_QVal_condition)
}

#' Get Protein Intensities after QValue NA filtering and medianpolish
#'
#' @export
#' @family workflow functions
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' rm(list=ls())
#' config <- spectronautDIAData250_config$clone(deep=T)
#' data <- spectronautDIAData250_analysis
#' tmp <-workflow_Q_NA_filtered_Hierarchy(data, config, hierarchy_level=2)
#' nrow(tmp$data)
#' config <- spectronautDIAData250_config$clone(deep=T)
#' res <-workflow_Q_NA_filtered_Hierarchy(data, config, hierarchy_level=1)
#' nrow(res$data)
#' res$newconfig
#' hierarchyCounts(res$data, res$newconfig)
workflow_Q_NA_filtered_Hierarchy <- function(data,
                                             config,
                                             percent = 60,
                                             hierarchy_level=1,
                                             factor_level=1){
  stat_input <- hierarchyCounts(data, config)
  data_NA <- removeLarge_Q_Values(data, config)
  data_NA <- summariseQValues(data_NA, config)
  data_NA_QVal <- data_NA %>% filter_at( "srm_QValueMin" , all_vars(. < config$parameter$qValThreshold )   )
  stat_qval <- hierarchyCounts(data_NA_QVal, config)
  resNACondition <- proteins_WithXPeptidesInCondition(data_NA_QVal,
                                                      config,
                                                      percent =percent,
                                                      factor_level = factor_level)
  data_NA_QVal_condition <- inner_join(resNACondition, data_NA_QVal )

  resDataLog <- LFQService::transform_work_intensity(data_NA_QVal_condition , config, log2)

  resDataLog <- applyToIntensityMatrix(resDataLog, config, robust_scale)
  figs3 <- applyToHierarchyBySample(resDataLog, config, medpolishPly, hierarchy_level = hierarchy_level)
  colnames(figs3$data[[1]])
  colnames(figs3$medpolishPly[[1]])

  protIntensity <- figs3 %>% select(config$table$hierarchyKeys()[1:hierarchy_level], medpolishPly) %>% unnest()

  newconfig <- make_reduced_hierarchy_config(config,
                                             workIntensity = "medpolish",
                                             hierarchy = config$table$hierarchy[1:hierarchy_level])
  return(list(data = protIntensity, newconfig = newconfig))
}

