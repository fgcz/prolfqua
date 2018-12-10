#' correlation preprocessing
#'
#' @export
#' @family workflow functions
#' @examples
#' library(LFQService)
#' rm(list=ls())
#' config <- spectronautDIAData250_config$clone(deep=TRUE)
#' config$parameter$min_nr_of_notNA  <- 20
#' data <- spectronautDIAData250_analysis
#' res <- workflow_correlation_preprocessing(data,config)
#' names(res)
#'
workflow_correlation_preprocessing <- function(data, config, minCorrelation = 0.7){
  stat_input <- hierarchyCounts(data, config)

  data_NA_QVal <- filter_byQValue(data, config)
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
#' config <- spectronautDIAData250_config$clone(deep=TRUE)
#' config$parameter$min_nr_of_notNA  <- 20
#' data <- spectronautDIAData250_analysis
#' res <- workflow_corr_filter_impute(data,config)
workflow_corr_filter_impute <- function(data,config, minCorrelation =0.6){
  stat_input <- hierarchyCounts(data, config)

  data_NA_QVal <- filter_byQValue(data, config)


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
#' config <- spectronautDIAData250_config$clone(deep=TRUE)
#' data <- spectronautDIAData250_analysis
#' hierarchyCounts(data, config)
#' tmp <-workflow_DIA_NA_preprocessing(data, config)
#' hierarchyCounts(tmp$data, config)
#' tmp <-workflow_DIA_NA_preprocessing(data, config, percent=70)
#' hierarchyCounts(tmp$data, config)
#' stopifnot(FALSE==(is.grouped_df(tmp$data)))
workflow_DIA_NA_preprocessing <- function(data,
                                          config,
                                          percent = 60,
                                          hierarchy_level = 2,
                                          factor_level = 1,
                                          min_peptides_protein = config$parameter$min_peptides_protein)
{
  stat_input <- hierarchyCounts(data, config)

  data_NA_QVal <- filter_byQValue(data, config)
  stat_qval <- hierarchyCounts(data_NA_QVal, config)

  resNACondition <- filter_factor_levels_by_missing(data_NA_QVal,
                                             config,
                                             percent = percent,
                                             factor_level = factor_level)

  stat_naFilter <- hierarchyCounts(resNACondition, config)
  protID <- summarizeHierarchy(resNACondition,config) %>%
    dplyr::filter(!!sym(paste0(config$table$hierarchyKeys()[hierarchy_level],"_n"))
                  >= min_peptides_protein)

  data_NA_QVal_condition <- protID %>%
    dplyr::select(config$table$hierarchyKeys()[1]) %>%
    inner_join(resNACondition)

  # Complete cases
  data_NA_QVal_condition <- completeCases( data_NA_QVal_condition , config)
  stat_peptidFitler <- hierarchyCounts(data_NA_QVal_condition, config)
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
#' config <- spectronautDIAData250_config$clone(deep=TRUE)
#' data <- spectronautDIAData250_analysis
#' tmp <-workflow_Q_NA_filtered_Hierarchy(data, config, hierarchy_level=2)
#' nrow(tmp$data)
#' config <- spectronautDIAData250_config$clone(deep=TRUE)
#' res <-workflow_Q_NA_filtered_Hierarchy(data, config, hierarchy_level=1)
#' nrow(res$data)
#' res$newconfig
#' hierarchyCounts(res$data, res$newconfig)
workflow_Q_NA_filtered_Hierarchy <- function(data,
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
  figs3 <- applyToHierarchyBySample(resDataLog, config, medpolishPly, hierarchy_level = hierarchy_level)
  protIntensity <- figs3 %>% dplyr::select(config$table$hierarchyKeys()[1:hierarchy_level], medpolishPly) %>% unnest()

  newconfig <- make_reduced_hierarchy_config(config,
                                             workIntensity = "medpolish",
                                             hierarchy = config$table$hierarchy[1:hierarchy_level])
  return(list(data = protIntensity, newconfig = newconfig))
}

