# Helper functions -----

#' Keep only those proteins with 2 IDENTIFIED peptides
#'
#' @export
#' @keywords internal
#' @examples
#'
#' library(LFQService)
#' library(tidyverse)
#'
#' resDataStart <- LFQServiceData::testData2954$resDataStart
#' config <-  LFQServiceData::testDataStart2954$config
#'
#' undebug(filter_proteins_by_peptide_count)
#' res <- filter_proteins_by_peptide_count(resDataStart, config)
#'
#' resDataStart <- LFQServiceData::skylineSRM_HL_data
#' config <-  LFQServiceData::skylineconfig_HL
#'
#' resDataStart <- setup_analysis(resDataStart , config)
#' res <- filter_proteins_by_peptide_count(resDataStart, config)
#'
filter_proteins_by_peptide_count <- function(resDataStart,
                                             config){

  # remove single hit wonders
  resDataStart <- nr_B_in_A(resDataStart,config)
  filteredPep <- dplyr::filter(resDataStart$data, !!sym(resDataStart$name) >= config$parameter$min_peptides_protein )

  return(list(data = filteredPep, config = config, name = resDataStart$data))
}

#' normalize data by log2 and robust scaling
#'
#' @export
#' @keywords internal
#' @examples
#' resDataStart <- LFQServiceData::testDataStart2954$resDataStart
#' config <-  LFQServiceData::testDataStart2954$config
#' filterPep <- LFQService:::filter_proteins_by_peptide_count( resDataStart ,  config )
#' normalize_log2_robscale(filterPep$data, filterPep$config)
#'
normalize_log2_robscale <- function(filteredPep, config){
  pepConfig <- config$clone(deep = TRUE)
  pepIntensityNormalized <- transform_work_intensity(filteredPep, pepConfig, log2)
  pepIntensityNormalized <- applyToIntensityMatrix(pepIntensityNormalized,
                                                   pepConfig,
                                                   .func = robust_scale)

  pepIntensityNormalized <- pepIntensityNormalized %>%
    dplyr::rename(transformedIntensity = pepConfig$table$getWorkIntensity())
  pepConfig$table$popWorkIntensity()
  pepConfig$table$setWorkIntensity("transformedIntensity")
  return(list(data = pepIntensityNormalized, config = pepConfig))
}


#' get the difference of two dataset
#'
#' @export
#' @keywords internal
#'
#' @examples
#'
#' resDataStart <- LFQServiceData::testData2954$resDataStart
#' config <-  LFQServiceData::testDataStart2954$config
#'
#' res <- filter_proteins_by_peptide_count(resDataStart, config)
#'
#' resDiff  <- difference(resDataStart, res$data, config)
#'
difference <- function(x, y, config){
  anti_join(x, y, by = config$table$idVars())
}
