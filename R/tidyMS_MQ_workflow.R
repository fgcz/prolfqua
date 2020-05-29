# Helper functions -----

#' Keep only those proteins with 2 IDENTIFIED peptides
#' @export
#'
filter_proteins_by_peptide_count <-  function(resDataStart,
                                             config){
  # filter_proteins_by_peptide_count renamed from .workflow_MQ_filter_peptides_V3
  config <- config$clone(deep = TRUE)
  nr_peptide_id <- config$parameter$min_peptides_protein
  # get proteins with more than 1 peptide before NA filtering.
  summaryH <- summarize_hierarchy(resDataStart, config)
  proteinsWith2Peptides <- summaryH %>% dplyr::filter( peptide_Id_n >= nr_peptide_id)
  filteredPep <- dplyr::inner_join(proteinsWith2Peptides, resDataStart, by = "protein_Id")
  return(list(data = filteredPep, config = config))
}

#' normlize data by log2 and robust scaling
#' @export
#' @examples
#' resDataStart <- LFQService::testDataStart2954$resDataStart
#' config <-  LFQService::testDataStart2954$config
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


#' Compare a dataset before and after filtering.
#'
#' @export
#' @examples
#' resDataStart <- LFQService::testData2954$resDataStart
#' config <-  LFQService::testDataStart2954$config
#' res <- filter_proteins_by_peptide_count(resDataStart, config)
#' summarize_filtering(list(data = resDataStart, config=config), res)
#' summarize_filtering(list(data = resDataStart, config=config), res, rm_one_hit=FALSE)
summarize_filtering <- function(startData,
                                endData,
                                rm_one_hit = TRUE){
  SUMMARY <- list()
  x3_start <- summarize_hierarchy(startData$data, startData$config)
  x3_start <- x3_start %>%
    dplyr::mutate(protein_with = dplyr::case_when(peptide_Id_n == 1 ~ "one",
                                                  peptide_Id_n > 1 ~ "two and more"))

  # Summarize filtered data - number of peptides with more than 2
  x3_filt <- summarize_hierarchy(endData$data, endData$config)
  x3_filt <- x3_filt %>%
    dplyr::mutate(protein_with = dplyr::case_when(peptide_Id_n == 1 ~ "one",
                                                  peptide_Id_n > 1 ~ "two and more"))

  SUMMARY$nrPeptidesPerProtein_start <- x3_start %>%
    dplyr::group_by( protein_with ) %>%
    dplyr::summarize(n = n() )
  SUMMARY$nrPeptidesPerProtein_filtered <- x3_filt %>% dplyr::group_by(protein_with) %>%
    dplyr::summarize(n = n())

  #  knitr::kable(res, caption = "nr of proteins with more than on peptide.")

  if (rm_one_hit) {
    x3_start %>% dplyr::select( -protein_with) %>% dplyr::filter(peptide_Id_n  > 1) -> x3_start
    x3_filt %>% dplyr::select( -protein_with) %>% dplyr::filter(peptide_Id_n  > 1) -> x3_filt
  }
  res <- dplyr::left_join(x3_start, x3_filt , by = "protein_Id", suffix = c(".start",".filt")) %>% arrange(peptide_Id_n.filt)
  SUMMARY$removed_proteins <- res %>% dplyr::filter(is.na(peptide_Id_n.filt))
  SUMMARY$removed_peptides <- dplyr::inner_join(SUMMARY$removed_proteins, startData$data, by = "protein_Id")
  return(SUMMARY)
}

# Workflow function ----

#' runs peptide data preprocessing for peptide level protein modelling
#' first runs peptide filtering
#' then it runs data noramlization (log2 transform intensities and apply robust z transformation)
#' But mainly also summarizes data filtering (reports which proteins were removed etc.)
#' @keywords internal
#' @param peptideFilterFunction can be filter_proteins_by_peptide_count
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' resDataStart <- LFQService::testData2954$resDataStart
#' config <-  LFQService::testDataStart2954$config
#' summarize_hierarchy(resDataStart, config)
#' res <- filter_proteins_by_peptide_count(resDataStart, config)
#' res <- normalize_log2_robscale(res$data, res$config)
#' summarize_filtering(list(data = resDataStart, config=config), res)
#' # resDataStart <- LFQService::skylineSRM_HL_data
#' # config <-  LFQService::skylineconfig_HL
#'
#' # res <- filter_proteins_by_peptide_count(resDataStart, config)

#'
#'
#'
#'
#'
workflow_MQ_protoV1 <- function(){
  error("Replace by call to:
        filter_proteins_by_peptide_count and
        normalize_log2_robscale
        ")
}

