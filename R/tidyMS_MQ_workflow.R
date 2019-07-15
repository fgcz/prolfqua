# Helper functions -----

#' Filter peptides for NA's within factors and keep only those proteins with 2 QUANTIFIED peptides
#'
#' @export
#' @examples
#' library(tidyverse)
#' testDataStart2954 <- LFQService::testDataStart2954
#' path <- "dummy_test"
#' dd <- LFQService:::.workflow_MQ_filter_peptides( testDataStart2954$resDataStart ,  testDataStart2954$config )
#' hierarchyCounts(dd$data, dd$config)
.workflow_MQ_filter_peptides <- function(resDataStart, config, percent = 50){
  config <- config$clone(deep = TRUE)
  resNACondition <- filter_factor_levels_by_missing(resDataStart,
                                                    config,
                                                    percent = percent
  )

  proteinsWith2Peptides <- summarizeHierarchy(resNACondition, config)
  proteinsWith2Peptides <- proteinsWith2Peptides %>% dplyr::filter( peptide_Id_n >= config$parameter$min_peptides_protein)
  filteredPep <- dplyr::inner_join(proteinsWith2Peptides, resNACondition)
  return(list(data=filteredPep, config=config))
}


#' Filter peptides for NA's within factors and keep only those proteins with 2 IDENTIFIED peptides
#' @export
#' @examples
#' library(tidyverse)
#' testDataStart2954 <- LFQService::testDataStart2954
#' dd <- LFQService:::.workflow_MQ_filter_peptides_V2( testDataStart2954$resDataStart ,  testDataStart2954$config )
#' hierarchyCounts(dd$data, dd$config)
.workflow_MQ_filter_peptides_V2 <-  function(resDataStart,
                                             config,
                                             percent = 50,
                                             nr_peptide_id = 2){


  config <- config$clone(deep = TRUE)

  # get proteins with more than 1 peptide before NA filtering.
  summaryH <- summarizeHierarchy(resDataStart, config)
  proteinsWith2Peptides <- summaryH %>% dplyr::filter( peptide_Id_n >= 2)
  # fitler for missingness
  resNACondition <- filter_factor_levels_by_missing(resDataStart,
                                                    config,
                                                    percent = percent)
  filteredPep <- dplyr::inner_join(proteinsWith2Peptides, resNACondition, by="protein_Id")
  return(list(data=filteredPep, config=config))
}

#' Keep only those proteins with 2 IDENTIFIED peptides
#' @export
.workflow_MQ_filter_peptides_V3 <-  function(resDataStart,
                                             config,
                                             percent = 50,
                                             nr_peptide_id = 2){


  config <- config$clone(deep = TRUE)

  # get proteins with more than 1 peptide before NA filtering.
  summaryH <- summarizeHierarchy(resDataStart, config)
  proteinsWith2Peptides <- summaryH %>% dplyr::filter( peptide_Id_n >= 2)


  filteredPep <- dplyr::inner_join(proteinsWith2Peptides, resDataStart, by="protein_Id")
  return(list(data=filteredPep, config=config))
}

#'
#'
#' @examples
#' resDataStart <- LFQService::testDataStart2954$resDataStart
#' config <-  LFQService::testDataStart2954$config
#' filterPep <- LFQService:::.workflow_MQ_filter_peptides_V2( resDataStart ,  config )
#' .workflow_MQ_normalize_log2_robscale(filterPep$data, filterPep$config)
.workflow_MQ_normalize_log2_robscale <- function(filteredPep, config){

  pepConfig <- config$clone(deep = TRUE)
  pepIntensityNormalized <- transform_work_intensity(filteredPep, pepConfig, log2)
  pepIntensityNormalized <- applyToIntensityMatrix(pepIntensityNormalized, pepConfig, .func = robust_scale)
  pepIntensityNormalized <- pepIntensityNormalized %>%
    dplyr::rename(transformedIntensity = pepConfig$table$getWorkIntensity())
  pepConfig$table$popWorkIntensity()
  pepConfig$table$setWorkIntensity("transformedIntensity")
  return(list(data = pepIntensityNormalized, config = pepConfig))
}

# Workflow function ----

#' runs data preprocessing for peptide level data based protein modelling
#' @export
#' @param peptideFilterFunction can be either .workflow_MQ_filter_peptides or .workflow_MQ_filter_peptides_V2
#' @examples
#' resDataStart <- LFQService::testData2954$resDataStart
#' config <-  LFQService::testDataStart2954$config
#' summarizeHierarchy(resDataStart, config)
#' path <- "dummy_test"
#' resultsV12954 <- LFQService::workflow_MQ_protoV1(resDataStart,
#'  config,
#'  path ,
#'  peptideFilterFunction = LFQService:::.workflow_MQ_filter_peptides_V3 )
#'
#'
#' LFQService:::.workflow_MQ_filter_peptides_V2( resDataStart ,  config )
#'
workflow_MQ_protoV1 <- function( resDataStart,
                                 config,
                                 path,
                                 peptideFilterFunction = LFQService:::.workflow_MQ_filter_peptides ){
  RESULTS <- list()
  RESULTS$path <- path
  config <- config$clone(deep=TRUE)
  RESULTS$config_resDataStart <- config
  RESULTS$resDataStart <- resDataStart

  filteredPep <- peptideFilterFunction( resDataStart , config )
  config <- filteredPep$config
  filteredPep <- filteredPep$data

  # do Normalization

  RESULTS$config_filteredPep <- config
  RESULTS$filteredPep <- filteredPep

  pepIntensityNormalized <- .workflow_MQ_normalize_log2_robscale(filteredPep, config)

  RESULTS$config_pepIntensityNormalized <- pepIntensityNormalized$config
  RESULTS$pepIntensityNormalized <- pepIntensityNormalized$data

  # Summarize number of peptides with more than 2
  x3_start <- summarizeHierarchy(RESULTS$resDataStart, RESULTS$config_resDataStart)
  x3_start <- x3_start %>% dplyr::mutate(protein_with = dplyr::case_when(peptide_Id_n == 1 ~ "one",
                                                                         peptide_Id_n > 1 ~ "two and more"))
  RESULTS$nrPeptidesPerProtein_start <- x3_start %>% dplyr::group_by(protein_with) %>%
    dplyr::summarize(n=n())



  # Summarize filtered data - number of peptides with more than 2
  x3_filt <- summarizeHierarchy(RESULTS$filteredPep, RESULTS$config_filteredPep)
  x3_filt <- x3_filt %>% dplyr::mutate(protein_with = dplyr::case_when(peptide_Id_n == 1 ~ "one",
                                                                       peptide_Id_n > 1 ~ "two and more"))
  RESULTS$nrPeptidesPerProtein_filtered <- x3_filt %>% dplyr::group_by(protein_with) %>%
    dplyr::summarize(n=n())

  #  knitr::kable(res, caption = "nr of proteins with more than on peptide.")

  x3_start %>% dplyr::select( -protein_with) %>% dplyr::filter(peptide_Id_n  > 1) -> x3_start
  x3_filt %>% dplyr::select( -protein_with) %>% dplyr::filter(peptide_Id_n  > 1) -> x3_filt

  res <- left_join(x3_start, x3_filt , by="protein_Id", suffix=c(".start",".filt")) %>% arrange(peptide_Id_n.filt)
  RESULTS$removed_proteins <- res %>% dplyr::filter(is.na(peptide_Id_n.filt))

  RESULTS$removed_peptides <- dplyr::inner_join(RESULTS$removed_proteins, RESULTS$resDataStart)
  ### PLOTTING
  return(RESULTS)
}

