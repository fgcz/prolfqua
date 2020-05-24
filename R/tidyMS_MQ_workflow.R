# Helper functions -----

#' Keep only those proteins with 2 IDENTIFIED peptides
#' @export
#'
.filter_proteins_by_peptide_count <-  function(resDataStart,
                                             config){
  # .filter_proteins_by_peptide_count renamed from .workflow_MQ_filter_peptides_V3
  config <- config$clone(deep = TRUE)
  nr_peptide_id <- config$parameter$min_peptides_protein
  # get proteins with more than 1 peptide before NA filtering.
  summaryH <- summarize_hierarchy(resDataStart, config)
  proteinsWith2Peptides <- summaryH %>% dplyr::filter( peptide_Id_n >= nr_peptide_id)
  filteredPep <- dplyr::inner_join(proteinsWith2Peptides, resDataStart, by = "protein_Id")
  return(list(data = filteredPep, config = config))
}

#'
#'
#' @examples
#' resDataStart <- LFQService::testDataStart2954$resDataStart
#' config <-  LFQService::testDataStart2954$config
#' filterPep <- LFQService:::.filter_proteins_by_peptide_count( resDataStart ,  config )
#' .normalize_log2_robscale(filterPep$data, filterPep$config)
#'
.normalize_log2_robscale <- function(filteredPep, config){
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

# Workflow function ----

#' runs peptide data preprocessing for peptide level protein modelling
#' first runs peptide filtering
#' then it runs data noramlization (log2 transform intensities and apply robust z transformation)
#' But mainly also summarizes data filtering (reports which proteins were removed etc.)
#'
#' @export
#' @param peptideFilterFunction can be .filter_proteins_by_peptide_count
#' @examples
#' resDataStart <- LFQService::testData2954$resDataStart
#' config <-  LFQService::testDataStart2954$config
#' summarize_hierarchy(resDataStart, config)
#' path <- "dummy_test"
#' resultsV12954 <- LFQService::workflow_MQ_protoV1(resDataStart,
#'  config,
#'  path ,
#'  peptideFilterFunction = LFQService:::.filter_proteins_by_peptide_count )
#'
#'
#' LFQService:::.filter_proteins_by_peptide_count( resDataStart ,  config )
#'
workflow_MQ_protoV1 <- function(resDataStart,
                                 config,
                                 path,
                                 peptideFilterFunction = LFQService:::.workflow_MQ_filter_peptides ){
  RESULTS <- list()
  RESULTS$path <- path
  config <- config$clone(deep = TRUE)
  RESULTS$config_resDataStart <- config
  RESULTS$resDataStart <- resDataStart

  filteredPep <- peptideFilterFunction( resDataStart , config )
  config <- filteredPep$config
  filteredPep <- filteredPep$data

  # do Normalization

  RESULTS$config_filteredPep <- config
  RESULTS$filteredPep <- filteredPep

  pepIntensityNormalized <- .normalize_log2_robscale(filteredPep,
                                                                 config)

  RESULTS$config_pepIntensityNormalized <- pepIntensityNormalized$config
  RESULTS$pepIntensityNormalized <- pepIntensityNormalized$data

  # Summarize number of peptides with more than 2
  x3_start <- summarize_hierarchy(RESULTS$resDataStart, RESULTS$config_resDataStart)
  x3_start <- x3_start %>%
    dplyr::mutate(protein_with = dplyr::case_when(peptide_Id_n == 1 ~ "one",
                                                  peptide_Id_n > 1 ~ "two and more"))

  RESULTS$nrPeptidesPerProtein_start <- x3_start %>% dplyr::group_by(protein_with) %>%
    dplyr::summarize(n=n())

  # Summarize filtered data - number of peptides with more than 2
  x3_filt <- summarize_hierarchy(RESULTS$filteredPep, RESULTS$config_filteredPep)
  x3_filt <- x3_filt %>%
    dplyr::mutate(protein_with = dplyr::case_when(peptide_Id_n == 1 ~ "one",
                                                  peptide_Id_n > 1 ~ "two and more"))
  RESULTS$nrPeptidesPerProtein_filtered <- x3_filt %>% dplyr::group_by(protein_with) %>%
    dplyr::summarize(n=n())

  #  knitr::kable(res, caption = "nr of proteins with more than on peptide.")

  x3_start %>% dplyr::select( -protein_with) %>% dplyr::filter(peptide_Id_n  > 1) -> x3_start
  x3_filt %>% dplyr::select( -protein_with) %>% dplyr::filter(peptide_Id_n  > 1) -> x3_filt

  res <- left_join(x3_start, x3_filt , by="protein_Id", suffix = c(".start",".filt")) %>% arrange(peptide_Id_n.filt)
  RESULTS$removed_proteins <- res %>% dplyr::filter(is.na(peptide_Id_n.filt))
  RESULTS$removed_peptides <- dplyr::inner_join(RESULTS$removed_proteins, RESULTS$resDataStart)
  ### PLOTTING
  return(RESULTS)
}

