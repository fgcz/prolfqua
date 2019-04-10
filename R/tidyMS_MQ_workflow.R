
.makeFigs <- function(filteredPep, config){
  factor_level <- config$table$factorLevel
  proteinIDsymbol <- sym(config$table$hierarchyKeys()[1])
  xnested <- filteredPep %>% group_by(UQ(proteinIDsymbol)) %>% nest()
  figs <- xnested %>% mutate(plot = map2(data, UQ(proteinIDsymbol) , plot_hierarchies_line, factor_level = factor_level, config ))
  figs <- figs %>% mutate(plotboxplot = map2(data, UQ(proteinIDsymbol) , plot_hierarchies_boxplot, config , factor_level = factor_level))
  return(figs)
}


.workflow_MQ_filter_peptides <- function(resDataStart, config, percent = 50){
  config <- config$clone(deep = TRUE)
  resNACondition <- filter_factor_levels_by_missing(resDataStart,
                                                    config,
                                                    percent = percent,
                                                    factor_level = config$table$factorLevel
  )

  resNACondition <- completeCases(resNACondition, config)
  filteredPep <- summarizeHierarchy(resNACondition, config)

  filteredPep <- inner_join(filteredPep, resNACondition, by="protein_Id", suffix = c(".NA_filt", ""))
  filteredPep <- filteredPep %>% dplyr::filter( peptide_Id_n.NA_filt >= config$parameter$min_peptides_protein)
  return(list(data=filteredPep, config=config))
}



.workflow_MQ_filter_peptides_V2 <- function(resDataStart, config, percent = 50){
  config <- config$clone(deep = TRUE)
  resNACondition <- filter_factor_levels_by_missing(resDataStart,
                                                    config,
                                                    percent = percent,
                                                    factor_level = config$table$factorLevel
  )

  resNACondition <- resNACondition %>% dplyr::select(protein_Id) %>% distinct() %>% inner_join(resDataStart , by="protein_Id")
  resNACondition <- completeCases(resNACondition, config)
  filteredPep <- summarizeHierarchy(resNACondition, config)

  filteredPep <- inner_join(filteredPep, resNACondition, by="protein_Id", suffix = c(".NA_filt", ""))
  filteredPep <- filteredPep %>% dplyr::filter( peptide_Id_n.NA_filt >= config$parameter$min_peptides_protein)
  return(list(data=filteredPep, config=config))
}



.workflow_MQ_normalize_log2_robscale <- function(filteredPep, config){
  pepConfig <- config$clone(deep = TRUE)
  pepIntensityNormalized <- transform_work_intensity(filteredPep, pepConfig, log2)
  pepIntensityNormalized <- applyToIntensityMatrix(pepIntensityNormalized, pepConfig, .func = robust_scale)
  pepIntensityNormalized <- pepIntensityNormalized %>% dplyr::rename(transformedIntensity = pepConfig$table$getWorkIntensity())
  pepConfig$table$popWorkIntensity()
  pepConfig$table$setWorkIntensity("transformedIntensity")
  return(list(data = pepIntensityNormalized, config = pepConfig))
}

#' median polish from normalized peptide intensities
#' @export
#' @examples
#' resultsV12954 <- LFQService::resultsV12954
#' res <- workflow_MQ_protein_quants(resultsV12954)
#' dim(res$protintensity)
#'
workflow_MQ_protein_quants <- function(results){
  pepIntensityNormalized <- results$pepIntensityNormalized
  config <- results$config_pepIntensityNormalized
  configProt <- config$clone(deep = TRUE)
  xnested <- pepIntensityNormalized %>% group_by_at(names(configProt$table$hkeysLevel())) %>% nest()
  protintensity <- LFQService::applyToHierarchyBySample(pepIntensityNormalized,configProt, medpolishPly,unnest = TRUE)
  configProt <- protintensity$newconfig
  protintensity <- protintensity$unnested
  return(list(data = protintensity, config = configProt))
}

#' runs data preprocessing for peptide level data based protein modelling
#' @export
#' @param peptideFilterFunction can be either .workflow_MQ_filter_peptides or .workflow_MQ_filter_peptides_V2
#' @examples
#' #testDataStart2954 <- readRDS("c:/Users/wolski/prog/LFQService/data/testDataStart2954.rds")
#' #usethis::use_data(testDataStart2954)
#'
#' testDataStart2954 <- LFQService::testDataStart2954
#' path <- "dummy_test"
#' resultsV12954 <- LFQService::workflow_MQ_protoV1(testDataStart2954$resDataStart, testDataStart2954$config, path ,
#'                                            peptideFilterFunction = LFQService:::.workflow_MQ_filter_peptides_V2 )
#' #usethis::use_data(resultsV12954)
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
  config <- pepIntensityNormalized$config
  pepIntensityNormalized <- pepIntensityNormalized$data
  RESULTS$config_pepIntensityNormalized <- config
  RESULTS$pepIntensityNormalized <- pepIntensityNormalized

  # Summarize number of peptides with more than 2
  x3_start <- summarizeHierarchy(RESULTS$resDataStart, RESULTS$config_resDataStart)
  x3_start <- x3_start %>% mutate(protein_with = case_when(peptide_Id_n == 1 ~ "one",
                                                           peptide_Id_n > 1 ~ "two and more"))
  RESULTS$nrPeptidesPerProtein_start <- x3_start %>% group_by(protein_with) %>% summarize(n=n())



  # Summarize filtered data - number of peptides with more than 2
  x3_filt <- summarizeHierarchy(RESULTS$filteredPep, RESULTS$config_filteredPep)
  x3_filt <- x3_filt %>% mutate(protein_with = case_when(peptide_Id_n == 1 ~ "one",
                                                         peptide_Id_n > 1 ~ "two and more"))
  RESULTS$nrPeptidesPerProtein_filtered <- x3_filt %>% group_by(protein_with) %>% summarize(n=n())

#  knitr::kable(res, caption = "nr of proteins with more than on peptide.")

  x3_start %>% dplyr::select( -protein_with) %>% dplyr::filter(peptide_Id_n  > 1) -> x3_start
  x3_filt %>% dplyr::select( -protein_with) %>% dplyr::filter(peptide_Id_n  > 1) -> x3_filt

  res <- left_join(x3_start, x3_filt , by="protein_Id", suffix=c(".start",".filt")) %>% arrange(peptide_Id_n.filt)
  RESULTS$removed_proteins <- res %>% dplyr::filter(is.na(peptide_Id_n.filt))

  RESULTS$removed_peptides <- inner_join(RESULTS$removed_proteins, RESULTS$resDataStart)
  ### PLOTTING
  return(RESULTS)
}

#' generates peptide level plots for all Proteins
#' @export
#'
workflow_MQ_figs_protoV1 <- function(RESULTS){
  figs_raw <- .makeFigs(RESULTS$filteredPep, RESULTS$config_filteredPep)
  figs_normalized <- .makeFigs(RESULTS$pepIntensityNormalized, RESULTS$config_pepIntensityNormalized)

  saveRDS(figs_raw, file=file.path(RESULTS$path, "figures_allProteins_RAW.Rda"))
  saveRDS(figs_normalized, file=file.path(RESULTS$path, "figures_allProteins_TRANS.Rda"))
}