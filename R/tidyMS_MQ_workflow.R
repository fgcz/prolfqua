
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
  filteredPep <- filteredPep %>% filter( peptide_Id_n.NA_filt >= config$parameter$min_peptides_protein)
  return(list(data=filteredPep, config=config))
}


.workflow_MQ_normalize_log2_robscale <- function(filteredPep, config){
  pepConfig <- config$clone(deep = TRUE)
  pepIntensityNormalized <- transform_work_intensity(filteredPep, pepConfig, log2)
  pepIntensityNormalized <- applyToIntensityMatrix(pepIntensityNormalized, pepConfig, .func = robust_scale)
  pepIntensityNormalized <- pepIntensityNormalized %>% rename(transformedIntensity = pepConfig$table$getWorkIntensity())
  pepConfig$table$popWorkIntensity()
  pepConfig$table$setWorkIntensity("transformedIntensity")
  return(list(data = pepIntensityNormalized, config = pepConfig))
}

#' median polish from normalized peptide intensities
#' @export
workflow_MQ_protein_quants <- function(pepIntensityNormalized , results){
  pepIntensityNormalized <- results$pepIntensityNormalized
  pepIntensityNormalized <- results$config_pepIntensityNormalized
  configProt <- config$clone(deep = TRUE)
  xnested <- pepIntensityNormalized %>% group_by_at(names(configProt$table$hkeysLevel())) %>% nest()
  protintensity <- LFQService::applyToHierarchyBySample(pepIntensityNormalized,configProt, medpolishPly,unnest = TRUE)
  configProt <- protintensity$newconfig
  protintensity <- protintensity$unnested
  return(list(data = protintensity, config = configProt))
}

#' runs data preprocessing for peptide level data based protein modelling
#' @export
#'
worklfow_MQ_protoV1 <- function(resDataStart, config, path){
  RESULTS <- list()
  RESULTS$path <- path
  config <- config$clone(deep=TRUE)
  RESULTS$config_resDataStart <- config
  RESULTS$resDataStart <- resDataStart

  filteredPep <- .workflow_MQ_filter_peptides( resDataStart , config )
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
