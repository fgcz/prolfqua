#' run the modelling using lmer and lm models
#'
#' @export
#'
application_run_modelling_V2 <- function(outpath,
                                         data,
                                         pepConfig,
                                         modelFunction,
                                         contrasts,
                                         modelling_dir="modelling_results_protein" ){
  assign("lfq_write_format", c("xlsx","html"), envir = .GlobalEnv)

  # create result structure
  modelling_path <- file.path(outpath, modelling_dir)
  if(!dir.exists(outpath)){
    dir.create(outpath)
  }
  if(!dir.exists(modelling_path)){
    dir.create(modelling_path)
  }


  #################################################
  ### Do missing value imputation

  res_contrasts_imputed <- workflow_missigness_impute_contrasts(data,
                                                                pepConfig,
                                                                contrasts)
  xx_imputed <- res_contrasts_imputed("long",what = "contrasts")

  ### make contrasts -----

  modellingResult_fun <- workflow_model_analyse(data,
                                                modelFunction,
                                                subject_Id = pepConfig$table$hkeysLevel())

  modellingResult <- modellingResult_fun()
  modellingResult_fun(modelling_path)

  #names(modellingResult)
  modelProteinF <- modellingResult$modellingResult$modelProtein
  res_contrasts <- workflow_contrasts_linfct_V2(modelProteinF,
                                                contrasts,
                                                pepConfig,
                                                prefix =  "contrasts",
                                                contrastfun = modelFunction$contrast_fun)


  # return(list(res_contrasts = res_contrasts, modellingResult_fun = modellingResult_fun))

  xx <- res_contrasts(modelling_path, columns = modelFunction$report_columns)


  merge_contrasts_results <- function(contrast_minimal,
                                      contrasts_imputed,
                                      subject_Id,
                                      modelFunction){
    res <- right_join(contrast_minimal, contrasts_imputed, by=c(subject_Id,"lhs" = "contrast"))
    res <- res %>% rename(contrast = lhs)
    res <- res %>% mutate(pseudo_estimate = case_when(is.na(estimate) ~ imputed, TRUE ~ estimate))
    res <- res %>% mutate(is_pseudo_estimate = case_when(is.na(estimate) ~ TRUE, TRUE ~ FALSE))

    for(column in modelFunction$report_columns){
      res <- res %>% mutate(!!sym(paste0("pseudo_",column)) := case_when(is.na(estimate) ~ 0, TRUE ~ !!sym(column)))
    }
    res <- res %>% dplyr::select(-imputed, -meanArea)
    return(res)
  }

  contrast_results <- merge_contrasts_results(xx$contrast_minimal,
                                              xx_imputed,
                                              subject_Id = pepConfig$table$hkeysLevel(),
                                              modelFunction = modelFunction)
  separate_hierarchy(contrast_results, config) -> filtered_dd

  lfq_write_table(filtered_dd, path = file.path(modelling_path, "foldchange_estimates.csv"))
  return(list(res_contrasts  = res_contrasts, modellingResult_fun = modellingResult_fun))
}


#' run the modelling using lmer and lm models - DEPRECATED use version V2
#'
#' @export
#'
application_run_modelling <- function(outpath,
                                      data,
                                      pepConfig,
                                      modelFunction,
                                      contrasts,
                                      modelling_dir="modelling_results_protein",
                                      DEBUG = FALSE){
  warning("DEPRECATED application_run_modelling! use application_run_modelling_V2 instead!")
  assign("lfq_write_format", c("xlsx","html"), envir = .GlobalEnv)

  # create result structure
  modelling_path <- file.path(outpath, modelling_dir)
  if(!dir.exists(outpath)){
    dir.create(outpath)
  }
  if(!dir.exists(modelling_path)){
    dir.create(modelling_path)
  }


  #################################################
  ### Do missing value imputation

  res_contrasts_imputed <- workflow_missigness_impute_contrasts(data,
                                                                pepConfig,
                                                                contrasts)

  ### make contrasts -----

  modellingResult_fun <- workflow_model_analyse(data,
                                                modelFunction,
                                                subject_Id = pepConfig$table$hkeysLevel())

  modellingResult <- modellingResult_fun()
  modellingResult_fun(modelling_path)

  #names(modellingResult)

  m <- get_complete_model_fit(modellingResult$modellingResult$modelProtein)

  #factor_contrasts <- linfct_factors_contrasts(m$linear_model[[1]])
  linfct <- linfct_from_model(m$linear_model[[1]], as_list = FALSE)
  linfct_A <- linfct_matrix_contrasts(linfct, contrasts)

  if(DEBUG){
    pdf(file.path(modelling_path,"Linear_functions.pdf"), width=18, height=10)
    quantable::imageWithLabels(t(linfct_A),
                               col = quantable::getBlueWhiteRed(),
                               marLeft = c(8,10,4.1,2.1))
    dev.off()
  }

  modelProteinF <- modellingResult$modellingResult$modelProtein %>%
    dplyr::filter(exists_lmer == TRUE)

  res_contrasts <- workflow_contrasts_linfct(modelProteinF,
                                             linfct_A,
                                             pepConfig,
                                             prefix =  "contrasts",
                                             contrastfun = modelFunction$contrast_fun)

  xx <- res_contrasts(modelling_path, columns = modelFunction$report_columns)
  xx_imputed <- res_contrasts_imputed("long",what = "contrasts")

  merge_contrasts_results <- function(contrast_minimal,
                                      xx_imputed,
                                      subject_Id,
                                      modelFunction){
    res <- right_join(contrast_minimal, xx_imputed, by=c(subject_Id,"lhs" = "contrast"))
    res <- res %>% rename(contrast = lhs)
    res <- res %>% mutate(pseudo_estimate = case_when(is.na(estimate) ~ imputed, TRUE ~ estimate))
    res <- res %>% mutate(is_pseudo_estimate = case_when(is.na(estimate) ~ TRUE, TRUE ~ FALSE))

    for(column in modelFunction$report_columns){
      res <- res %>% mutate(!!sym(paste0("pseudo_",column)) := case_when(is.na(estimate) ~ 0, TRUE ~ !!sym(column)))
    }
    res <- res %>% dplyr::select(-imputed, -meanArea)
    return(res)
  }

  contrast_results <- merge_contrasts_results(xx$contrast_minimal, xx_imputed,
                                              subject_Id = pepConfig$table$hkeysLevel(), modelFunction = modelFunction)
  separate_hierarchy(contrast_results, config) -> filtered_dd

  lfq_write_table(filtered_dd, path = file.path(modelling_path, "foldchange_estimates.csv"))
}

#' add Annotation to an MQ output
#' @export
#'
application_set_up_MQ_run <- function(outpath,
                                      inputMQfile,
                                      inputAnntation,
                                      config,
                                      id_extractor = function(df){fgczgseaora::get_UniprotID_from_fasta_header(df, idcolumn = "top_protein")},
                                      qcdir = "qc_results"){
  assign("lfq_write_format", c("xlsx"), envir = .GlobalEnv)
  config$table$hierarchy[["protein_Id"]] <- c("top_protein","protein.group.id", "UniprotID")
  # create result structure
  qc_path <- file.path(outpath, qcdir )


  if(!dir.exists(outpath)){
    dir.create(outpath)
  }

  message(qc_path,"\n")
  if(!dir.exists(qc_path)){
    dir.create(qc_path)
  }


  ## read the data

  resPepProtAnnot <- tidyMQ_modificationSpecificPeptides(inputMQfile)

  {
    height <- length(unique(resPepProtAnnot$raw.file))/2 * 300
    png(file.path(qc_path, "retention_time_plot.png"), height = height, width=1200)
    resPepProtVis <- resPepProtAnnot %>% dplyr::filter(mod.peptide.intensity > 4)
    tmp <- ggplot(resPepProtVis, aes(x = retention.time, y= log2(mod.peptide.intensity))) + geom_point(alpha=1/20, size=0.3) + facet_wrap(~raw.file, ncol=2 )
    print(tmp)
    dev.off()
  }

  peptidestxt <- tidyMQ_Peptides(inputMQfile)
  resPepProtAnnot <- tidyMQ_from_modSpecific_to_peptide(resPepProtAnnot, peptidestxt)
  resPepProtAnnot <- tidyMQ_top_protein_name(resPepProtAnnot)



  # add annotation
  annotation <- readxl::read_xlsx(inputAnntation)
  noData <- annotation[!annotation$raw.file %in% resPepProtAnnot$raw.file,]


  if(nrow(noData)){
    message("some files in annotation have no measurements")
    message(paste(noData,collapse = " "))
  }

  measSamples <- unique(resPepProtAnnot$raw.file)
  noAnnot <- measSamples[! measSamples%in% annotation$raw.file]

  if(length(noAnnot) > 0 ){
    message("some measured samples have no annotation!")
    message(paste(noAnnot,collapse = " "))
  }

  resPepProtAnnot <- inner_join(annotation, resPepProtAnnot, by= "raw.file")
  ####

  ###  Setup analysis ####


  resPepProtAnnot <- resPepProtAnnot %>% dplyr::filter(reverse !=  TRUE)

  resPepProtAnnot$isotope = "light"
  resPepProtAnnot <- id_extractor(resPepProtAnnot)
  resDataStart <- setup_analysis(resPepProtAnnot, config)

  resDataStart <- remove_small_intensities( resDataStart, config, threshold = 4 ) %>%
    completeCases(config)
  return(list(data = resDataStart,config=config, qc_path = qc_path))
}




#' preprocess peptide data, compute protein data, store results in qc_path folder
#' @export
#'
application_summarize_data_pep_to_prot <- function(data,
                                                   config,
                                                   qc_path,
                                                   DEBUG= FALSE,
                                                   WRITE_PROTS=TRUE){
  assign("lfq_write_format", c("xlsx"), envir = .GlobalEnv)

  results <- LFQService::workflow_MQ_protoV1(
    data,
    config,
    outpath,
    peptideFilterFunction = LFQService:::.workflow_MQ_filter_peptides_V3 )

  wideFRAME <- LFQService::toWideConfig(results$pepIntensityNormalized,
                                        results$config_pepIntensityNormalized)


  lfq_write_table(separate_hierarchy(wideFRAME$data,
                                     results$config_pepIntensityNormalized),
                  file.path(qc_path, "peptide_intensities.csv"))


  if(!DEBUG){
    LFQService::render_MQSummary_rmd(data,
                                     config$clone(deep=TRUE),
                                     pep=TRUE,
                                     workdir = ".",
                                     dest_path = qc_path,
                                     dest_file_name="peptide_intensities_qc",
                                     format = "html")
  }


  ### PROTEIN QUANTIFICATION ####
  protintensity <- medpolish_protein_quants( results$pepIntensityNormalized,
                                             results$config_pepIntensityNormalized )

  unnest <- protintensity("unnest")
  quants_write(unnest$data, unnest$config, qc_path, na_fraction = 0.5)


  figs <- protintensity("plot")

  if(WRITE_PROTS){
    pdf(file.path(qc_path, "protein_intensities_inference_figures.pdf"))
    lapply(figs$plot, print)
    dev.off()
  }


  # render protein quantification reprot
  if(!DEBUG){
    LFQService::render_MQSummary_rmd(protintensity("unnest")$data,
                                     protintensity("unnest")$config$clone(deep=TRUE),
                                     pep=FALSE,
                                     workdir = ".",
                                     dest_path = qc_path,
                                     dest_file_name="protein_intensities_qc",
                                     format="html"
    )
  }
  return(list(results  = results, prot_results = protintensity))
}


#' DEPRECATED
#' @export
application_summarize_data <- function(data,
                                       config,
                                       qc_path,
                                       DEBUG= FALSE,
                                       WRITE_PROTS=TRUE){
  message("application_summarize_data is deprecated\n")
  stop("use application_summarize_data_pep_to_prot\n")
}



#' preprocess peptide data, compute protein data, store results in qc_path folder
#' @export
#'
application_summarize_compound <- function(data,
                                           config,
                                           qc_path,
                                           DEBUG= FALSE,
                                           WRITE_PROTS=TRUE,
                                           prefix = c("peptide", "compound")){
  prefix <- match.arg(prefix)
  assign("lfq_write_format", c("xlsx"), envir = .GlobalEnv)


  results <- LFQService:::.workflow_MQ_normalize_log2_robscale(data, config)

  wideFRAME <- LFQService::toWideConfig(results$data,
                                        results$config)


  #return(wideFRAME)

  lfq_write_table(separate_hierarchy(wideFRAME$data,
                                     results$config),
                  path = file.path(qc_path, paste0(prefix, "_intensities.csv")))


  if(!DEBUG){
    LFQService::render_MQSummary_rmd(results$data,
                                     results$config$clone(deep=TRUE),
                                     pep=TRUE,
                                     workdir = ".",
                                     dest_path = qc_path,
                                     dest_file_name = paste0(prefix, "_intensities_qc"),
                                     format = "html")
  }


  quants_write(results$data, results$config, qc_path, prefix = prefix)



  if(WRITE_PROTS){
    figs <- plot_hierarchies_line_df(results$data, results$config )
    pdf(file.path(qc_path, paste0(prefix, "_intensities_.pdf")))
    lapply(figs, print)
    dev.off()
  }


  return(list(data  = results$data, config = results$config))
}
#################################################
### Do missing value imputation


