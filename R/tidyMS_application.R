.columnsImputed <- function(all, contrasts){
  get_sides <- function(contrast){
    bb <- str_split(contrast,"[-+]")
    bb <- gsub("[ `]", "", bb[[1]])
    return(bb)
  }

  res <- NULL

  for(i in 1:length(contrasts)){
    cname <- names(contrasts)[i]
    cc <- get_sides(contrasts[i])

    tmp <- all %>% dplyr::filter(contrast %in% c(cname,cc))
    tmp <- tmp %>% dplyr::select(-meanArea) %>% tidyr::spread(contrast , imputed)

    tmp <- tmp %>% add_column(lhs = cname,.after = 1)
    tmp <- tmp %>% add_column(c1_name = cc[1],.after = 2)
    tmp <- tmp %>% add_column(c2_name = cc[2],.after = 3)
    tmp <- tmp %>% dplyr::rename(c1 = !!sym(cc[1]), c2 = !!sym(cc[2]), estimate = !!sym(cname))
    res <- dplyr::bind_rows(res,tmp)
  }
  return(res)
}

# merges contrasts and imputed contrasts
.makeResult_contrasts <- function(contrasts_xx,
                                  contrasts_xx_imputed,
                                  subject_Id,
                                  remove_imputed = TRUE ){

  contrast_results <- right_join( contrasts_xx$contrast_minimal,
                                  contrasts_xx_imputed,
                                  by= c(subject_Id,
                                        "lhs", "c1_name", "c2_name"), suffix = c("","_imputed"))
  contrast_results <- dplyr::rename(contrast_results, contrast = lhs) #
  contrast_results <- contrast_results %>%
    mutate(pseudo_estimate = case_when(is.na(estimate) ~ estimate_imputed, TRUE ~ estimate))

  contrast_results <- contrast_results %>%
    mutate(is_pseudo_estimate = case_when(is.na(estimate) ~ TRUE, TRUE ~ FALSE))
  if(remove_imputed){
    contrast_results <- contrast_results %>%
      mutate(c1 = case_when(is.na(estimate) ~ c1_imputed, TRUE ~ c1))
    contrast_results <- contrast_results %>%
      mutate(c2 = case_when(is.na(estimate) ~ c2_imputed, TRUE ~ c2))
    contrast_results <- contrast_results %>%
      dplyr::select(-contains("_imputed"))
  }

  separate_hierarchy(contrast_results, config) -> filtered_dd
  return(filtered_dd)
}

#' run the modelling using lmer and lm models
#'
#' @export
#'
application_run_modelling_V2 <- function(outpath,
                                         data,
                                         pepConfig,
                                         modelFunction,
                                         contrasts,
                                         modelling_dir="modelling_results_protein" ,
                                         remove_imputed = TRUE,
                                         do_not_report = "",
                                         DEBUG= FALSE

){
  assign("lfq_write_format", c("xlsx","html"), envir = .GlobalEnv)

  # create result structure
  modelling_path <- file.path(outpath, modelling_dir)
  if(!dir.exists(outpath)){
    dir.create(outpath)
  }
  if(!dir.exists(modelling_path)){
    dir.create(modelling_path)
  }

  ### make modelling  -----
  modellingResult_fun <- workflow_model_analyse(data,
                                                modelFunction,
                                                subject_Id = pepConfig$table$hkeysLevel())
  #################################################
  ### Do missing value imputation
  res_contrasts_imputed <- workflow_missigness_impute_contrasts(data,
                                                                pepConfig,
                                                                contrasts)
  contrasts_xx_imputed <- res_contrasts_imputed("long",what = "all")
  contrasts_xx_imputed <- .columnsImputed(contrasts_xx_imputed,
                                contrasts = contrasts[setdiff(names(contrasts) , do_not_report)])

  #### Compute contrasts from model ####
  modellingResult <-  modellingResult_fun()
  modelProteinF <- modellingResult$modellingResult$modelProtein
  res_contrasts <- workflow_contrasts_linfct_V2(modelProteinF,
                                                contrasts,
                                                pepConfig,
                                                modelName = modelFunction$model_name,
                                                prefix =  "contrasts",
                                                contrastfun = modelFunction$contrast_fun)

  res_fun <- function(do=c("result","write"),DEBuG = FALSE){
    do <- match.arg(do)
    if(DEBUG){
      list(modelFunction = modelFunction,
           imputed = contrasts_xx_imputed,
           remove_imputed = remove_imputed,
           subject_Id = pepConfig$table$hkeysLevel(),
           modelling_path = modelling_path,
           modellingResult_fun = modellingResult_fun,
           res_contrasts = res_contrasts
      )
    }

    if(do == "result"){
      result_table <- .makeResult_contrasts(res_contrasts(columns = modelFunction$report_columns)
                                            ,contrasts_xx_imputed,
                                            pepConfig$table$hkeysLevel())
      return(result_table)
    }else if(do == "write"){
      modellingResult_fun(modelling_path)
      res_contrasts(modelling_path, columns = modelFunction$report_columns)
      lfq_write_table(filtered_dd, path = file.path(modelling_path, "foldchange_estimates.csv"))
    }
  }
  return(res_fun)
}

#' add Annotation to an MQ output
#'
#' for an usage example see run_script lfq_mixed_model_inference
#'
#' @export
#'
application_set_up_MQ_run <- function(outpath,
                                      inputMQfile,
                                      inputAnnotation,
                                      config,
                                      id_extractor = function(df){fgczgseaora::get_UniprotID_from_fasta_header(df, idcolumn = "top_protein")},
                                      qcdir = "qc_results",
                                      use = c("peptides", "modificationSpecificPeptides" )){
  peptides <- match.arg(use)

  assign("lfq_write_format", c("xlsx"), envir = .GlobalEnv)
  if(!is.null(id_extractor)){
    config$table$hierarchy[["protein_Id"]] <- c(config$table$hierarchy[["protein_Id"]], "UniprotID")
  }
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


  mod_peptides_available <- "modificationSpecificPeptides.txt" %in% unzip(inputMQfile, list = TRUE)$Name
  resPepProtAnnot <- NULL
  if(mod_peptides_available){
    resPepProtAnnot <- tidyMQ_modificationSpecificPeptides(inputMQfile)
    { # create visualization for modified peptide sequence
      height <- length(unique(resPepProtAnnot$raw.file))/2 * 300
      png(file.path(qc_path, "retention_time_plot.png"), height = height, width=1200)
      resPepProtVis <- resPepProtAnnot %>% dplyr::filter(mod.peptide.intensity > 4)

      tmp <- ggplot(resPepProtVis, aes(x = retention.time, y= log2(mod.peptide.intensity))) +
        geom_point(alpha=1/20, size=0.3) +
        facet_wrap(~raw.file, ncol=2 )

      print(tmp)
      dev.off()
    }
  }

  peptides_available <- "peptides.txt" %in% unzip(inputMQfile, list = TRUE)$Name
  if(peptides_available & peptides == "peptides"){
    resPepProtAnnot <- tidyMQ_Peptides(inputMQfile)
  }
  else if((peptides_available & mod_peptides_available) & peptides == "modificationSpecificPeptides"){
    peptidestxt <- tidyMQ_Peptides(inputMQfile)
    resPepProtAnnot <- tidyMQ_from_modSpecific_to_peptide(resPepProtAnnot, peptidestxt)
  }else if(!peptides_available & mod_peptides_available){
    warning("no peptides.txt found!! working with modificationSpecificPeptides")
    config$table$workIntensity <- "mod.peptide.intensity"
    config$table$hierarchy[["peptide_Id"]] <- c("sequence", "modifications", "mod.peptide.id")
  }else{
    warning("peptides_available : " ,peptides_available, "peptides : ", peptides)
  }

  resPepProtAnnot <- tidyMQ_top_protein_name(resPepProtAnnot)
  resPepProtAnnot <- resPepProtAnnot %>% dplyr::filter(reverse !=  TRUE)
  resPepProtAnnot$isotope = "light"

  if(!is.null(id_extractor)){
    resPepProtAnnot <- id_extractor(resPepProtAnnot)
  }

  {# add annotation
    if(is.character(inputAnnotation)){
      annotation <- readxl::read_xlsx(inputAnnotation)
    }else{
      annotation <- inputAnnotation
    }
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
    ###  Setup analysis ####
  }


  resDataStart <- setup_analysis(resPepProtAnnot, config)
  resDataStart <- remove_small_intensities( resDataStart, config, threshold = 4 ) %>%
    complete_cases(config)
  return(list(data = resDataStart,
              config=config,
              qc_path = qc_path))
}




#' preprocess peptide data, compute protein data, store results in qc_path folder
#' @export
#'
application_summarize_data_pep_to_prot <- function(data,
                                                   config,
                                                   qc_path,
                                                   DEBUG= FALSE,
                                                   WRITE_PROTS=TRUE){
  message("deprecated use data_pep_to_prot instead")
  res_fun <- data_pep_to_prot(data,
                          config,
                          qc_path)

  if(!DEBUG){res_fun("render")}
  if(WRITE_PROTS){res_fun("plotprot")}
  res_fun("pepwrite")
  res_fun("protwrite")
  return(res_fun(DEBUG=TRUE))
}

#' data_pep_to_prot
#' @export
data_pep_to_prot <- function(data,
                             config,
                             qc_path){
  #qc_path <- qc_path
  results <- LFQService::workflow_MQ_protoV1(
    data,
    config,
    outpath,
    peptideFilterFunction = LFQService:::.workflow_MQ_filter_peptides_V3 )

  ### PROTEIN QUANTIFICATION ####
  protintensity_fun <- medpolish_protein_quants( results$pepIntensityNormalized,
                                                 results$config_pepIntensityNormalized )


  ### generate return value
  res_fun <- function(do = c("render","plotprot","pepwrite","protwrite"),DEBUG=FALSE){
    do <- match.arg(do)
    if(DEBUG){
      return(list(qc_path = qc_path, results = results, protintensity_fun = protintensity_fun ))
    }
    assign("lfq_write_format", c("xlsx"), envir = .GlobalEnv)

    if(do == "render"){
      LFQService::render_MQSummary_rmd(results$filteredPep,
                                       results$config_filteredPep$clone(deep=TRUE),
                                       pep=TRUE,
                                       workdir = ".",
                                       dest_path = qc_path,
                                       dest_file_name="peptide_intensities_qc",
                                       format = "html")
      LFQService::render_MQSummary_rmd(unnestProt$data,
                                       unnestProt$config$clone(deep=TRUE),
                                       pep=FALSE,
                                       workdir = ".",
                                       dest_path = qc_path,
                                       dest_file_name="protein_intensities_qc",
                                       format="html"
      )
    }else if(do == "plotprot"){
      figs <- protintensity_fun("plot")
      pdf(file.path(qc_path, "protein_intensities_inference_figures.pdf"))
      lapply(figs$plot, print)
      dev.off()

    }else if(do =="protwrite"){
      unnestProt <- protintensity_fun("unnest")
      quants_write(unnestProt$data, unnestProt$config, qc_path, na_fraction = 0.3)
    }else if(do == "pepwrite"){
      wideFRAMEPeptide <- LFQService::toWideConfig(results$pepIntensityNormalized,
                                                   results$config_pepIntensityNormalized)
      lfq_write_table(separate_hierarchy(wideFRAMEPeptide$data,
                                         results$config_pepIntensityNormalized),
                      file.path(qc_path, "peptide_intensities.csv"))
    }
  }
  return(res_fun)
}





#' Used for metabolomics data analysis.
#'
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
  quants_write(results$data, results$config, qc_path, prefix = prefix)

  wideFRAME <- LFQService::toWideConfig(results$data,
                                        results$config)
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


