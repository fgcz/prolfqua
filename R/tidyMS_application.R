.columnsImputed <- function(all, contrasts) {
  getAST <- function(ee) purrr::map_if(as.list(ee), is.call, getAST)

  get_sides <- function(contrast, all_variables) {
    ast_list <- getAST(rlang::parse_expr(contrast))
    ast_array <- array(as.character(unlist(ast_list)))
    bb <- intersect(gsub("`","",ast_array),all_variables)
    return(bb)
  }

  all_variables <- c(names(contrasts), unique(all$contrast))


  res <- NULL

  for (i in 1:length(contrasts)) {
    cname <- names(contrasts)[i]
    cc <- get_sides(contrasts[i], all_variables)
    if (length(cc) != 2) {
      message("there are ", length(cc) , "> 2 elements")
      next;
    }

    tmp <- all %>% dplyr::filter( contrast %in% c(cname,cc) )
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
.makeResult_contrasts <- function(contrast_minimal,
                                  contrasts_xx_imputed,
                                  subject_Id,
                                  config,
                                  remove_imputed = TRUE ) {

  contrast_results <- dplyr::right_join( contrast_minimal,
                                         contrasts_xx_imputed,
                                         by = c(subject_Id,
                                                "lhs", "c1_name", "c2_name"), suffix = c("","_imputed"))

  contrast_results <- dplyr::rename(contrast_results, contrast = lhs) #
  contrast_results <- contrast_results %>%
    dplyr::mutate(pseudo_estimate = dplyr::case_when(is.na(estimate) ~ estimate_imputed, TRUE ~ estimate))
  contrast_results <- contrast_results %>%
    dplyr::mutate(is_pseudo_estimate = dplyr::case_when(is.na(estimate) ~ TRUE, TRUE ~ FALSE))

  if (remove_imputed) {
    contrast_results <- contrast_results %>%
      dplyr::mutate(c1 = dplyr::case_when(is.na(estimate) ~ c1_imputed, TRUE ~ c1))
    contrast_results <- contrast_results %>%
      dplyr::mutate(c2 = dplyr::case_when(is.na(estimate) ~ c2_imputed, TRUE ~ c2))
    contrast_results <- contrast_results %>%
      dplyr::select(-dplyr::contains("_imputed"))
  }

  separate_hierarchy(contrast_results, config) -> filtered_dd
  return(filtered_dd)
}

#' run the modelling using lmer and lm models
#' @param do_not_report contrasts not to report
#' @export
#' @examples
#'
application_run_modelling_V2 <- function(outpath,
                                         data,
                                         config,
                                         modelFunction,
                                         contrasts,
                                         modelling_dir = "modelling_results_protein" ,
                                         remove_imputed = TRUE,
                                         do_not_report = "",
                                         DEBUG= FALSE)
{
  assign("lfq_write_format", c("xlsx","html"), envir = .GlobalEnv)
  # create result structure
  modelling_path <- file.path(outpath, modelling_dir)
  if (!dir.exists(outpath)) {
    dir.create(outpath)
  }
  if (!dir.exists(modelling_path)) {
    dir.create(modelling_path)
  }

  ### make modeling  -----
  modellingResult_fun <- workflow_model_analyse(data,
                                                modelFunction,
                                                subject_Id = config$table$hkeysDepth())

  ### Do missing value imputation
  res_contrasts_imputed <- workflow_missigness_impute_contrasts(data,
                                                                config,
                                                                contrasts)
  contrasts_xx_imputed <- res_contrasts_imputed("long",what = "all")
  contrasts_xx_imputed <- .columnsImputed(contrasts_xx_imputed,
                                          contrasts = contrasts[setdiff(names(contrasts) ,
                                                                        do_not_report)])

  #### Compute contrasts from model ####

  modellingResult <-  modellingResult_fun()
  modelProteinF <- modellingResult$modellingResult$modelProtein
  res_contrasts <- workflow_contrasts_linfct_V2(modelProteinF,
                                                contrasts,
                                                config,
                                                modelName = modelFunction$model_name,
                                                prefix =  "contrasts",
                                                contrastfun = modelFunction$contrast_fun)

  # RESULT FUNCTION
  res_fun <- function(do = c("result",
                             "write_modelling",
                             "write_contrasts"),
                      DEBUG = FALSE,
                      remove_imputed = TRUE) {

    do <- match.arg(do)
    if (DEBUG) {
      res <- list(modelFunction = modelFunction,
                  imputed = contrasts_xx_imputed,
                  remove_imputed = remove_imputed,
                  subject_Id = config$table$hkeysDepth(),
                  modelling_path = modelling_path,
                  modellingResult_fun = modellingResult_fun,
                  res_contrasts = res_contrasts
      )
      return(res)
    }

    if (do == "result") {
      contrast_minimal <- res_contrasts(columns = modelFunction$report_columns)$contrast_minimal
      result_table <- .makeResult_contrasts(contrast_minimal
                                            ,contrasts_xx_imputed,
                                            config$table$hkeysDepth(),
                                            config,
                                            remove_imputed = remove_imputed)

      return(result_table)
    } else if (do == "write_modelling") {
      modellingResult_fun(modelling_path)
    } else if (do == "write_contrasts") {
      filtered_dd <- res_contrasts(modelling_path, columns = modelFunction$report_columns)
      result_table <- .makeResult_contrasts(filtered_dd$contrast_minimal
                                            ,contrasts_xx_imputed,
                                            config$table$hkeysDepth(),
                                            config)
      lfq_write_table(result_table, path = file.path(modelling_path, "foldchange_estimates.csv"))
      return(result_table)
    }
  }
  return(res_fun)
}




#' add Annotation to an MQ output
#'
#' for an usage example see run_script lfq_mixed_model_inference
#' @param use use either peptide txt or mod specific peptide
#' @export
#'
application_add_annotation <- function(intensityData,
                                       inputAnnotation,
                                       fileName = "raw.file") {
  assign("lfq_write_format", c("xlsx"), envir = .GlobalEnv)
  ## read the data
  {# add annotation
    if ( is.character(inputAnnotation) ) {
      annotation <- readxl::read_xlsx(inputAnnotation)
    } else {
      annotation <- inputAnnotation
    }
    noData <- annotation[!annotation[[fileName]] %in% intensityData[[fileName]],]
    if (nrow(noData)) {
      message("some files in annotation have no measurements")
      message(paste(noData, collapse = " "))
    }
    measSamples <- unique(intensityData[[fileName]])
    noAnnot <- measSamples[!measSamples %in% annotation[[fileName]] ]
    if (length(noAnnot) > 0 ) {
      message("some measured samples have no annotation!")
      message(paste(noAnnot,collapse = " "))
    }
    resPepProtAnnot <- inner_join(annotation, intensityData, by = fileName)
    ###  Setup analysis ####
  }
  return(resPepProtAnnot)
}



#' Used for metabolomics data analysis.
#'
#' @export
#'
application_summarize_compound <- function(data,
                                           config,
                                           qc_path,
                                           prefix = c("ms")) {
  qc_apth <- qc_path
  prefix <- match.arg(prefix)
  assign("lfq_write_format", c("xlsx"), envir = .GlobalEnv)

  results <- LFQService:::normalize_log2_robscale(data, config)


  res_fun <- function(do = c("plot", "write", "render", "print_compounds", "data"),
                      DEBUG = FALSE){
    do <- match.arg(do)
    if (DEBUG) {
      return(list(qc_path = qc_path, prefix = prefix, results = results ))
    }

    if ( do == "plot") {
      quants_write(results$data, results$config, qc_path)
    }else if (do == "render") {
      LFQService::render_MQSummary_rmd(results$data,
                                       results$config$clone(deep = TRUE),
                                       pep = TRUE,
                                       workdir = ".",
                                       dest_path = qc_path,
                                       dest_file_name = paste0(prefix, "_intensities_qc"),
                                       format = "html")
    }else if (do == "print_compounds") {
      figs <- plot_hierarchies_line_df(results$data, results$config )
      pdf(file.path(qc_path, paste0(prefix, "_intensities_inference_figures.pdf")))
      lapply(figs, print)
      dev.off()
    }else if (do == "data") {
      return(list(data  = results$data, config = results$config))
    }

  }
  return(res_fun)

}

#################################################
### Do missing value imputation


