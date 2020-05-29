#' run the modelling using lmer and lm models - DEPRECATED use version V2
#'
#' @keywords internal
#' @noRd
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
                                                subject_Id = pepConfig$table$hkeysDepth())

  modellingResult <- modellingResult_fun()
  modellingResult_fun(modelling_path)

  #names(modellingResult)

  m <- get_complete_model_fit(modellingResult$modellingResult$modelProtein)

  #factor_contrasts <- linfct_factors_contrasts(m$linear_model[[1]])
  linfct <- linfct_from_model(m$linear_model[[1]], as_list = FALSE)
  linfct_A <- linfct_matrix_contrasts(linfct, contrasts)

  if (DEBUG) {
    pdf(file.path(modelling_path,"Linear_functions.pdf"), width = 18, height = 10)
    #quantable::imageWithLabels(t(linfct_A),
    #                           col = quantable::getBlueWhiteRed(),
    #                           marLeft = c(8,10,4.1,2.1))
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
    res <- right_join(contrast_minimal, xx_imputed, by = c(subject_Id,"lhs" = "contrast"))
    res <- res %>% dplyr::rename(contrast = lhs)
    res <- res %>% dplyr::mutate(pseudo_estimate = case_when(is.na(estimate) ~ imputed, TRUE ~ estimate))
    res <- res %>% dplyr::mutate(is_pseudo_estimate = case_when(is.na(estimate) ~ TRUE, TRUE ~ FALSE))

    for (column in modelFunction$report_columns) {
      res <- res %>% dplyr::mutate(!!sym(paste0("pseudo_",column)) := case_when(is.na(estimate) ~ estimate,
                                                                                TRUE ~ !!sym(column)))
    }
    res <- res %>% dplyr::select(-imputed, -meanArea)
    return(res)
  }

  contrast_results <- merge_contrasts_results(xx$contrast_minimal, xx_imputed,
                                              subject_Id = pepConfig$table$hkeysDepth(), modelFunction = modelFunction)
  separate_hierarchy(contrast_results, config) -> filtered_dd

  lfq_write_table(filtered_dd, path = file.path(modelling_path, "foldchange_estimates.csv"))
}


#' Do contrast
#'
#' @keywords internal
#' @noRd
workflow_contrasts_linfct <- function(models,
                                      contrasts,
                                      config,
                                      modelName = "Model",
                                      prefix = "Contrasts",
                                      contrastfun = LFQService::my_contest )
{
  warning("DEPRECATE workflow_contrasts_linfct!\n use workflow_contrasts_linfct_V2")
  if (class(contrasts) == "matrix") {
    linfct_A <- contrasts
  }else{
    models <- models %>% dplyr::filter(exists_lmer == TRUE)
    m <- get_complete_model_fit(models)
    linfct <- linfct_from_model(m$linear_model[[1]], as_list = FALSE)
    linfct_A <- linfct_matrix_contrasts(linfct, contrasts)
  }

  subject_Id <- config$table$hkeysDepth()
  contrast_result <- contrasts_linfct(models,
                                      linfct_A,
                                      subject_Id = subject_Id,
                                      contrastfun = contrastfun )

  contrast_result <- moderated_p_limma_long(contrast_result)
  subject_Id <- subject_Id
  prefix <- prefix
  modelName <- modelName

  res_fun <- function(path = NULL, columns = c("p.value",
                                               "p.value.adjusted",
                                               "moderated.p.value",
                                               "moderated.p.value.adjusted"),
                      DEBUG = FALSE){
    if (DEBUG) {
      return(list(contrast_result = contrast_result,
                  modelName = modelName,
                  config = config,
                  prefix = prefix,
                  subject_Id = subject_Id,
                  columns = columns
      ))
    }

    visualization <- contrasts_linfct_vis(contrast_result,
                                          modelName ,
                                          prefix = prefix,
                                          subject_Id = subject_Id,
                                          columns = columns
    )

    relevant_columns <- c("lhs",
                          "sigma",
                          "df",
                          "isSingular",
                          "estimate",
                          "conf.low",
                          "conf.high") # other relevant columns.
    contrast_minimal <- contrast_result %>% dplyr::select(subject_Id, relevant_columns, columns )

    contrasts_wide <- pivot_model_contrasts_2_Wide(contrast_minimal,
                                                   subject_Id = subject_Id,
                                                   columns = c("estimate", columns))

    if (!is.null(path)) {
      if (FALSE) {
        contrasts_linfct_write(contrast_minimal,
                               config,
                               path = path,
                               modelName = modelName,
                               prefix = prefix,
                               columns = c("estimate", columns))
      }
      contrasts_linfct_vis_write(visualization, path = path, format = "pdf")
      contrasts_linfct_vis_write(visualization, path = path, format = "html")
    }

    res <- list(contrast_result = contrast_result,
                contrast_minimal = contrast_minimal,
                contrasts_wide = contrasts_wide,
                visualization = visualization,
                modelName = modelName,
                prefix = prefix)

    invisible(res)
  }
  return( res_fun )
}


#' preprocess peptide data, compute protein data, store results in qc_path folder
#'
#' @keywords internal
#' @noRd
application_summarize_data_pep_to_prot <- function(data,
                                                   config,
                                                   qc_path,
                                                   DEBUG= FALSE,
                                                   WRITE_PROTS=TRUE) {

  message("deprecated use data_pep_to_prot instead")
  res_fun <- data_pep_to_prot(data,
                              config,
                              qc_path)

  if (!DEBUG) {res_fun("render")}
  if (WRITE_PROTS) {res_fun("plotprot")}
  res_fun("pepwrite")
  res_fun("protwrite")
  return(res_fun(DEBUG = TRUE))
}
