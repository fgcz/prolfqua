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
      res <- res %>% mutate(!!sym(paste0("pseudo_",column)) := case_when(is.na(estimate) ~ estimate, TRUE ~ !!sym(column)))
    }
    res <- res %>% dplyr::select(-imputed, -meanArea)
    return(res)
  }

  contrast_results <- merge_contrasts_results(xx$contrast_minimal, xx_imputed,
                                              subject_Id = pepConfig$table$hkeysLevel(), modelFunction = modelFunction)
  separate_hierarchy(contrast_results, config) -> filtered_dd

  lfq_write_table(filtered_dd, path = file.path(modelling_path, "foldchange_estimates.csv"))
}
