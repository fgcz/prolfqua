#' Likelihood ratio test
#' @family modelling
#' @export
#' @param modelProteinF table with models (see build model)
#' @param modelName name of model
#' @param modelProteinF_Int reduced model
#' @param modelName_Int name of reduced model
#' @param subject_Id subject id typically Assession or protein_Id
#' @param path default NULL, set to a directory if you need to write diagnostic plots.
#' @examples
#' data_2Factor <- prolfqua::sim_lfq_data_2Factor_config(
#'  Nprot = 200,
#'  with_missing = TRUE,
#'  weight_missing = 2)
#'
#' pMerged <- LFQData$new(data_2Factor$data, data_2Factor$config)
#'
#' pMerged$config$table$get_response()
#' pMerged$factors()
#'
#' formula_condition_and_Batches <-
#'   prolfqua::strategy_lm("abundance ~ Treatment + Background")
#' modCB <- prolfqua::build_model(
#'   pMerged$data,
#'   formula_condition_and_Batches,
#'   subject_Id = pMerged$config$table$hierarchy_keys() )
#'
#' formula_condition <-
#'   prolfqua::strategy_lm("abundance ~ Treatment")
#' modC <- prolfqua::build_model(
#'   pMerged$data,
#'   formula_condition,
#'   subject_Id = pMerged$config$table$hierarchy_keys() )
#'
#' tmp <- LR_test(modCB$modelDF, "modCB", modC$modelDF, "modB")
#' hist(tmp$likelihood_ratio_test.pValue)
#'
LR_test <- function(modelProteinF,
                    modelName,
                    modelProteinF_Int,
                    modelName_Int,
                    subject_Id = "protein_Id",
                    path = NULL
){
  # Model Comparison
  reg <- dplyr::inner_join(dplyr::select(modelProteinF, !!sym(subject_Id), "linear_model"),
                           dplyr::select(modelProteinF_Int, !!sym(subject_Id), "linear_model") , by = subject_Id)

  reg <- reg |> dplyr::mutate(modelComparisonLikelihoodRatioTest = map2(!!sym("linear_model.x"),
                                                                        !!sym("linear_model.y"),
                                                                        .likelihood_ratio_test ))
  likelihood_ratio_test_result <- reg |>
    dplyr::select(!!sym(subject_Id), .data$modelComparisonLikelihoodRatioTest) |>
    tidyr::unnest(cols = c("modelComparisonLikelihoodRatioTest"))
  likelihood_ratio_test_result <- likelihood_ratio_test_result |>
    dplyr::rename(likelihood_ratio_test.pValue = .data$modelComparisonLikelihoodRatioTest)


  if (!is.null(path)) {
    fileName <- paste("hist_LRT_", modelName, "_", modelName_Int, ".pdf", sep = "")
    fileName <- file.path(path, fileName)
    message("writing figure : " , fileName , "\n")
    pdf(fileName)
    par(mfrow = c(2, 1))
    hist(likelihood_ratio_test_result$likelihood_ratio_test.pValue,
         breaks = 20)
    plot(ecdf(
      likelihood_ratio_test_result$likelihood_ratio_test.pValue
    ))
    abline(v = c(0.01, 0.05), col = c(3, 2))
    dev.off()
  }

  return(likelihood_ratio_test_result)
}


#' Build protein models from data
#'
#'
#'
#' @param data data - a data frame
#' @param modelFunction model function
#' @param subject_Id grouping variable
#' @param modelName model name
#' @return
#' a object of class \code{\link{Model}}
#' @family modelling
#' @seealso \code{\link{model_analyse}}, \code{\link{strategy_lmer}} \code{\link{strategy_lm}}
#'
#' @export
#' @examples
#' D <- prolfqua::sim_lfq_data_peptide_config(Nprot = 20, weight_missing = 0.1)
#' D$data$abundance |> is.na() |> sum()
#' D <- prolfqua::sim_lfq_data_peptide_config(Nprot = 20, weight_missing = 0.1, seed =3)
#' D$data$abundance |> is.na() |> sum()
#' modelName <- "f_condtion_r_peptide"
#' formula_randomPeptide <-
#'   strategy_lmer("abundance  ~ group_ + (1 | peptide_Id) + (1 | sampleName)",
#'    model_name = modelName)
#'
#'
#' mod <- prolfqua::build_model(
#'  D$data,
#'  formula_randomPeptide,
#'  modelName = modelName,
#'  subject_Id = D$config$table$hierarchy_keys_depth())
#' aovtable <- mod$get_anova()
#'
#' mod <- prolfqua::build_model(
#'  LFQData$new(D$data, D$config),
#'  formula_randomPeptide,
#'  modelName = modelName)
#' model_summary(mod)
#'
#'
build_model <- function(data,
                        model_strategy,
                        subject_Id = if ("LFQData" %in% class(data)) {data$subject_Id()} else {"protein_Id"},
                        modelName = model_strategy$model_name){

  dataX <- if ("LFQData" %in% class(data)) { data$data }else{ data }
  modellingResult <- model_analyse(dataX,
                                   model_strategy,
                                   modelName = modelName,
                                   subject_Id = subject_Id)
  return( Model$new(modelDF = modellingResult$modelProtein,
                    model_strategy = model_strategy,
                    modelName = modellingResult$modelName,
                    subject_Id = subject_Id))
}





#' Summarize modelling and error reporting
#' @param mod model table see \code{\link{build_model}}
#' @keywords internal
#' @family modelling
#' @export
model_summary <- function(mod){
  res <- list()
  res$exists <- table(mod$modelDF$exists_lmer)
  res$isSingular <- table(mod$modelDF$isSingular)
  return(res)

}
