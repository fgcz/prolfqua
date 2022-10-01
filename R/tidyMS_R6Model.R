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
#' data_Yeast2Factor <- prolfqua::old2new(prolfqua::prolfqua_data("data_Yeast2Factor"))
#' pMerged <- LFQData$new(data_Yeast2Factor$data, data_Yeast2Factor$config)
#'
#' pMerged$data$Run_ID <- as.numeric(pMerged$data$Run_ID)
#' pMerged$config$table$get_work_intensity()
#' pMerged$factors()
#'
#' formula_condition_and_Batches <-
#'   prolfqua::strategy_lm("transformedIntensity ~ condition_ + batch_")
#' modCB <- prolfqua::build_model(
#'   pMerged$data,
#'   formula_condition_and_Batches,
#'   subject_Id = pMerged$config$table$hierarchyKeys() )
#'
#' formula_condition <-
#'   prolfqua::strategy_lm("transformedIntensity ~ condition_")
#' modC <- prolfqua::build_model(
#'   pMerged$data,
#'   formula_condition,
#'   subject_Id = pMerged$config$table$hierarchyKeys() )
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
#' # library(tidyverse)
#' D <- old2new(prolfqua_data('data_ionstar')$normalized())
#' D$data <- dplyr::filter(D$data ,protein_Id %in% sample(protein_Id, 100))
#'
#' modelName <- "f_condtion_r_peptide"
#' formula_randomPeptide <-
#'   strategy_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id) + (1 | sampleName)",
#'    model_name = modelName)
#'
#'
#' pepIntensity <- D$data
#' config <- D$config
#'
#'
#'
#' mod <- prolfqua:::build_model(
#'  pepIntensity,
#'  formula_randomPeptide,
#'  modelName = modelName,
#'  subject_Id = config$table$hkeysDepth())
#'
#' mod <- prolfqua:::build_model(
#'  LFQData$new(pepIntensity, config),
#'  formula_randomPeptide,
#'  modelName = modelName)
#' model_summary(mod)
#'
build_model <- function(data,
                        modelFunction,
                        subject_Id = if ("LFQData" %in% class(data)) {data$subject_Id()} else {"protein_Id"},
                        modelName = modelFunction$model_name){

  dataX <- if ("LFQData" %in% class(data)) { data$data }else{ data }
  modellingResult <- model_analyse(dataX,
    modelFunction,
    modelName = modelName,
    subject_Id = subject_Id)
  return( Model$new(modelDF = modellingResult$modelProtein,
                    modelFunction = modelFunction,
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
