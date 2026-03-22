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
#' pMerged$config$get_response()
#' pMerged$factors()
#'
#' formula_condition_and_Batches <-
#'   prolfqua::strategy_lm("abundance ~ Treatment + Background")
#' modCB <- prolfqua::build_model(
#'   pMerged$data,
#'   formula_condition_and_Batches,
#'   subject_Id = pMerged$config$hierarchy_keys() )
#'
#' formula_condition <-
#'   prolfqua::strategy_lm("abundance ~ Treatment")
#' modC <- prolfqua::build_model(
#'   pMerged$data,
#'   formula_condition,
#'   subject_Id = pMerged$config$hierarchy_keys() )
#'
#' tmp <- LR_test(modCB$modelDF, "modCB", modC$modelDF, "modB")
#' hist(tmp$likelihood_ratio_test.pValue)
#'
LR_test <- function(
  modelProteinF,
  modelName,
  modelProteinF_Int,
  modelName_Int,
  subject_Id = "protein_Id",
  path = NULL
) {
  # Model Comparison
  reg <- dplyr::inner_join(
    dplyr::select(modelProteinF, !!sym(subject_Id), "linear_model"),
    dplyr::select(modelProteinF_Int, !!sym(subject_Id), "linear_model"),
    by = subject_Id
  )

  reg <- reg |>
    dplyr::mutate(
      modelComparisonLikelihoodRatioTest = map2(
        !!sym("linear_model.x"),
        !!sym("linear_model.y"),
        .likelihood_ratio_test
      )
    )
  likelihood_ratio_test_result <- reg |>
    dplyr::select(!!sym(subject_Id), dplyr::all_of("modelComparisonLikelihoodRatioTest")) |>
    tidyr::unnest(cols = c("modelComparisonLikelihoodRatioTest"))
  likelihood_ratio_test_result <- likelihood_ratio_test_result |>
    dplyr::rename(likelihood_ratio_test.pValue = .data$modelComparisonLikelihoodRatioTest)

  if (!is.null(path)) {
    fileName <- paste("hist_LRT_", modelName, "_", modelName_Int, ".pdf", sep = "")
    fileName <- file.path(path, fileName)
    message("writing figure : ", fileName, "\n")
    pdf(fileName)
    par(mfrow = c(2, 1))
    hist(likelihood_ratio_test_result$likelihood_ratio_test.pValue, breaks = 20)
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
#' @param data data - a data frame or LFQData object
#' @param model_strategy model strategy object (e.g. from strategy_lmer or strategy_lm)
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
#'  subject_Id = D$config$hierarchy_keys_depth())
#' aovtable <- mod$get_anova()
#'
#' mod <- prolfqua::build_model(
#'  LFQData$new(D$data, D$config),
#'  formula_randomPeptide,
#'  modelName = modelName)
#' model_summary(mod)
#'
#'
build_model <- function(
  data,
  model_strategy,
  subject_Id = if ("LFQData" %in% class(data)) {
    data$subject_Id()
  } else {
    "protein_Id"
  },
  modelName = model_strategy$model_name
) {
  dataX <- if ("LFQData" %in% class(data)) data$data else data
  modellingResult <- model_analyse(dataX, model_strategy, modelName = modelName, subject_Id = subject_Id)
  return(Model$new(
    modelDF = modellingResult$modelDF,
    model_strategy = model_strategy,
    modelName = modellingResult$modelName,
    subject_Id = subject_Id
  ))
}


#' Build protein models with LOD imputation for failed fits
#'
#' Fits per-protein models, then re-fits failed/singular proteins after
#' imputing missing values with the limit of detection (LOD) and clamping.
#' Covariance is borrowed from successful fits so that variance is not
#' underestimated by the constant imputation.
#'
#' @param lfqdata LFQData object (aggregated to protein level)
#' @param model_strategy model strategy object (e.g. from strategy_lm)
#' @param modelName model name (default appends "Imputed")
#' @param lod numeric limit of detection; if NULL, auto-computed from data
#' @param borrow_method "sigma" borrows scalar sigma and uses per-protein
#'   (X'X)^-1; "vcov" borrows element-wise median of full vcov matrices
#' @param df_method "observed" uses max(n_observed - p, 1);
#'   "borrowed" uses median df from successful fits
#' @return a object of class \code{\link{Model}}
#' @family modelling
#' @seealso \code{\link{build_model}}, \code{\link{impute_refit_singular}}
#' @export
#' @examples
#' istar <- sim_lfq_data_protein_config(Nprot = 30, weight_missing = 0.5)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$rename_response("transformedIntensity")
#' strat <- strategy_lm(paste(lfqdata$config$get_response(), "~ group_"))
#' mod <- build_model_impute(lfqdata, strat)
#'
build_model_impute <- function(
  lfqdata,
  model_strategy,
  modelName = paste0(model_strategy$model_name, "Imputed"),
  lod = NULL,
  borrow_method = c("sigma", "vcov"),
  df_method = c("observed", "borrowed")
) {
  borrow_method <- match.arg(borrow_method)
  df_method <- match.arg(df_method)
  subject_Id <- lfqdata$subject_Id()
  response <- lfqdata$config$get_response()

  modellingResult <- model_analyse(
    lfqdata$data,
    model_strategy,
    modelName = modelName,
    subject_Id = subject_Id
  )

  if (is.null(lod)) {
    mh <- MissingHelpers$new(lfqdata$data, lfqdata$config)
    lod <- mh$get_LOD()
  }

  modellingResult$modelDF <- impute_refit_singular(
    modellingResult$modelDF,
    model_strategy,
    lod = lod,
    response = response,
    borrow_method = borrow_method,
    df_method = df_method
  )

  return(Model$new(
    modelDF = modellingResult$modelDF,
    model_strategy = model_strategy,
    modelName = modellingResult$modelName,
    subject_Id = subject_Id
  ))
}


#' Summarize modelling and error reporting
#' @param mod model table see \code{\link{build_model}}
#' @keywords internal
#' @family modelling
#' @export
#' @examples
#' D <- prolfqua::sim_lfq_data_peptide_config(Nprot = 20, weight_missing = 0.1)
#' formula_rp <- strategy_lmer("abundance ~ group_ + (1 | peptide_Id) + (1 | sampleName)")
#' mod <- prolfqua::build_model(
#'   LFQData$new(D$data, D$config), formula_rp)
#' res <- model_summary(mod)
#' stopifnot(is.list(res))
#' stopifnot(all(c("exists", "isSingular") %in% names(res)))
model_summary <- function(mod) {
  res <- list()
  res$exists <- table(mod$modelDF$exists_lmer)
  res$isSingular <- table(mod$modelDF$isSingular)
  return(res)
}
