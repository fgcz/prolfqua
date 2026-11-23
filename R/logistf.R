#' Complete LFQData and encode observed intensities as a binary response
#' @param lfqdata [LFQData] object with the data and its configuration
#' @noRd
.prepare_detection_lfqdata <- function(lfqdata) {
  if (!inherits(lfqdata, "LFQData")) {
    abort_bad_argument("lfqdata", "be an LFQData object")
  }
  lfq_missing <- lfqdata$get_copy()
  lfq_missing$complete_cases()
  lfq_missing$set_data(.encode_bin_resp(lfq_missing))
  lfq_missing$set_config_value("bin_resp", "bin_resp")
  lfq_missing
}


#' Build Firth logistic model for aggregated LFQData
#'
#' Encodes missingness as a binary response and fits the Firth logistic backend
#' used by the missingness model path in \code{prolfquapp}.
#'
#' @param lfqdata aggregated \code{\link{LFQData}} object
#' @param modelstr model formula string without the response variable
#' @param control convergence controls, see \code{\link{strategy_logistf}}
#'   (e.g. \code{"~ group_"})
#' @return a \code{\link{ModelFirth}} object
#' @export
#' @family modelling
#' @examplesIf requireNamespace("logistf", quietly = TRUE)
#' istar <- sim_lfq_data_protein_config(
#'   n_proteins = 10, with_missing = TRUE, weight_missing = 0.5, seed = 3)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' mod <- build_model_glm_protein(lfqdata, "~ group_")
#' head(mod$get_coefficients())
build_model_glm_protein <- function(lfqdata, modelstr, control = list()) {
  .assert_aggregated_facade_input(lfqdata, "build_model_glm_protein")
  lfq_missing <- .prepare_detection_lfqdata(lfqdata)
  formula <- paste(lfq_missing$get_config()$bin_resp, modelstr)
  .build_model_logistf(lfq_missing, formula, control = control)
}


#' Fit Firth logistic regression to binary detection data
#'
#' Fits Firth's logistic regression per subject to the binary detection
#' response of aggregated data, one row per subject and sample.
#' @param data aggregated LFQData with a binary detection response
#' @param formula model formula with the binary response
#' @param control convergence controls, see \code{\link{strategy_logistf}}
#' @return a [ModelFirth] object holding one Firth logistic model per subject
#' @noRd
.build_model_logistf <- function(data, formula, control = list()) {
  if (!inherits(data, "LFQData")) {
    abort_bad_argument("data", "be an LFQData object")
  }
  .assert_aggregated_facade_input(data, "build_model_logistf")
  strategy <- strategy_logistf(formula, control = control)
  models <- .model_analyse(
    data$data_long(),
    strategy,
    model_name = "logistf",
    label = "firth",
    subject_id = data$subject_id()
  )
  models$strategy <- strategy
  model <- ModelFirth$new(models, subject_id = data$subject_id())
  .record_model_fits(model$get_provenance_log(), models$model_df, data$subject_id(), strategy, model$model_name)
  model$get_provenance_log()$report()
  model
}


#' Firth's logistic regression strategy (R6 class)
#'
#' Encapsulates everything needed to fit per-protein Firth's bias-reduced
#' logistic regression via \code{\link[logistf]{logistf}} and extract contrasts.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' strat <- StrategyLogistf$new("bin_resp ~ condition")
#' strat$formula
StrategyLogistf <- R6::R6Class(
  "StrategyLogistf",
  inherit = StrategyBase,
  public = list(
    #' @description Create a new StrategyLogistf
    #' @param modelstr model formula string
    #' @param model_name name of model
    #' @param control named list of convergence controls, see \code{\link{strategy_logistf}}
    initialize = function(modelstr, model_name = "logistf", control = list()) {
      super$initialize(modelstr, model_name, control = control)
    },

    #' @description Get residual degrees of freedom
    #' @param model fitted model
    df_residual = function(model) model$n - length(model$coefficients),

    #' @description Get residual standard error (always 1 for logistic)
    #' @param model fitted model
    sigma = function(model) 1
  ),
  private = list(
    # Collapse identical predictor rows into frequency weights.
    prepare = function(x) {
      rlang::check_installed("logistf", reason = "for Firth logistic regression.")
      predictor_vars <- all.vars(update(self$formula, . ~ .))
      x |>
        dplyr::group_by(dplyr::across(dplyr::all_of(predictor_vars))) |>
        dplyr::summarize(Freq = dplyr::n(), .groups = "drop")
    },
    control_names = function() {
      rlang::check_installed("logistf", reason = "for Firth logistic regression.")
      names(formals(logistf::logistf.control))
    },
    fit = function(DFT) {
      logistf::logistf(
        self$formula,
        data = DFT,
        weights = Freq,
        pl = FALSE,
        control = do.call(logistf::logistf.control, self$control)
      )
    }
  )
)


#' Create Firth's logistic regression strategy
#'
#' Convenience wrapper that creates a \code{\link{StrategyLogistf}} object.
#' @export
#' @rdname strategy
#' @family modelling
#' @return a \code{\link{StrategyLogistf}} object
#' @examplesIf requireNamespace("logistf", quietly = TRUE)
#' tmp <- strategy_logistf("bin_resp ~ condition", model_name = "parallel design")
#' tmp$formula
#'
#' istar <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 10, with_missing = TRUE,
#'   weight_missing = 0.5, seed = 3)
#' lfq <- LFQData$new(istar$data, istar$config)
#' detection <- lfq$data_long()
#' detection$bin_resp <- as.integer(!is.na(detection[[lfq$response()]]))
#' config <- lfq$get_config()
#' config$bin_resp <- "bin_resp"
#' istar <- LFQData$new(detection, config)
#' df <- istar$summarize_hierarchy()
#' df2 <- df[df[[ncol(df)]] > 1, ]
#' istar2 <- istar$get_subset(df2)
#' istar2$data_long() |>
#'   dplyr::group_by(protein_Id) |>
#'   tidyr::nest() -> nestProtein
#' modelFunction <- strategy_logistf("bin_resp ~ group_ + peptide_Id",
#'   model_name = "random_example")
#' modelFunction$model_fun(nestProtein$data[[1]])
#' modelFunction$model_fun(nestProtein$data[[4]])
#' @param control for \code{strategy_logistf()}: arguments of \code{\link[logistf]{logistf.control}}, e.g.
#'   \code{list(maxit = 50)}; default: logistf's defaults (25 iterations, convergence criteria 1e-5 on the
#'   log likelihood, the score and the coefficients), suited to log odds of detection.
strategy_logistf <- function(modelstr, model_name = "logistf", control = list()) {
  StrategyLogistf$new(modelstr, model_name, control = control)
}
