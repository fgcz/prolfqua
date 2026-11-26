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


#' Weight binary detection rows by identification confidence
#'
#' An observed row (detection 1) with identification confidence c and feature-level confidence c_f becomes a
#' detection with weight c and a non-detection with weight c_f - c: a doubtful detection counts only partly. A
#' missing row becomes a non-detection with weight c_f, so a feature that is probably a false identification counts
#' little either way. Rows with weight 0 are dropped; with confidence 1 everywhere, or `confidence` NULL, the rows
#' are unchanged, with weight 1.
#' @param data long-format data with one binary detection row per subject and sample
#' @param bin_resp binary detection column
#' @param confidence identification confidence column, or NULL to count every row once
#' @param confidence_feature feature-level identification confidence column
#' @return `data` with the column `.detection_weight`, observed rows split in two
#' @noRd
.detection_weight_rows <- function(data, bin_resp, confidence, confidence_feature) {
  if (is.null(confidence)) {
    data$.detection_weight <- 1
    return(data)
  }
  observed <- data[data[[bin_resp]] == 1, , drop = FALSE]
  absent <- data[data[[bin_resp]] == 0, , drop = FALSE]
  detected <- observed
  detected$.detection_weight <- observed[[confidence]]
  missed <- observed
  missed[[bin_resp]] <- 0L
  missed$.detection_weight <- observed[[confidence_feature]] - observed[[confidence]]
  absent$.detection_weight <- absent[[confidence_feature]]
  dplyr::filter(dplyr::bind_rows(detected, missed, absent), .data$.detection_weight > 0)
}

#' Fit Firth logistic regression to binary detection data
#'
#' Fits Firth's logistic regression per subject to the binary detection
#' response of aggregated data, one row per subject and sample, weighted by identification confidence (see
#' `.detection_weight_rows()`).
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
  config <- data$get_config()
  rows <- .detection_weight_rows(
    data$data_long(),
    config$bin_resp,
    .used_confidence(config),
    config$identification_confidence_feature
  )
  strategy <- strategy_logistf(formula, control = control, weights = ".detection_weight")
  models <- .model_analyse(
    rows,
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
#' @examplesIf requireNamespace("logistf", quietly = TRUE)
#' strat <- StrategyLogistf$new("bin_resp ~ condition")
#' strat$formula
StrategyLogistf <- R6::R6Class(
  "StrategyLogistf",
  inherit = StrategyBase,
  public = list(
    #' @field weights optional column name with per-row weights; rows with identical predictors and response are
    #'   collapsed into one row whose frequency weight is the sum of their weights
    weights = NULL,

    #' @description Create a new StrategyLogistf
    #' @param modelstr model formula string
    #' @param model_name name of model
    #' @param control named list of convergence controls, see \code{\link{strategy_logistf}}
    #' @param weights optional column name with per-row weights; NULL counts every row once
    initialize = function(modelstr, model_name = "logistf", control = list(), weights = NULL) {
      super$initialize(modelstr, model_name, control = control)
      self$weights <- weights
    },

    #' @description Degrees of freedom of the Wald tests: `Inf`, so contrasts use the normal distribution (z), as
    #'   the coefficient table, `confint()` and logistf itself do
    #' @param model fitted model
    df_residual = function(model) Inf,

    #' @description Check if model is singular (NA coefficients or fewer than two observations beyond the
    #'   coefficients)
    #' @param model fitted model
    isSingular = function(model) anyNA(stats::coefficients(model)) || model$n - length(model$coefficients) < 2,

    #' @description Get residual standard error (always 1 for logistic)
    #' @param model fitted model
    sigma = function(model) 1
  ),
  private = list(
    # Collapse identical predictor rows into frequency weights: the number of rows, or the sum of their weights.
    prepare = function(x) {
      predictor_vars <- all.vars(update(self$formula, . ~ .))
      x$.row_weight <- if (is.null(self$weights)) 1L else x[[self$weights]]
      x |>
        dplyr::group_by(dplyr::across(dplyr::all_of(predictor_vars))) |>
        dplyr::summarize(Freq = sum(.data$.row_weight), .groups = "drop")
    },
    check_packages = function() rlang::check_installed("logistf", reason = "for Firth logistic regression."),
    control_names = function() names(formals(logistf::logistf.control)),
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
strategy_logistf <- function(modelstr, model_name = "logistf", control = list(), weights = NULL) {
  StrategyLogistf$new(modelstr, model_name, control = control, weights = weights)
}
