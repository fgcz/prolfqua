# StrategyLimpa -----

#' R6 class for limpa modelling strategy
#'
#' Consumed by \code{\link{build_model_limpa}}. Configures the
#' \code{limpa::voomaLmFitWithImputation} call: formula for the design matrix,
#' trend/robust flags for \code{\link[limma]{eBayes}}, and optional vooma
#' plotting parameters.
#'
#' @export
#' @family modelling
StrategyLimpa <- R6::R6Class(
  "StrategyLimpa",
  public = list(
    #' @field formula model formula
    formula = NULL,
    #' @field model_name name of model
    model_name = character(),
    #' @field trend logical, passed to \code{\link[limma]{eBayes}}
    trend = FALSE,
    #' @field robust logical, passed to \code{\link[limma]{eBayes}}
    robust = FALSE,
    #' @field plot logical, plot the vooma mean-variance trend
    plot = FALSE,
    #' @field span lowess smoother span for vooma trend (NULL = auto)
    span = NULL,
    #' @description Create a new StrategyLimpa
    #' @param modelstr model formula as string (e.g. "abundance ~ group_")
    #' @param model_name name of model
    #' @param trend logical, passed to \code{\link[limma]{eBayes}}
    #' @param robust logical, passed to \code{\link[limma]{eBayes}}
    #' @param plot logical, plot the vooma mean-variance trend
    #' @param span lowess smoother span (NULL = auto)
    initialize = function(modelstr, model_name = "limpa", trend = FALSE, robust = FALSE, plot = FALSE, span = NULL) {
      self$formula <- as.formula(modelstr)
      self$model_name <- model_name
      self$trend <- trend
      self$robust <- robust
      self$plot <- plot
      self$span <- span
    }
  )
)

#' Create limpa modelling strategy
#'
#' Wrapper that returns a \code{\link{StrategyLimpa}} R6 object.
#'
#' @param modelstr model formula as string (e.g. "abundance ~ group_")
#' @param model_name name of model
#' @param trend logical, passed to \code{\link[limma]{eBayes}}
#' @param robust logical, passed to \code{\link[limma]{eBayes}}
#' @param plot logical, plot the vooma mean-variance trend
#' @param span lowess smoother span (NULL = auto)
#' @export
#' @family modelling
strategy_limpa <- function(modelstr, model_name = "limpa", trend = FALSE, robust = FALSE, plot = FALSE, span = NULL) {
  StrategyLimpa$new(modelstr, model_name = model_name, trend = trend, robust = robust, plot = plot, span = span)
}


# build_model_limpa -----

#' Build limpa vooma model from aggregated LFQData
#'
#' Takes protein-level LFQData produced by \code{\link{AggregateLimpa}} (which
#' contains standard error and observation count columns) and fits a
#' limpa vooma model via \code{limpa::voomaLmFitWithImputation}. The standard
#' errors from DPC quantification are used as a bivariate predictor for the
#' vooma precision weights, and the observation counts identify imputed entries
#' for degree-of-freedom correction.
#'
#' Returns a \code{\link{ModelLimma}} object (since the output is a standard
#' \code{MArrayLM}), which can be passed directly to
#' \code{\link{ContrastsLimma}}.
#'
#' @param lfqdata \code{\link{LFQData}} object from \code{\link{AggregateLimpa}}.
#'   Must have \code{config$opt_se} and \code{config$nr_children} set.
#' @param strategy \code{\link{StrategyLimpa}} object
#' @param modelName model name (default from strategy)
#' @return \code{\link{ModelLimma}} object
#' @export
#' @family modelling
build_model_limpa <- function(lfqdata, strategy, modelName = strategy$model_name) {
  if (!requireNamespace("limpa", quietly = TRUE)) {
    stop(
      "Package 'limpa' is required for build_model_limpa. ",
      "Install it from Bioconductor: BiocManager::install('limpa')"
    )
  }

  se_col <- lfqdata$config$opt_se
  nr_col <- lfqdata$config$nr_children
  stopifnot(
    "LFQData must have config$opt_se set (from AggregateLimpa)" = length(se_col) > 0 && nchar(se_col) > 0,
    "LFQData must have config$nr_children set" = length(nr_col) > 0 && nchar(nr_col) > 0
  )

  setup <- .lfqdata_to_elist(lfqdata, strategy$formula)

  # Extract SE matrix for vooma predictor
  se_wide <- lfqdata$to_wide(as.matrix = TRUE, value = se_col)
  predictor <- log(se_wide$data + 1e-6)

  # Extract n_observations for imputed flag
  n_obs_wide <- lfqdata$to_wide(as.matrix = TRUE, value = nr_col)
  imputed <- (n_obs_wide$data == 0)

  # Fit vooma model with imputation-aware precision weights
  fit_args <- list(
    y = setup$elist,
    design = setup$design,
    predictor = predictor,
    imputed = imputed,
    plot = strategy$plot
  )
  if (!is.null(strategy$span)) {
    fit_args$span <- strategy$span
  }
  fit <- do.call(limpa::voomaLmFitWithImputation, fit_args)

  ModelLimma$new(
    fit = fit,
    design = setup$design,
    formula = strategy$formula,
    subject_Id = setup$subject_Id,
    modelName = modelName,
    rowdata = setup$rowdata,
    trend = strategy$trend,
    robust = strategy$robust,
    dummy_model = setup$dummy_model,
    p.adjust = prolfqua::adjust_p_values
  )
}
