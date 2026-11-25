# StrategyLimpa -----

#' R6 class for limpa modelling strategy
#'
#' Consumed by \code{\link{build_model_limpa}}. Configures the
#' \code{limpa::voomaLmFitWithImputation} call: formula for the design matrix,
#' trend/robust flags for \code{\link[limma]{eBayes}}, and optional vooma
#' plotting parameters.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' strat <- StrategyLimpa$new("abundance ~ group_")
#' strat$formula
#' strat$model_name
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
    #' @field weights value column names whose product gives the prior weights, or NULL
    weights = NULL,
    #' @description Create a new StrategyLimpa
    #' @param modelstr,model_name,trend,robust,plot,span,weights see \code{\link{strategy_limpa}}
    initialize = function(
      modelstr,
      model_name = "limpa",
      trend = FALSE,
      robust = FALSE,
      plot = FALSE,
      span = NULL,
      weights = NULL
    ) {
      self$formula <- .as_model_formula(modelstr)
      self$model_name <- .assert_string(model_name)
      self$trend <- .assert_flag(trend)
      self$robust <- .assert_flag(robust)
      self$plot <- .assert_flag(plot)
      self$span <- span
      self$weights <- weights
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
#' @param weights value column names (e.g. the identification confidence) whose product is passed as
#'   \code{prior.weights}; cells without a value (imputed) get weight 1. Default \code{NULL} (no prior weights).
#' @return a [StrategyLimpa] object holding the formula, model name, `eBayes()` options, plot flag, span and weights.
#' @export
#' @family modelling
#' @examples
#' strat <- strategy_limpa("abundance ~ group_")
#' strat$formula
#' strat$model_name
strategy_limpa <- function(
  modelstr,
  model_name = "limpa",
  trend = FALSE,
  robust = FALSE,
  plot = FALSE,
  span = NULL,
  weights = NULL
) {
  StrategyLimpa$new(
    modelstr,
    model_name = model_name,
    trend = trend,
    robust = robust,
    plot = plot,
    span = span,
    weights = weights
  )
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
#' @param lfqdata \code{\link{LFQData}} object. Must have
#'   \code{config$nr_children} set (used to flag imputed observations via
#'   \code{nr_children == 0}). \code{config$opt_se} is optional: when set
#'   (e.g. from \code{\link{AggregateLimpa}}), the SE column is used as a
#'   vooma precision-weight predictor; otherwise vooma is fit without an
#'   external predictor.
#' @param strategy \code{\link{StrategyLimpa}} object
#' @param model_name model name (default from strategy)
#' @return \code{\link{ModelLimma}} object
#' @export
#' @family modelling
#' @examples
#' if (requireNamespace("limpa", quietly = TRUE)) {
#'   istar <- sim_lfq_data_peptide_config(n_proteins = 10)
#'   lfq <- LFQData$new(istar$data, istar$config)
#'   lfq <- lfq$get_Transformer()$log2()$lfq
#'   lfq <- AggregateLimpa$new(lfq, "protein")$aggregate()
#'   strat <- strategy_limpa("limpa ~ group_")
#'   mod <- build_model_limpa(lfq, strat)
#'   mod$model_name
#' }
build_model_limpa <- function(lfqdata, strategy, model_name = strategy$model_name) {
  rlang::check_installed("limpa", reason = "for build_model_limpa().")
  .assert_model_input(lfqdata)

  se_col <- lfqdata$get_config()$opt_se
  nr_col <- lfqdata$nr_children_col()
  if (!(length(nr_col) > 0 && nchar(nr_col) > 0)) {
    abort_invalid_config("`lfqdata` must have `config$nr_children` set.")
  }

  setup <- .lfqdata_to_elist(lfqdata, strategy$formula)

  predictor <- NULL
  if (length(se_col) > 0 && nchar(se_col) > 0) {
    se_wide <- lfqdata$data_wide(as.matrix = TRUE, value = se_col)
    predictor <- log(se_wide$data + 1e-6)
  }

  n_obs_wide <- lfqdata$data_wide(as.matrix = TRUE, value = nr_col)
  imputed <- (n_obs_wide$data == 0)

  # Limpa's voomaLmFitWithImputation requires NA-free input. AggregateLimpa
  # pre-imputes; other aggregators (e.g. medpolish) do not. Impute any
  # remaining NAs via limpa::imputeByExpTilt and mark them as imputed.
  na_mask <- is.na(setup$elist$E)
  if (any(na_mask)) {
    setup$elist <- limpa::imputeByExpTilt(setup$elist)
    imputed <- imputed | na_mask
  }

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
  if (!is.null(strategy$weights)) {
    # imputed cells have no identification of their own and get weight 1
    prior_weights <- .resolve_weights(lfqdata, strategy, setup$annotation)
    prior_weights[is.na(prior_weights)] <- 1
    fit_args$prior.weights <- prior_weights
  }
  fit <- do.call(limpa::voomaLmFitWithImputation, fit_args)
  .new_model_limma(fit, setup, strategy, model_name, fit$EList$E, fit$EList$weights, !imputed)
}
