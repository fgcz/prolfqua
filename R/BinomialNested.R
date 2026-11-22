# Binomial detection model ----

#' Count detected and undetected child features per parent and sample
#' @param lfqdata [LFQData] object with the data and its configuration
#' @noRd
.summarize_detection_counts <- function(lfqdata) {
  if (!inherits(lfqdata, "LFQData")) {
    abort_bad_argument("lfqdata", "be an LFQData object")
  }
  bin_resp <- lfqdata$get_config()$bin_resp
  grouping <- unique(c(
    lfqdata$subject_id(),
    lfqdata$sample_name(),
    lfqdata$factor_keys()
  ))

  counts <- lfqdata$data_long() |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping))) |>
    dplyr::summarise(
      detected = sum(.data[[bin_resp]]),
      available = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(undetected = .data$available - .data$detected)

  if (!isTRUE(all(counts$detected >= 0L) && all(counts$undetected >= 0L))) {
    abort_prolfqua("Internal error: detection counts are negative.")
  }
  counts
}

#' Validate the symmetric binomial pseudo-count
#' @param prior_count non-negative symmetric detection pseudo-count
#' @noRd
.validate_binomial_prior_count <- function(prior_count) {
  candidate <- if (is.numeric(prior_count)) prior_count else NA_real_
  is_valid <- all(c(
    length(candidate) == 1L,
    !anyNA(candidate),
    is.finite(candidate),
    candidate >= 0
  ))
  if (!is_valid) {
    abort_bad_argument("prior_count", "be one non-negative number")
  }
  prior_count
}


#' Quasibinomial detection-count strategy
#'
#' Fits detected and undetected child-feature counts for each parent feature
#' using \code{\link[stats]{glm}} with a quasibinomial family. The symmetric
#' pseudo-count stabilizes fits under complete separation; it is not equivalent
#' to Firth's bias-reducing penalty.
#'
#' @return An R6 class generator.
#' @export
#' @examples
#' dat <- data.frame(
#'   group_ = factor(rep(c("A", "B"), each = 4)),
#'   detected = c(1, 2, 1, 3, 4, 5, 3, 5),
#'   undetected = c(4, 3, 4, 2, 1, 0, 2, 0)
#' )
#' strategy <- StrategyBinomial$new("~ group_")
#' fit <- strategy$model_fun(dat)
#' coefficients(fit)
#' @srrstats {G5.8d} Tests reject inappropriate hierarchy, rank deficiency, insufficient residual df and too few PCA
#'   samples.
StrategyBinomial <- R6::R6Class(
  "StrategyBinomial",
  inherit = StrategyBase,
  public = list(
    #' @field prior_count symmetric pseudo-count added to both outcomes
    prior_count = NULL,

    #' @description Create a quasibinomial count strategy. \code{sigma()} is the
    #'   Pearson residual scale used by \code{vcov()}.
    #' @param modelstr right-hand-side model formula, for example
    #'   \code{"~ group_"}
    #' @param prior_count non-negative symmetric pseudo-count
    #' @param model_name model identity
    #' @param control named list of convergence controls, see \code{\link{strategy_binomial}}
    initialize = function(modelstr, prior_count = 0.1, model_name = "binomial_nested", control = list()) {
      prior_count <- .validate_binomial_prior_count(prior_count)
      super$initialize(paste("cbind(.detected, .undetected)", modelstr), model_name, control = control)
      self$prior_count <- prior_count
    },

    #' @description Never flags: a zero Pearson dispersion (identical detection pattern within groups) is bounded
    #'   below by the binomial variance (\code{binomial_bound}), which keeps the test defined
    #' @param sigma Pearson residual scales
    #' @param scale unused
    zero_residual_variance = function(sigma, scale) rep(FALSE, length(sigma))
  ),
  private = list(
    control_defaults = function() list(epsilon = 1e-8, maxit = 25),
    # controls of the inner glm fit
    glm_control = function() self$control,
    prepare = function(x) {
      x$.detected <- x$detected + self$prior_count
      x$.undetected <- x$undetected + self$prior_count
      x
    },
    fit = function(x) {
      model <- stats::glm(
        self$formula,
        data = x,
        family = stats::quasibinomial(),
        control = do.call(stats::glm.control, private$glm_control())
      )
      if (anyNA(stats::coef(model)) || model$rank < length(stats::coef(model))) {
        abort_prolfqua("The binomial design is rank-deficient.")
      }
      if (stats::df.residual(model) < 2L) {
        abort_prolfqua("The binomial model needs at least two residual degrees of freedom.")
      }
      model
    }
  )
)


#' Create a quasibinomial detection-count strategy
#'
#' @param modelstr right-hand-side model formula
#' @param prior_count non-negative symmetric pseudo-count
#' @param model_name model identity
#' @param control \code{epsilon} (default 1e-8) and \code{maxit} (default 25) of \code{\link[stats]{glm.control}}:
#'   the relative deviance change that stops IRLS and its iteration limit; relative, so independent of the
#'   number of detections, and glm's own defaults.
#' @return A \code{\link{StrategyBinomial}} object.
#' @export
#' @examples
#' strategy <- strategy_binomial("~ group_", prior_count = 0.1)
#' strategy$formula
strategy_binomial <- function(modelstr, prior_count = 0.1, model_name = "binomial_nested", control = list()) {
  StrategyBinomial$new(modelstr, prior_count, model_name, control = control)
}


#' Firth-penalized quasibinomial detection-count strategy
#'
#' Fits detected and undetected child-feature counts like
#' \code{\link{StrategyBinomial}}, but replaces the fixed pseudo-count with
#' Firth's bias reduction: each count row gets \eqn{h/2} added to both outcomes,
#' where \eqn{h} is its hat value at the fitted probabilities, iterated to
#' convergence. The coefficients are Firth's penalized maximum-likelihood
#' estimates, finite under complete separation. \code{vcov()} and \code{sigma()}
#' are those of the quasibinomial fit to the adjusted counts, so the dispersion
#' is estimated from the scatter of the counts between samples.
#'
#' @return An R6 class generator.
#' @export
#' @examples
#' dat <- data.frame(
#'   group_ = factor(rep(c("A", "B"), each = 4)),
#'   detected = c(0, 0, 0, 0, 4, 5, 3, 5),
#'   undetected = c(5, 5, 5, 5, 1, 0, 2, 0)
#' )
#' strategy <- StrategyBinomialFirth$new("~ group_")
#' fit <- strategy$model_fun(dat)
#' coefficients(fit)
StrategyBinomialFirth <- R6::R6Class(
  "StrategyBinomialFirth",
  inherit = StrategyBinomial,
  public = list(
    #' @description Create a Firth-penalized quasibinomial count strategy.
    #' @param modelstr right-hand-side model formula, for example
    #'   \code{"~ group_"}
    #' @param model_name model identity
    #' @param control named list of convergence controls, see \code{\link{strategy_binomial_firth}}
    initialize = function(modelstr, model_name = "binomial_firth_nested", control = list()) {
      super$initialize(modelstr, prior_count = 0, model_name = model_name, control = control)
    }
  ),
  private = list(
    control_defaults = function() list(maxit = 100, tol = 1e-8),
    # the inner glm fits keep glm's defaults; the controls apply to the Firth iteration
    glm_control = function() list(),
    fit = function(x) {
      trials <- x$detected + x$undetected
      hat <- rep(1, nrow(x))
      previous <- NULL
      for (iteration in seq_len(self$control$maxit)) {
        x$.detected <- x$detected + hat / 2
        x$.undetected <- x$undetected + hat / 2
        model <- super$fit(x)
        if (!is.null(previous) && max(abs(stats::coef(model) - previous)) < self$control$tol) {
          model$firth_iterations <- iteration
          return(model)
        }
        previous <- stats::coef(model)
        design <- stats::model.matrix(model)
        weight <- trials * stats::fitted(model) * (1 - stats::fitted(model))
        hat <- weight * rowSums((design %*% solve(crossprod(design, weight * design))) * design)
      }
      abort_prolfqua("The Firth iteration did not converge.")
    }
  )
)


#' Create a Firth-penalized quasibinomial detection-count strategy
#'
#' @param modelstr right-hand-side model formula
#' @param model_name model identity
#' @param control \code{maxit} (default 100) and \code{tol} (default 1e-8): iteration limit of the Firth
#'   adjustment and the largest coefficient change (log odds) that stops it; log odds of detection range over
#'   about -10 to 10, so 1e-8 is far below any reported precision.
#' @return A \code{\link{StrategyBinomialFirth}} object.
#' @export
#' @examples
#' strategy <- strategy_binomial_firth("~ group_")
#' strategy$formula
strategy_binomial_firth <- function(modelstr, model_name = "binomial_firth_nested", control = list()) {
  StrategyBinomialFirth$new(modelstr, model_name, control = control)
}
