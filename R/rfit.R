# Rfit rank-based regression backend ----
#
# `Rfit::rfit()` fits rank-based linear models with an `lm`-like formula
# interface. The fitted object exposes `coef()` and `vcov()` but, unlike `lm`,
# carries no model frame or terms and returns an unnamed covariance matrix.
# `linfct_from_model()` and `compute_contrast()` rely on all three. The strategy
# therefore augments each fit with its model frame and terms and tags it with
# the `rfit_prolfqua` subclass, for which the S3 methods below restore named
# `vcov()` and supply the rank scale / residual df. With that glue the backend
# reuses the classic `build_model()` / `compute_contrast()` Wald path unchanged.

#' Named variance-covariance matrix for an augmented rfit fit
#'
#' `Rfit::vcov.rfit` returns an unnamed matrix; `compute_contrast()` indexes the
#' covariance by coefficient name, so this method restores the dimnames from
#' `coef()`.
#' @param object an `rfit_prolfqua` fit
#' @param ... passed to the next method
#' @return the rfit covariance matrix with row/column names set to the
#'   coefficient names.
#' @exportS3Method stats::vcov
#' @family modelling
#' @examples
#' fit <- strategy_rfit("Sepal.Length ~ Species")$model_fun(iris)
#' vcov(fit)
vcov.rfit_prolfqua <- function(object, ...) {
  v <- NextMethod()
  nm <- names(stats::coef(object))
  dimnames(v) <- list(nm, nm)
  v
}


#' Residual degrees of freedom for an augmented rfit fit
#'
#' @param object an `rfit_prolfqua` fit
#' @param ... ignored
#' @return number of observations minus number of coefficients.
#' @exportS3Method stats::df.residual
#' @family modelling
#' @examples
#' fit <- strategy_rfit("Sepal.Length ~ Species")$model_fun(iris)
#' df.residual(fit)
df.residual.rfit_prolfqua <- function(object, ...) {
  length(object$residuals) - length(stats::coef(object))
}


#' Scale estimate for an augmented rfit fit
#'
#' Returns the rank-based scale `tauhat`. This stands in for the residual
#' standard deviation of an `lm`; it does not affect the Wald std.error or
#' p-value (those come from `vcov()`), but it is reported and used as the per-
#' protein variance estimate if results are passed to limma/eBayes moderation.
#' @param object an `rfit_prolfqua` fit
#' @param ... ignored
#' @return the rank-based scale estimate `tauhat`.
#' @exportS3Method stats::sigma
#' @family modelling
#' @examples
#' fit <- strategy_rfit("Sepal.Length ~ Species")$model_fun(iris)
#' sigma(fit)
sigma.rfit_prolfqua <- function(object, ...) {
  object$tauhat
}


#' Rank-based regression strategy (R6 class)
#'
#' Encapsulates everything needed to fit per-protein rank-based linear models
#' via \code{\link[Rfit]{rfit}} and extract Wald contrasts. The fit is augmented
#' (model frame, terms, \code{rfit_prolfqua} subclass) so it satisfies the same
#' \code{coef()} / \code{vcov()} / \code{terms()} contract as \code{lm}.
#'
#' Unlike \code{\link{StrategyLM}}, \code{rfit} takes no observation weights, so
#' \code{nr_children} weighting is not supported by this backend.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' strat <- StrategyRfit$new("Intensity ~ condition", model_name = "parallel design")
#' strat$model_fun(get_formula = TRUE)
StrategyRfit <- R6::R6Class(
  "StrategyRfit",
  public = list(
    #' @field formula model formula
    formula = NULL,
    #' @field model_name name of model
    model_name = NULL,
    #' @field report_columns columns to report
    report_columns = NULL,
    #' @field is_mixed always FALSE for rfit
    is_mixed = FALSE,
    #' @field anova_df list with anova function and column names
    anova_df = NULL,

    #' @description Create a new StrategyRfit
    #' @param modelstr model formula string
    #' @param model_name name of model
    #' @param report_columns columns to report
    initialize = function(
      modelstr,
      model_name = "rfit",
      report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value", "moderated.p.value.adjusted")
    ) {
      self$formula <- as.formula(modelstr)
      self$model_name <- model_name
      self$report_columns <- report_columns
      self$anova_df <- get_anova_df(test = "F")
    },

    #' @description Fit rfit to one protein's data, augmenting the fit so the
    #'   classic contrast path can introspect it.
    #' @param x data.frame for one protein
    #' @param pb optional progress bar
    #' @param get_formula if TRUE, return formula instead of fitting
    model_fun = function(x, pb, get_formula = FALSE) {
      if (get_formula) {
        return(self$formula)
      }
      if (!missing(pb)) {
        pb$tick()
      }
      if (!requireNamespace("Rfit", quietly = TRUE)) {
        stop("Package 'Rfit' is required for the rfit backend. Install it with install.packages('Rfit').")
      }
      tryCatch(
        {
          fit <- Rfit::rfit(self$formula, data = x)
          # Rfit silently zeros unestimable coefficients in rank-deficient
          # designs (unlike lm which uses NA), and the resulting vcov()
          # then fails with a Cholesky error. Detect and fail the fit so
          # the protein is reported in get_missing() rather than crashing
          # the contrast loop.
          if (!is.null(fit$qrx1$rank) && fit$qrx1$rank < ncol(fit$x)) {
            return(.error_handler(simpleError("rfit design is rank-deficient")))
          }
          fit$model <- stats::model.frame(self$formula, data = x)
          fit$terms <- stats::terms(self$formula, data = x)
          class(fit) <- c("rfit_prolfqua", class(fit))
          fit
        },
        error = .error_handler
      )
    },

    #' @description Check if model is singular (NA coefficients or df < 2)
    #' @param model fitted model
    isSingular = function(model) {
      if (any(is.na(coefficients(model)))) {
        return(TRUE)
      }
      df <- self$df_residual(model)
      if (is.na(df) || df < 2) {
        return(TRUE)
      }
      FALSE
    },

    #' @description Compute contrasts from fitted model
    #' @param ... passed to \code{\link{compute_contrast}}
    contrast_fun = function(...) compute_contrast(...),

    #' @description Get residual degrees of freedom
    #' @param model fitted model
    df_residual = function(model) length(model$residuals) - length(stats::coef(model)),

    #' @description Get the rank-based scale estimate
    #' @param model fitted model
    sigma = function(model) model$tauhat
  )
)


#' Create rank-based regression strategy
#'
#' Convenience wrapper that creates a \code{\link{StrategyRfit}} object.
#' @rdname strategy
#' @export
#' @param modelstr model formula
#' @param model_name name of model
#' @param report_columns columns to report
#' @family modelling
#' @return a \code{\link{StrategyRfit}} object
#' @examples
#' tmp <- strategy_rfit("Intensity ~ condition", model_name = "parallel design")
#' tmp$model_fun(get_formula = TRUE)
strategy_rfit <- function(
  modelstr,
  model_name = "rfit",
  report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value", "moderated.p.value.adjusted")
) {
  StrategyRfit$new(modelstr, model_name, report_columns)
}
