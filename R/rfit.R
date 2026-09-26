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
#' strat$formula
StrategyRfit <- R6::R6Class(
  "StrategyRfit",
  inherit = StrategyBase,
  public = list(
    #' @description Create a new StrategyRfit. Residual df and scale come from
    #'   \code{df.residual()} and \code{sigma()} of the augmented fit.
    #' @param modelstr model formula string
    #' @param model_name name of model
    initialize = function(modelstr, model_name = "rfit") {
      super$initialize(modelstr, model_name)
    }
  ),
  private = list(
    prepare = function(x) {
      if (!requireNamespace("Rfit", quietly = TRUE)) {
        stop("Package 'Rfit' is required for the rfit backend. Install it with install.packages('Rfit').")
      }
      x
    },
    # Fit rfit, augmenting the fit so the classic contrast path can introspect it.
    fit = function(x) {
      fit <- Rfit::rfit(self$formula, data = x)
      # When the rank fit is no better than the intercept-only model,
      # `Rfit::rfit()` falls back to `bhat0 <- c(median(y), rep(0, p))`,
      # an UNNAMED coefficient vector at full QR rank. The names matter:
      # `linfct_from_model()` -> `.model_coeff_matrix()` copies them onto
      # the coefficient matrix columns, and the contrast path then crashes
      # in `tibble::as_tibble()` ("Columns must be named") when they are
      # NULL. The design matrix always carries the correct column names,
      # so restore them. (This is distinct from the rank-deficient case
      # below, which keeps full rank and so is not caught by that guard.)
      if (is.null(names(fit$coefficients))) {
        names(fit$coefficients) <- colnames(fit$x)
      }
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
    }
  )
)


#' Create rank-based regression strategy
#'
#' Convenience wrapper that creates a \code{\link{StrategyRfit}} object.
#' @rdname strategy
#' @export
#' @family modelling
#' @return a \code{\link{StrategyRfit}} object
#' @examples
#' tmp <- strategy_rfit("Intensity ~ condition", model_name = "parallel design")
#' tmp$formula
strategy_rfit <- function(modelstr, model_name = "rfit") {
  StrategyRfit$new(modelstr, model_name)
}
