# Creating models from configuration ----

#' Residual degrees of freedom for rlm objects
#'
#' \code{stats::df.residual} returns \code{NA} for \code{\link[MASS]{rlm}}
#' objects. This S3 method computes weighted residual df instead.
#'
#' @param object an \code{rlm} object
#' @param ... ignored
#' @return numeric scalar
#' @keywords internal
#' @exportS3Method
df.residual.rlm <- function(object, ...) {
  sum(object$w) - object$rank
}

#' Residual scale estimate for rlm objects
#'
#' \code{stats::sigma} returns \code{NA} for \code{\link[MASS]{rlm}} objects.
#' This S3 method computes a weighted residual scale estimate instead.
#'
#' @param object an \code{rlm} object
#' @param ... ignored
#' @return numeric scalar
#' @keywords internal
#' @exportS3Method
sigma.rlm <- function(object, ...) {
  sqrt(sum(object$w * object$resid^2) / (sum(object$w) - object$rank))
}

#' Base class for per-subject model strategies
#'
#' Holds the fields, constructor and default methods shared by \code{\link{StrategyLM}},
#' \code{\link{StrategyRLM}}, \code{\link{StrategyLmer}}, \code{\link{StrategyRfit}},
#' \code{\link{StrategyLogistf}} and \code{\link{StrategyBinomial}}. A subclass implements the private
#' \code{fit(x)} method and overrides only the methods that differ.
#'
#' @return An R6 class generator.
#' @keywords internal
#' @srrstats {RE3.0} Failed fits and warnings raised while fitting or aggregating, including non-convergence of rlm,
#'   lmer, Firth and medpolish, give one `prolfqua_fit_failed` warning per step naming the subjects and the reason.
#' @srrstats {RE4.1} Strategies expose a formula and fitting specification without fitting a model.
StrategyBase <- R6::R6Class(
  "StrategyBase",
  public = list(
    #' @field formula model formula
    formula = NULL,
    #' @field model_name name of model
    model_name = NULL,
    #' @field anova_df ANOVA extractor, see \code{\link{AnovaExtractor}}
    anova_df = NULL,
    #' @field control convergence controls passed to the fitting backend (empty for non-iterative fits)
    control = NULL,

    #' @description Create a new strategy
    #' @param modelstr model formula string
    #' @param model_name name of model
    #' @param control named list of convergence controls merged over the strategy's defaults; names the
    #'   backend does not accept raise a classed error
    initialize = function(modelstr, model_name = "Model", control = list()) {
      self$formula <- .as_model_formula(modelstr)
      self$model_name <- model_name
      self$anova_df <- .get_anova_df(test = "F")
      allowed <- if (length(control) > 0L) private$control_names() else names(private$control_defaults())
      self$control <- .merge_control(control, private$control_defaults(), allowed)
    },

    #' @description Fit the model to one subject's data. A failed fit returns the error message.
    #' @param x data.frame for one subject
    #' @param pb optional progress bar
    model_fun = function(x, pb) {
      if (!missing(pb)) {
        pb$tick()
      }
      x <- private$prepare(x)
      tryCatch(private$fit(x), error = .error_handler)
    },

    #' @description Check if model is singular (NA coefficients or df < 2)
    #' @param model fitted model
    isSingular = function(model) {
      df <- self$df_residual(model)
      anyNA(coefficients(model)) || is.na(df) || df < 2
    },

    #' @description Compute contrasts from fitted model
    #' @param ... \code{m} (fitted model), \code{linfct} (contrast weights on its coefficients) and
    #'   \code{confint}, the confidence level (default 0.95)
    contrast_fun = function(...) .compute_estimable_contrast(...),

    #' @description Get residual degrees of freedom
    #' @param model fitted model
    df_residual = function(model) stats::df.residual(model),

    #' @description Get residual standard error
    #' @param model fitted model
    sigma = function(model) stats::sigma(model),

    #' @description Flag residual scales indistinguishable from zero (constant or noiseless response); their
    #'   tests are set to NA
    #' @param sigma residual scales, one per contrast row
    #' @param scale response level of the same rows (average abundance)
    zero_residual_variance = function(sigma, scale) .zero_residual_variance(sigma, scale)
  ),
  private = list(
    # convergence controls of the backend and their defaults; non-iterative fits have none
    control_defaults = function() list(),
    control_names = function() names(private$control_defaults()),
    # Input preparation that runs outside the error handler of model_fun().
    prepare = function(x) x,
    fit = function(x) abort_prolfqua("`fit()` is not implemented for this class.")
  )
)

#' Linear mixed-effects model strategy (R6 class)
#'
#' Encapsulates everything needed to fit per-protein linear mixed-effects models
#' via \code{\link[lmerTest]{lmer}} and extract contrasts.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examplesIf requireNamespace("lme4", quietly = TRUE) && requireNamespace("lmerTest", quietly = TRUE)
#' istar <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 10, with_missing = FALSE)
#' istar <- prolfqua::LFQData$new(istar$data, istar$config)
#' istar$set_data(istar$data_long() |> dplyr::group_by(protein_Id) |>
#'   dplyr::mutate(abundanceC = abundance - mean(abundance)) |> dplyr::ungroup())
#' strat <- StrategyLmer$new("abundanceC ~ group_ + (1|peptide_Id)",
#'   model_name = "random_example")
#' mod <- build_model(istar, strat)
#' sum(mod$model_df$has_model_fit)
StrategyLmer <- R6::R6Class(
  "StrategyLmer",
  inherit = StrategyBase,
  public = list(
    #' @description Check if model is singular
    #' @param model fitted model
    isSingular = function(model) lme4::isSingular(model),

    #' @description Compute contrasts from fitted model
    #' @param ... \code{model} (fitted model), \code{linfct} (contrast weights on its coefficients) and
    #'   \code{ddf}, the denominator degrees of freedom: \code{"Satterthwaite"} (default) or \code{"Kenward-Roger"}
    contrast_fun = function(...) .compute_lmer_contrast(...)
  ),
  private = list(
    prepare = function(x) {
      rlang::check_installed(c("lme4", "lmerTest"), reason = "for mixed models.")
      x
    },
    control_names = function() {
      rlang::check_installed("lme4", reason = "for mixed models.")
      names(formals(lme4::lmerControl))
    },
    fit = function(x) lmerTest::lmer(self$formula, data = x, control = do.call(lme4::lmerControl, self$control))
  )
)


#' Create linear mixed-effects model strategy
#'
#' Convenience wrapper that creates a \code{\link{StrategyLmer}} object.
#' @rdname strategy
#' @export
#' @family modelling
#' @return a \code{\link{StrategyLmer}} object
#' @examples
#' modelFunction <- strategy_lmer("abundanceC ~ group_ + (1|peptide_Id)",
#'   model_name = "random_example")
#' modelFunction$formula
#' @param control for \code{strategy_lmer()}: arguments of \code{\link[lme4]{lmerControl}}, e.g.
#'   \code{list(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))}; default: lme4's defaults (nloptwrap,
#'   gradient tolerance 0.002), which suit per-protein models of a few dozen observations.
strategy_lmer <- function(modelstr, model_name = "Model", control = list()) {
  StrategyLmer$new(modelstr, model_name, control = control)
}

#' Linear model strategy (R6 class)
#'
#' Encapsulates everything needed to fit per-protein linear models and extract
#' contrasts: the formula, model fitting function, singularity check, contrast
#' computation, ANOVA, and residual statistics.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' strat <- StrategyLM$new("Intensity ~ condition", model_name = "parallel design")
#' strat$formula
#' strat$weights
StrategyLM <- R6::R6Class(
  "StrategyLM",
  inherit = StrategyBase,
  public = list(
    #' @field weights optional character string naming a column in the data
    #'   containing per-observation weights, passed to \code{\link[stats]{lm}}.
    weights = NULL,

    #' @description Create a new StrategyLM
    #' @param modelstr model formula string
    #' @param model_name name of model
    #' @param weights optional character string naming a column in the data
    #'   containing per-observation weights
    initialize = function(modelstr, model_name = "Model", weights = NULL) {
      super$initialize(modelstr, model_name)
      self$weights <- weights
    }
  ),
  private = list(
    fit = function(x) {
      if (is.null(self$weights)) {
        return(lm(self$formula, data = x))
      }
      eval(bquote(lm(.(self$formula), data = x, weights = .(as.name(self$weights)))))
    }
  )
)


#' Create linear model strategy
#'
#' Convenience wrapper that creates a \code{\link{StrategyLM}} object.
#' @rdname strategy
#' @export
#' @param modelstr model formula
#' @param model_name name of model
#' @param weights optional character string naming a column in the data
#'   containing per-observation weights, passed to \code{\link[stats]{lm}}.
#'   Default \code{NULL} (unweighted).
#' @family modelling
#' @return a \code{\link{StrategyLM}} object
#' @examples
#' tmp <- strategy_lm("Intensity ~ condition", model_name = "parallel design")
#' tmp$formula
#' tmp$weights
strategy_lm <- function(modelstr, model_name = "Model", weights = NULL) {
  StrategyLM$new(modelstr, model_name, weights)
}


#' Robust linear model strategy (R6 class)
#'
#' Encapsulates everything needed to fit per-protein robust linear models
#' via \code{\link[MASS]{rlm}} and extract contrasts.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' strat <- StrategyRLM$new("Intensity ~ condition", model_name = "parallel design")
#' strat$formula
StrategyRLM <- R6::R6Class(
  "StrategyRLM",
  inherit = StrategyBase,
  public = list(
    #' @description Get the robust scale estimate used by \code{vcov()}
    #'
    #' Returns \code{model$s}, the robust scale \code{MASS::rlm} builds its
    #' \code{vcov()} / \code{summary()} standard errors from. This must match the
    #' scale embedded in \code{vcov()} so that \code{ContrastsModerated}'s
    #' \code{sigma / sqrt(var.post)} rescaling stays coherent; \code{stats::sigma()}
    #' returns the ordinary-residual scale instead and would distort moderated
    #' \code{rlm} statistics.
    #' @param model fitted model
    sigma = function(model) model$s
  ),
  private = list(
    control_defaults = function() list(maxit = 20, acc = 1e-4),
    fit = function(x) {
      ls <- stats::lm(self$formula, data = x)
      response <- stats::model.response(stats::model.frame(ls))
      if (!.zero_residual_variance(stats::sigma(ls), mean(response))) {
        return(MASS::rlm(self$formula, data = x, method = "M", maxit = self$control$maxit, acc = self$control$acc))
      }
      # A constant or noiseless response has no robust scale: rlm's IRLS divides residuals of floating-point
      # noise by a scale that can underflow to zero, and on some platforms all weights and coefficients turn
      # NA. Fit rlm to an all-zero response (exact zero scale) and carry over the least-squares estimates.
      zero <- x
      zero[[deparse(self$formula[[2]])]] <- 0
      fit <- MASS::rlm(self$formula, data = zero, method = "M", maxit = self$control$maxit, acc = self$control$acc)
      fit$coefficients <- stats::coef(ls)
      fit$fitted.values <- stats::fitted(ls)
      fit$residuals <- stats::residuals(ls)
      fit$model <- ls$model
      fit
    }
  )
)


#' Create robust linear model strategy
#'
#' Convenience wrapper that creates a \code{\link{StrategyRLM}} object.
#' @rdname strategy
#' @export
#' @family modelling
#' @return a \code{\link{StrategyRLM}} object
#' @examples
#' tmp <- strategy_rlm("Intensity ~ condition", model_name = "parallel design")
#' tmp$formula
#' strategy_rlm("Intensity ~ condition", control = list(maxit = 50))$control
#' @param control for \code{strategy_rlm()}: \code{maxit} (default 20) and \code{acc} (default 1e-4), the
#'   iteration limit and the relative change of the IRLS residuals that stops \code{\link[MASS]{rlm}}; the
#'   criterion is relative, so it does not depend on the intensity scale, and 1e-4 is far below the precision of
#'   log2 intensities.
strategy_rlm <- function(modelstr, model_name = "Model", control = list()) {
  StrategyRLM$new(modelstr, model_name, control = control)
}


#' R6 class for extracting ANOVA results as a data frame
#'
#' Wraps the ANOVA extraction logic and associated column names.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' extractor <- AnovaExtractor$new(test = "F")
#' fit <- stats::lm(Sepal.Length ~ Species, data = iris)
#' res <- extractor$fun(fit)
#' stopifnot("factor" %in% colnames(res))
AnovaExtractor <- R6::R6Class(
  "AnovaExtractor",
  public = list(
    #' @field test one case-sensitive test choice: "F", "Chisq", "LRT", "Rao" or "Cp".
    test = "F",
    #' @field col_pval p-value column name in ANOVA output; empty for "Cp", which has no p-value.
    col_pval = character(),
    #' @description Create a new AnovaExtractor
    #' @param test one case-sensitive test choice: "F", "Chisq", "LRT", "Rao" or "Cp".
    #'   The fitted model's ANOVA method must support the chosen test.
    initialize = function(test = "F") {
      self$test <- .arg_match_prolfqua(test, c("F", "Chisq", "LRT", "Rao", "Cp"), "test")
      self$col_pval <- switch(self$test, F = "Pr..F.", Cp = character(), "Pr..Chi.")
    },
    #' @description Extract ANOVA table from a fitted model
    #' @param x a fitted model
    #' @return data.frame with factor column and ANOVA statistics
    fun = function(x) {
      x <- anova(x, test = self$test)
      colnames(x) <- make.names(colnames(x))
      x <- data.frame(factor = rownames(x), x)
      x
    }
  )
)

#' Create an ANOVA extractor returning a data frame
#' @param test one case-sensitive test choice: "F", "Chisq", "LRT", "Rao" or "Cp".
#'   The fitted model's ANOVA method must support the chosen test.
#' @return An \code{\link{AnovaExtractor}} object.
#' @noRd
.get_anova_df <- function(test = "F") {
  AnovaExtractor$new(test = test)
}
