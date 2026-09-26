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
StrategyBase <- R6::R6Class(
  "StrategyBase",
  public = list(
    #' @field formula model formula
    formula = NULL,
    #' @field model_name name of model
    model_name = NULL,
    #' @field anova_df ANOVA extractor, see \code{\link{AnovaExtractor}}
    anova_df = NULL,

    #' @description Create a new strategy
    #' @param modelstr model formula string
    #' @param model_name name of model
    initialize = function(modelstr, model_name = "Model") {
      self$formula <- as.formula(modelstr)
      self$model_name <- model_name
      self$anova_df <- get_anova_df(test = "F")
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
    #' @param ... passed to \code{\link{compute_contrast}}
    contrast_fun = function(...) compute_contrast(...),

    #' @description Get residual degrees of freedom
    #' @param model fitted model
    df_residual = function(model) stats::df.residual(model),

    #' @description Get residual standard error
    #' @param model fitted model
    sigma = function(model) stats::sigma(model)
  ),
  private = list(
    # Input preparation that runs outside the error handler of model_fun().
    prepare = function(x) x,
    fit = function(x) stop("fit is not implemented")
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
#' @examples
#' istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 10, with_missing = FALSE)
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
    #' @param ... passed to \code{\link{compute_lmer_contrast}}
    contrast_fun = function(...) compute_lmer_contrast(...)
  ),
  private = list(
    fit = function(x) lmerTest::lmer(self$formula, data = x)
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
strategy_lmer <- function(modelstr, model_name = "Model") {
  StrategyLmer$new(modelstr, model_name)
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
    fit = function(x) MASS::rlm(self$formula, data = x, method = "M")
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
strategy_rlm <- function(modelstr, model_name = "Model") {
  StrategyRLM$new(modelstr, model_name)
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
    #' @field test statistical test type ("F", "Chisq", etc.)
    test = "F",
    #' @field col_pval p-value column name in ANOVA output
    col_pval = character(),
    #' @description Create a new AnovaExtractor
    #' @param test statistical test type
    initialize = function(test = "F") {
      self$test <- test
      self$col_pval <- paste0("Pr..", substr(test, 1, 3), ".")
    },
    #' @description Extract ANOVA table from a fitted model
    #' @param x a fitted model
    #' @return data.frame with factor column and ANOVA statistics
    fun = function(x) {
      x <- anova(x, test = self$test)
      colnames(x) <- make.names(colnames(x))
      x <- data.frame(factor = rownames(x), x)
      return(x)
    }
  )
)

#' anova returning dataframe
#' @keywords internal
#' @family modelling
#' @export
#' @examples
#' x <- get_anova_df(test = "F")
#' x <- get_anova_df(test = "Chisq")
get_anova_df <- function(test = "F") {
  AnovaExtractor$new(test = test)
}
