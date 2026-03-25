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
#' @export
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
#' @export
sigma.rlm <- function(object, ...) {
  sqrt(sqrt(sum(object$w * object$resid^2) / (sum(object$w) - object$rank)))
}

#' Linear mixed-effects model strategy (R6 class)
#'
#' Encapsulates everything needed to fit per-protein linear mixed-effects models
#' via \code{\link[lmerTest]{lmer}} and extract contrasts.
#'
#' @export
#' @family modelling
#' @examples
#' istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 10, with_missing = FALSE)
#' istar <- prolfqua::LFQData$new(istar$data, istar$config)
#' istar$data <- istar$data |> dplyr::group_by(protein_Id) |>
#'   dplyr::mutate(abundanceC = abundance - mean(abundance)) |> dplyr::ungroup()
#' strat <- StrategyLmer$new("abundanceC ~ group_ + (1|peptide_Id)",
#'   model_name = "random_example")
#' mod <- build_model(istar, strat)
#' sum(mod$modelDF$has_model_fit)
StrategyLmer <- R6::R6Class(
  "StrategyLmer",
  public = list(
    #' @field formula model formula
    formula = NULL,
    #' @field model_name name of model
    model_name = NULL,
    #' @field report_columns columns to report
    report_columns = NULL,
    #' @field is_mixed always TRUE for lmer
    is_mixed = TRUE,
    #' @field anova_df list with anova function and column names
    anova_df = NULL,

    #' @description Create a new StrategyLmer
    #' @param modelstr model formula string
    #' @param model_name name of model
    #' @param report_columns columns to report
    initialize = function(
      modelstr,
      model_name = "Model",
      report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value", "moderated.p.value.adjusted")
    ) {
      self$formula <- as.formula(modelstr)
      self$model_name <- model_name
      self$report_columns <- report_columns
      self$anova_df <- get_anova_df(test = "F")
    },

    #' @description Fit lmer to one protein's data
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
      tryCatch(lmerTest::lmer(self$formula, data = x), error = .error_handler)
    },

    #' @description Check if model is singular
    #' @param model fitted model
    isSingular = function(model) lme4::isSingular(model),

    #' @description Compute contrasts from fitted model
    #' @param ... passed to \code{\link{compute_lmer_contrast}}
    contrast_fun = function(...) compute_lmer_contrast(...),

    #' @description Get residual degrees of freedom
    #' @param model fitted model
    df_residual = function(model) stats::df.residual(model),

    #' @description Get residual standard error
    #' @param model fitted model
    sigma = function(model) stats::sigma(model)
  )
)


#' Create linear mixed-effects model strategy
#'
#' Convenience wrapper that creates a \code{\link{StrategyLmer}} object.
#' @rdname strategy
#' @export
#' @param modelstr model formula
#' @param model_name name of model
#' @param report_columns columns to report
#' @family modelling
#' @return a \code{\link{StrategyLmer}} object
#' @examples
#' modelFunction <- strategy_lmer("abundanceC ~ group_ + (1|peptide_Id)",
#'   model_name = "random_example")
#' modelFunction$model_fun(get_formula = TRUE)
strategy_lmer <- function(
  modelstr,
  model_name = "Model",
  report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value", "moderated.p.value.adjusted")
) {
  StrategyLmer$new(modelstr, model_name, report_columns)
}

#' Linear model strategy (R6 class)
#'
#' Encapsulates everything needed to fit per-protein linear models and extract
#' contrasts: the formula, model fitting function, singularity check, contrast
#' computation, ANOVA, and residual statistics.
#'
#' @export
#' @family modelling
#' @examples
#' strat <- StrategyLM$new("Intensity ~ condition", model_name = "parallel design")
#' strat$model_fun(get_formula = TRUE)
#' strat$weights
StrategyLM <- R6::R6Class(
  "StrategyLM",
  public = list(
    #' @field formula model formula
    formula = NULL,
    #' @field model_name name of model
    model_name = NULL,
    #' @field report_columns columns to report
    report_columns = NULL,
    #' @field weights optional character string naming a column in the data
    #'   containing per-observation weights, passed to \code{\link[stats]{lm}}.
    weights = NULL,
    #' @field is_mixed always FALSE for lm
    is_mixed = FALSE,
    #' @field anova_df list with anova function and column names
    anova_df = NULL,

    #' @description Create a new StrategyLM
    #' @param modelstr model formula string
    #' @param model_name name of model
    #' @param report_columns columns to report
    #' @param weights optional character string naming a column in the data
    #'   containing per-observation weights
    initialize = function(
      modelstr,
      model_name = "Model",
      report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value", "moderated.p.value.adjusted"),
      weights = NULL
    ) {
      self$formula <- as.formula(modelstr)
      self$model_name <- model_name
      self$report_columns <- report_columns
      self$weights <- weights
      self$anova_df <- get_anova_df(test = "F")
    },

    #' @description Fit lm to one protein's data
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
      if (!is.null(self$weights)) {
        call <- bquote(lm(.(self$formula), data = x, weights = .(as.name(self$weights))))
        tryCatch(eval(call), error = .error_handler)
      } else {
        tryCatch(lm(self$formula, data = x), error = .error_handler)
      }
    },

    #' @description Check if model is singular
    #' @param model fitted model
    isSingular = function(model) isSingular_lm(model),

    #' @description Compute contrasts from fitted model
    #' @param ... passed to \code{\link{compute_contrast}}
    contrast_fun = function(...) compute_contrast(...),

    #' @description Get residual degrees of freedom
    #' @param model fitted model
    df_residual = function(model) stats::df.residual(model),

    #' @description Get residual standard error
    #' @param model fitted model
    sigma = function(model) stats::sigma(model)
  )
)


#' Create linear model strategy
#'
#' Convenience wrapper that creates a \code{\link{StrategyLM}} object.
#' @rdname strategy
#' @export
#' @param modelstr model formula
#' @param model_name name of model
#' @param report_columns columns to report
#' @param weights optional character string naming a column in the data
#'   containing per-observation weights, passed to \code{\link[stats]{lm}}.
#'   Default \code{NULL} (unweighted).
#' @family modelling
#' @return a \code{\link{StrategyLM}} object
#' @examples
#' tmp <- strategy_lm("Intensity ~ condition", model_name = "parallel design")
#' tmp$model_fun(get_formula = TRUE)
#' tmp$weights
strategy_lm <- function(
  modelstr,
  model_name = "Model",
  report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value", "moderated.p.value.adjusted"),
  weights = NULL
) {
  StrategyLM$new(modelstr, model_name, report_columns, weights)
}


#' Robust linear model strategy (R6 class)
#'
#' Encapsulates everything needed to fit per-protein robust linear models
#' via \code{\link[MASS]{rlm}} and extract contrasts.
#'
#' @export
#' @family modelling
#' @examples
#' strat <- StrategyRLM$new("Intensity ~ condition", model_name = "parallel design")
#' strat$model_fun(get_formula = TRUE)
StrategyRLM <- R6::R6Class(
  "StrategyRLM",
  public = list(
    #' @field formula model formula
    formula = NULL,
    #' @field model_name name of model
    model_name = NULL,
    #' @field report_columns columns to report
    report_columns = NULL,
    #' @field is_mixed always FALSE for rlm
    is_mixed = FALSE,
    #' @field anova_df list with anova function and column names
    anova_df = NULL,

    #' @description Create a new StrategyRLM
    #' @param modelstr model formula string
    #' @param model_name name of model
    #' @param report_columns columns to report
    initialize = function(
      modelstr,
      model_name = "Model",
      report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value", "moderated.p.value.adjusted")
    ) {
      self$formula <- as.formula(modelstr)
      self$model_name <- model_name
      self$report_columns <- report_columns
      self$anova_df <- get_anova_df(test = "F")
    },

    #' @description Fit rlm to one protein's data
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
      tryCatch(MASS::rlm(self$formula, data = x, method = "M"), error = .error_handler)
    },

    #' @description Check if model is singular (NA coefficients or df < 2)
    #' @param model fitted model
    isSingular = function(model) {
      if (any(is.na(coefficients(model)))) {
        return(TRUE)
      }
      df <- df.residual(model)
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
    df_residual = function(model) stats::df.residual(model),

    #' @description Get residual standard error
    #' @param model fitted model
    sigma = function(model) stats::sigma(model)
  )
)


#' Create robust linear model strategy
#'
#' Convenience wrapper that creates a \code{\link{StrategyRLM}} object.
#' @rdname strategy
#' @export
#' @param modelstr model formula
#' @param model_name name of model
#' @param report_columns columns to report
#' @family modelling
#' @return a \code{\link{StrategyRLM}} object
#' @examples
#' tmp <- strategy_rlm("Intensity ~ condition", model_name = "parallel design")
#' tmp$model_fun(get_formula = TRUE)
strategy_rlm <- function(
  modelstr,
  model_name = "Model",
  report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value", "moderated.p.value.adjusted")
) {
  StrategyRLM$new(modelstr, model_name, report_columns)
}


#' R6 class for extracting ANOVA results as a data frame
#'
#' Wraps the ANOVA extraction logic and associated column names.
#'
#' @export
#' @family modelling
AnovaExtractor <- R6::R6Class(
  "AnovaExtractor",
  public = list(
    #' @field test statistical test type ("F", "Chisq", etc.)
    test = "F",
    #' @field col_pval p-value column name in ANOVA output
    col_pval = character(),
    #' @field col_fdr FDR column name
    col_fdr = character(),
    #' @description Create a new AnovaExtractor
    #' @param test statistical test type
    initialize = function(test = "F") {
      self$test <- test
      self$col_pval <- paste0("Pr..", substr(test, 1, 3), ".")
      self$col_fdr <- paste0("FDR.Pr..", substr(test, 1, 3), ".")
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
