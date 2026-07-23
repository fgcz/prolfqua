# Binomial detection model ----

.summarize_detection_counts <- function(lfqdata) {
  stopifnot("LFQData" %in% class(lfqdata))
  bin_resp <- lfqdata$get_config()$bin_resp
  grouping <- unique(c(
    lfqdata$subject_id(),
    lfqdata$sample_name(),
    lfqdata$isotope_label(),
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

  stopifnot(all(counts$detected >= 0L))
  stopifnot(all(counts$undetected >= 0L))
  counts
}

.validate_binomial_prior_count <- function(prior_count) {
  candidate <- if (is.numeric(prior_count)) prior_count else NA_real_
  is_valid <- all(c(
    length(candidate) == 1L,
    !anyNA(candidate),
    is.finite(candidate),
    candidate >= 0
  ))
  if (!is_valid) {
    stop("`prior_count` must be one non-negative number.", call. = FALSE)
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
StrategyBinomial <- R6::R6Class(
  "StrategyBinomial",
  public = list(
    #' @field formula quasibinomial model formula
    formula = NULL,
    #' @field model_name model identity
    model_name = NULL,
    #' @field report_columns result columns supported by the strategy
    report_columns = NULL,
    #' @field is_mixed always FALSE
    is_mixed = FALSE,
    #' @field anova_df ANOVA extractor
    anova_df = NULL,
    #' @field prior_count symmetric pseudo-count added to both outcomes
    prior_count = NULL,

    #' @description Create a quasibinomial count strategy.
    #' @param modelstr right-hand-side model formula, for example
    #'   \code{"~ group_"}
    #' @param prior_count non-negative symmetric pseudo-count
    #' @param model_name model identity
    #' @param report_columns result columns supported by the strategy
    initialize = function(
      modelstr,
      prior_count = 0.1,
      model_name = "binomial_nested",
      report_columns = c(
        "statistic",
        "p.value",
        "p.value.adjusted",
        "moderated.p.value",
        "moderated.p.value.adjusted"
      )
    ) {
      prior_count <- .validate_binomial_prior_count(prior_count)
      self$formula <- stats::as.formula(paste("cbind(.detected, .undetected)", modelstr))
      self$model_name <- model_name
      self$report_columns <- report_columns
      self$anova_df <- get_anova_df(test = "F")
      self$prior_count <- prior_count
    },

    #' @description Fit one parent's detection counts.
    #' @param x parent-by-sample count data
    #' @param pb optional progress reporter
    #' @param get_formula if TRUE, return the model formula without fitting
    model_fun = function(x, pb, get_formula = FALSE) {
      if (get_formula) {
        return(self$formula)
      }
      if (!missing(pb)) {
        pb$tick()
      }

      x$.detected <- x$detected + self$prior_count
      x$.undetected <- x$undetected + self$prior_count
      tryCatch(
        {
          model <- stats::glm(self$formula, data = x, family = stats::quasibinomial())
          if (anyNA(stats::coef(model)) || model$rank < length(stats::coef(model))) {
            stop("binomial design is rank-deficient")
          }
          if (stats::df.residual(model) < 2L) {
            stop("binomial model requires at least two residual degrees of freedom")
          }
          model
        },
        error = .error_handler
      )
    },

    #' @description Check whether the model is singular.
    #' @param model fitted quasibinomial model
    isSingular = function(model) is_singular_lm(model),

    #' @description Compute linear contrasts.
    #' @param ... passed to \code{\link{compute_contrast}}
    contrast_fun = function(...) compute_contrast(...),

    #' @description Return residual degrees of freedom.
    #' @param model fitted quasibinomial model
    df_residual = function(model) stats::df.residual(model),

    #' @description Return the Pearson residual scale used by \code{vcov()}.
    #' @param model fitted quasibinomial model
    sigma = function(model) stats::sigma(model)
  )
)


#' Create a quasibinomial detection-count strategy
#'
#' @param modelstr right-hand-side model formula
#' @param prior_count non-negative symmetric pseudo-count
#' @param model_name model identity
#' @param report_columns result columns supported by the strategy
#' @return A \code{\link{StrategyBinomial}} object.
#' @export
#' @examples
#' strategy <- strategy_binomial("~ group_", prior_count = 0.1)
#' strategy$model_fun(get_formula = TRUE)
strategy_binomial <- function(
  modelstr,
  prior_count = 0.1,
  model_name = "binomial_nested",
  report_columns = c(
    "statistic",
    "p.value",
    "p.value.adjusted",
    "moderated.p.value",
    "moderated.p.value.adjusted"
  )
) {
  StrategyBinomial$new(modelstr, prior_count, model_name, report_columns)
}
