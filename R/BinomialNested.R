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
    initialize = function(modelstr, prior_count = 0.1, model_name = "binomial_nested") {
      prior_count <- .validate_binomial_prior_count(prior_count)
      super$initialize(paste("cbind(.detected, .undetected)", modelstr), model_name)
      self$prior_count <- prior_count
    }
  ),
  private = list(
    prepare = function(x) {
      x$.detected <- x$detected + self$prior_count
      x$.undetected <- x$undetected + self$prior_count
      x
    },
    fit = function(x) {
      model <- stats::glm(self$formula, data = x, family = stats::quasibinomial())
      if (anyNA(stats::coef(model)) || model$rank < length(stats::coef(model))) {
        stop("binomial design is rank-deficient")
      }
      if (stats::df.residual(model) < 2L) {
        stop("binomial model requires at least two residual degrees of freedom")
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
#' @return A \code{\link{StrategyBinomial}} object.
#' @export
#' @examples
#' strategy <- strategy_binomial("~ group_", prior_count = 0.1)
#' strategy$formula
strategy_binomial <- function(modelstr, prior_count = 0.1, model_name = "binomial_nested") {
  StrategyBinomial$new(modelstr, prior_count, model_name)
}
