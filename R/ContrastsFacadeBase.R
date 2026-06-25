# ContrastsFacadeBase ----

#' Base class for contrast analysis facades
#'
#' Holds the fields and delegating methods shared by every facade in
#' \code{ContrastsFacades.R} and \code{ContrastsChildToParentFacades.R}.
#' Subclasses implement the unique pipeline wiring in \code{initialize()} (build
#' the model, build the inner contrast object, store \code{self$.lfqdata} /
#' \code{self$.contrast_names}) and set \code{self$facade_name} to their registry
#' key. The limma family additionally sets \code{self$.drop_na_diff <- TRUE} to
#' drop rows whose fold change could not be estimated.
#'
#' \code{get_contrasts()} delegates to the inner contrast object and then stamps
#' the facade key into \code{modelName} via \code{.stamp_facade_identity()};
#' rescue/imputation state is carried in the \code{estimate_type} column set by
#' the inner contrast classes.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @keywords internal
ContrastsFacadeBase <- R6::R6Class(
  "ContrastsFacadeBase",
  inherit = ContrastsInterface,
  public = list(
    #' @field model fitted model object
    model = NULL,
    #' @field contrast inner contrast object the facade delegates to
    contrast = NULL,
    #' @field .lfqdata stored reference to the input LFQData (used by get_missing)
    .lfqdata = NULL,
    #' @field .contrast_names names of the requested contrasts
    .contrast_names = NULL,
    #' @field facade_name registry key stamped into the modelName column
    facade_name = NULL,
    #' @field .drop_na_diff drop rows with NA diff (limma family); default FALSE
    .drop_na_diff = FALSE,
    #' @description get contrast results, stamped with the facade key
    #' @param ... passed to the inner contrast object's get_contrasts
    get_contrasts = function(...) {
      res <- .stamp_facade_identity(self$contrast$get_contrasts(...), self$facade_name)
      if (isTRUE(self$.drop_na_diff) && "diff" %in% colnames(res)) {
        res <- res[!is.na(res$diff), , drop = FALSE]
      }
      res
    },
    #' @description get subject x contrast pairs that could not be estimated
    get_missing = function() {
      .compute_missing(self$.lfqdata, self$.contrast_names, self$get_contrasts())
    },
    #' @description get \code{\link{ContrastsPlotter}}
    #' @param ... passed to the inner contrast object's get_Plotter
    get_Plotter = function(...) self$contrast$get_Plotter(...),
    #' @description convert results to wide format
    #' @param ... passed to the inner contrast object's to_wide
    to_wide = function(...) self$contrast$to_wide(...)
  )
)
