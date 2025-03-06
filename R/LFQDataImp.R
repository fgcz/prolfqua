# LFQDataImp ----
#'
#' Decorates LFQData with methods to impute missing values
#'
#' @export
#' @family LFQData
#'
#' @examples
LFQDataImp <- R6::R6Class(
  "LFQDataImp",
  public = list(
    #' @field lfq LFQData
    lfq = NULL,
    initialize = function(lfq, prefix = "protein"){
      self$lfq = lfq$clone(deep = TRUE)
      self$prefix = prefix
    })
)
