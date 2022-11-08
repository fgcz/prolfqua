# ContrastsTable -----

#'
#' holds results when contrasts are added.
#'
#' @export
#' @family modelling
#' @examples
#'
#' bb <- prolfqua_data("data_ionstar")$normalized()
#' bb$config <- old2new(bb$config)
#' configur <- bb$config$clone(deep = TRUE)
#' configur$table$hierarchyDepth <- 2
#' data <- bb$data
#' lfqdata <- LFQData$new(data, configur)
#' Contrasts <- c(
#'   "dilution.b-a" = "dilution.b - dilution.a",
#'   "dilution.c-e" = "dilution.c - dilution.b"
#' )
#' csi <- ContrastsSimpleImpute$new(lfqdata, contrasts = Contrasts)
#' ctr <- csi$get_contrasts()
#' csi$subject_Id
#' xcx <- ContrastsTable$new(ctr, subject_Id = csi$subject_Id, modelName = "TableTest")
#' xcx$get_contrasts()
#' xcx$get_Plotter()$volcano()
#' stopifnot(is.null(xcx$get_contrast_sides()))
#' stopifnot(is.null(xcx$get_linfct()))
#' stopifnot(ncol(xcx$to_wide()) == 8)
#'
ContrastsTable <- R6::R6Class(
  "ContrastsTable",
  inherit = ContrastsInterface,
  public = list(
    #' @field contrast_result contrast results
    contrast_result = NULL,
    #' @field modelName model name
    modelName = character(),
    #' @field subject_Id default protein_Id
    subject_Id = character(),
    #' @description
    #' intitialize
    #' @param contrastsdf data.frame
    #' @param subject_Id default protein_Id
    #' @param modelName default ContrastTable
    initialize = function(contrastsdf,
                          subject_Id = "protein_Id",
                          modelName = "ContrastTable") {
      self$contrast_result <- contrastsdf
      self$subject_Id <- subject_Id
      self$modelName <- modelName
    },
    #' @description
    #' return sides of contrast
    #' @return data.frame
    get_contrast_sides = function() {
      NULL
    },
    #' @description not implemented
    get_linfct = function() {
      NULL
    },

    #' @description
    #' get contrasts
    #' @seealso \code{\link{summary_ROPECA_median_p.scaled}}
    #' @param all should all columns be returned (default FALSE)
    #' @param global use a global linear function (determined by get_linfct)
    get_contrasts = function(all = FALSE) {
      self$contrast_result
    },
    #' @description
    #' get \code{\link{ContrastsPlotter}}
    #' @param FCthreshold fold change threshold
    #' @param FDRthreshold fdr threshold
    #' @return \code{\link{ContrastsPlotter}}
    #'
    get_Plotter = function(FCthreshold = 1, FDRthreshold = 0.1) {
      res <- ContrastsPlotter$new(
        self$contrast_result,
        subject_Id = self$subject_Id,
        fcthresh = FCthreshold,
        volcano = list(
          list(score = "p.value", xlim = c(0, 1, 0.05)),
          list(score = "FDR", thresh = FDRthreshold)
        ),
        histogram = list(
          list(score = "p.value", xlim = c(0, 1, 0.05)),
          list(score = "FDR", xlim = c(0, 1, 0.05))
        ),
        modelName = "modelName",
        diff = "diff",
        contrast = "contrast"
      )
      return(res)
    },
    #' @description convert to wide format
    #' @param columns value column default beta.based.significance
    #' @return data.frame
    to_wide = function(columns = c("p.value", "FDR")) {
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(contrast_minimal,
        subject_Id = self$subject_Id,
        columns = c("diff", columns),
        contrast = "contrast"
      )
      return(contrasts_wide)
    }
  )
)
