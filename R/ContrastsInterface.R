# ContrastsInterface ----
#' Base class for all Contrasts classes
#' @export
ContrastsInterface <- R6::R6Class(
  "ContrastsInterface",
  public = list(
    #' @description
    #' get table with sides of the contrast
    get_contrast_sides = function(){stop("get_contrast_sides not implemented.")},
    #' @description
    #' get table with contrast results (similar to limma topTable function)
    get_contrasts = function(){stop("get_contrasts not implemented.")},
    #' @description
    #' initialize plotter
    get_Plotter = function(){stop("get_Plotter not implemented.")},
    #' @description
    #' create wide representation of data.
    to_wide = function(){stop("to_wide not implemented.")},
    #' @description
    #' column description
    column_description = function() {
      description <- c(
        "modelName" = "type of model",
        "contrast" = "name of difference e.g. group1_vs_group2",
        "avgAbd" = "mean abundance value of protein in all samples",
        "diff" = "difference among conditions",
        "FDR" = "false discovery rate",
        "statistic" = "t-statistics",
        "std.error" = "standard error",
        "df" = "degrees of freedom",
        "p.value" = "p-value",
        "conf.low" = "lower value of 95 confidence interval",
        "conf.high" = "high value of 95 confidence interval",
        "sigma" = "residual standard deviation of linear model (needed for empirical Bayes variance shrinkage)."
      )
      description <- data.frame(column_name = names(description), description = description)
      return(description)
    }
  )
)

.requiredContrastColumns <- c("contrast",
                              "diff",
                              "FDR",
                              "df",
                              "sigma",
                              "statistic",
                              "p.value",
                              "conf.low",
                              "conf.high")







# Merge contrasts ----
#' add contrast results from two different functions. Typically used with Contrast and Cotnrast simple imputed.
#'
#' @param prefer contrasts to use preferentially
#' @param add contrasts to add from if missing in prefer
#' @param modelName name of the merged model default "mergedModel"
#' @export
#' @family modelling
#'
addContrastResults <- function(prefer, add, modelName = "mergedModel"){
  cA <- prefer$get_contrasts()
  cB <- add$get_contrasts()

  stopifnot(length(setdiff(colnames(cA), colnames(cB))) == 0)

  cA <- dplyr::filter(cA,!is.na(.data$statistic))
  moreID <- setdiff(distinct(select(cB, c(prefer$subject_Id, "contrast"))),
                    distinct(select(cA, c(add$subject_Id, "contrast"))))
  more <- inner_join(moreID, cB )

  sameID <- select(cA, c(add$subject_Id, "contrast"))
  same <- inner_join(sameID, cB)

  merged <- bind_rows(cA, more)

  if (prefer$modelName == add$modelName) {
    prefermodelName <- paste0(prefer$modelName, "_prefer")
    addmodelName <- paste0(add$modelName, "_add")
    cA$modelName <- prefermodelName
    more$modelName <- addmodelName
  } else {
    prefermodelName <- prefer$modelName
    addmodelName <- add$modelName
  }

  merged$modelName <- factor(merged$modelName,
                             levels = c(levels(factor(cA$modelName)), addmodelName))

  merged <- ContrastsTable$new(merged,
                               subject_Id = prefer$subject_Id,
                               modelName = paste0(prefermodelName, "_", addmodelName))
  more <- ContrastsTable$new(more, subject_Id = prefer$subject_Id, modelName = addmodelName)
  same <-  ContrastsTable$new(same, subject_Id = prefer$subject_Id, modelName = addmodelName)
  return(list(merged = merged, more = more, same = same))
}


