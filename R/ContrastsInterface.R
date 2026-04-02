# ContrastsInterface ----
#' Base class for all Contrasts classes
#' @export
#' @examples
#' int <- ContrastsInterface$new()
#' testthat::expect_error(int$get_contrast_sides())
#' testthat::expect_error(int$get_contrasts())
#' testthat::expect_error(int$get_missing())
#' testthat::expect_error(int$get_Plotter())
#' testthat::expect_error(int$to_wide())
#' int$column_description()
ContrastsInterface <- R6::R6Class(
  "ContrastsInterface",
  public = list(
    #' @description
    #' get table with sides of the contrast
    get_contrast_sides = function() {
      stop("get_contrast_sides not implemented.")
    },
    #' @description
    #' get table with contrast results (similar to limma topTable function)
    get_contrasts = function() {
      stop("get_contrasts not implemented.")
    },
    #' @description
    #' initialize plotter
    get_Plotter = function() {
      stop("get_Plotter not implemented.")
    },
    #' @description
    #' create wide representation of data.
    to_wide = function() {
      stop("to_wide not implemented.")
    },
    #' @description
    #' get protein × contrast pairs that could not be estimated.
    #' Returns a data.frame with hierarchy columns and contrast, or 0 rows.
    get_missing = function() {
      stop("get_missing not implemented.")
    },
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


# Merge contrasts ----
#' Merge contrast results coming from two different model.
#'
#' Typically used with results of \code{\link{Contrasts}} and \code{\link{ContrastsMissing}}
#'
#' @param prefer contrasts to use preferentially
#' @param add contrasts to add from if missing in prefer
#' @param modelName name of the merged model default "mergedModel"
#' @export
#' @family modelling
#'
merge_contrasts_results <- function(prefer, add, modelName = "mergedModel") {
  c_a <- prefer$get_contrasts()
  c_b <- add$get_contrasts()

  if (length(colnames(c_a)) < length(colnames(c_b))) {
    c_b <- dplyr::select(c_b, colnames(c_a))
  }

  c_a <- dplyr::filter(c_a, !is.na(.data$statistic))
  more_id <- setdiff(
    distinct(select(c_b, c(prefer$subject_Id, "contrast"))),
    distinct(select(c_a, c(add$subject_Id, "contrast")))
  )
  more <- inner_join(more_id, c_b)

  same_id <- select(c_a, c(add$subject_Id, "contrast"))
  same <- inner_join(same_id, c_b)

  merged <- bind_rows(c_a, more)

  if (prefer$modelName == add$modelName) {
    prefer_model_name <- paste0(prefer$modelName, "_prefer")
    add_model_name <- paste0(add$modelName, "_add")
    c_a$modelName <- prefer_model_name
    more$modelName <- add_model_name
  } else {
    prefer_model_name <- prefer$modelName
    add_model_name <- add$modelName
  }

  merged$modelName <- factor(merged$modelName, levels = c(levels(factor(c_a$modelName)), add_model_name))

  merged <- ContrastsTable$new(
    merged,
    subject_Id = prefer$subject_Id,
    modelName = paste0(prefer_model_name, "_", add_model_name)
  )
  more <- ContrastsTable$new(more, subject_Id = prefer$subject_Id, modelName = add_model_name)
  same <- ContrastsTable$new(same, subject_Id = prefer$subject_Id, modelName = add_model_name)
  return(list(merged = merged, more = more, same = same))
}
