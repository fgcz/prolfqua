# ModelInterface -----
#' R6 interface class representing modelling result
#'
#' @return An R6 class generator.
#' @export
#' @examples
#'
#' mi <- ModelInterface$new()
#'
#' testthat::expect_error(mi$get_coefficients())
#' testthat::expect_error(mi$get_anova())
#' testthat::expect_error(mi$coef_histogram())
#' testthat::expect_error(mi$coef_volcano())
#' testthat::expect_error(mi$coef_pairs())
#' testthat::expect_error(mi$anova_histogram())
#'
#'
ModelInterface <- R6::R6Class(
  "ModelInterface",
  public = list(
    #' @description
    #' return model coefficients
    #' @return data.frame
    get_coefficients = function() {
      stop("get_coefficients not implmeneted")
    },
    #' @description
    #' perform ANOVA analysis
    #' @return data.frame
    get_anova = function() {
      stop("get_anova not implmeneted")
    },
    #' @description
    #' histogram of coefficient p-values
    #' @return list with \code{plot} (ggplot) and file \code{name}
    coef_histogram = function() {
      histogram_coeff <- ggplot(private$coef_table(), aes(x = !!sym(private$coef_pvalue), group = factor)) +
        geom_histogram(breaks = seq(0, 1, by = 0.05)) +
        facet_wrap(~factor)
      list(plot = histogram_coeff, name = paste0("Coef_Histogram_", self$model_name, ".pdf"))
    },
    #' @description
    #' volcano plot of non-intercept coefficients
    #' @return list with \code{plot} (ggplot) and file \code{name}
    coef_volcano = function() {
      volcano_plot <- private$coef_table() |>
        dplyr::filter(factor != "(Intercept)") |>
        prolfqua::multigroup_volcano(
          effect = "Estimate",
          significance = private$coef_pvalue,
          contrast = "factor",
          label = "subject_id",
          xintercept = c(-1, 1),
          colour = private$volcano_colour
        )
      list(plot = volcano_plot, name = paste0(private$volcano_prefix, self$model_name, ".pdf"))
    },
    #' @description
    #' coefficient estimates in wide format, one column per coefficient, for a pairs plot
    #' @return list with \code{plot} (data.frame) and file \code{name}
    coef_pairs = function() {
      for_pairs <- private$coef_table() |>
        dplyr::select(all_of(c("subject_id", "factor", "Estimate"))) |>
        tidyr::pivot_wider(names_from = "factor", values_from = "Estimate")
      list(plot = for_pairs, name = paste0("Coef_Pairsplot_", self$model_name, ".pdf"))
    },
    #' @description
    #' histogram of ANOVA p-values or FDR
    #' @param what show either "p.value" or "FDR"
    #' @return list with \code{plot} (ggplot) and file \code{name}
    anova_histogram = function(what = c("p.value", "FDR")) {
      what <- match.arg(what)
      histogram_anova <- ggplot(self$get_anova(), aes(x = !!sym(what), group = factor)) +
        geom_histogram(breaks = seq(0, 1, by = 0.05)) +
        facet_wrap(~factor)
      list(plot = histogram_anova, name = paste0("Anova_p.values_", self$model_name, ".pdf"))
    }
  ),
  private = list(
    # coefficient p-value column of get_coefficients()
    coef_pvalue = "Pr...t..",
    # colour column of the coefficient volcano plot
    volcano_colour = "isSingular",
    volcano_prefix = "Coef_volcano_plot_",
    coef_table = function() tidyr::unite(self$get_coefficients(), "subject_id", self$subject_id)
  )
)
