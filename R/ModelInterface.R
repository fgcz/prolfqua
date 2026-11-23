# ModelInterface -----
#' R6 interface class representing modelling result
#'
#' @return An R6 class generator.
#' @export
#' @examples
#'
#' mi <- ModelInterface$new()
#'
#' abstract <- c("get_coefficients", "get_anova", "coef_histogram", "coef_volcano", "coef_pairs",
#'   "anova_histogram")
#' for (m in abstract) {
#'   stopifnot(inherits(try(mi[[m]](), silent = TRUE), "try-error"))
#' }
#'
#'
ModelInterface <- R6::R6Class(
  "ModelInterface",
  public = list(
    #' @description
    #' problems recorded while fitting: failed fits, fit warnings (e.g. non-convergence), singular fits, fits refit
    #' at the limit of detection, missing responses dropped
    #' @return tibble with columns step, event, n, ids (list), reason, failure, response
    provenance = function() {
      self$get_provenance_log()$entries()
    },
    #' @description
    #' internal: the model's provenance record
    #' @return Provenance object
    get_provenance_log = function() {
      if (is.null(private$.provenance)) {
        private$.provenance <- Provenance$new()
      }
      private$.provenance
    },
    #' @description
    #' return model coefficients
    #' @return data.frame
    get_coefficients = function() {
      abort_prolfqua("`get_coefficients()` is not implemented for this class.")
    },
    #' @description
    #' perform ANOVA analysis
    #' @return data.frame
    get_anova = function() {
      abort_prolfqua("`get_anova()` is not implemented for this class.")
    },
    #' @description
    #' histogram of coefficient p-values
    #' @return list with \code{plot} (ggplot) and file \code{name}
    coef_histogram = function() {
      histogram_coeff <- ggplot(private$coef_table(), aes(x = .data[[private$coef_pvalue]], group = factor)) +
        theme_prolfqua() +
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
        .multigroup_volcano_labelled(
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
      histogram_anova <- ggplot(self$get_anova(), aes(x = .data[[what]], group = factor)) +
        theme_prolfqua() +
        geom_histogram(breaks = seq(0, 1, by = 0.05)) +
        facet_wrap(~factor)
      list(plot = histogram_anova, name = paste0("Anova_p.values_", self$model_name, ".pdf"))
    }
  ),
  private = list(
    .provenance = NULL,
    # coefficient p-value column of get_coefficients()
    coef_pvalue = "Pr...t..",
    # colour column of the coefficient volcano plot
    volcano_colour = "isSingular",
    volcano_prefix = "Coef_volcano_plot_",
    coef_table = function() tidyr::unite(self$get_coefficients(), "subject_id", dplyr::all_of(self$subject_id))
  )
)
