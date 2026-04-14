# Model -----
#' R6 class representing modelling result
#'
#' @export
#' @family modelling
#' @examples
#'
#'
#'
#'
#' istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 10, with_missing = TRUE,
#'   weight_missing = 0.5, seed = 3)
#' istar$data <- prolfqua::encode_bin_resp(LFQData$new(istar$data, istar$config))
#' istar$config$bin_resp <- "bin_resp"
#' tmp <- LFQData$new(istar$data, istar$config)
#' formula <- paste0(tmp$get_config()$bin_resp , "~ group_")
#' mod <- build_model_logistf(tmp, formula)
#' tmp <- mod$get_coefficients()
#' mod$coef_histogram()
#' mod$coef_pairs()
#' mod$get_anova()
#' mod$coef_volcano()
#' mod$anova_histogram()
#' mod$write_coef_figures(tempdir())
#'
#' istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 10, with_missing = TRUE,
#'   weight_missing = 0.5, seed = 3)
#' istar$data <- prolfqua::encode_bin_resp(LFQData$new(istar$data, istar$config))
#' istar$config$bin_resp <- "bin_resp"
#' tmp <- LFQData$new(istar$data, istar$config)
#' formula <- paste0(tmp$get_config()$bin_resp , "~ group_")
#' mod <- build_model_logistf(tmp, formula)
#' tmp <- mod$get_coefficients()
#' stopifnot(nrow(tmp) == 30)
#' mod$coef_histogram()
#' mod$coef_pairs()
#' mod$get_anova()
#' mod$coef_volcano()
#' mod$anova_histogram()
#' mod$write_coef_figures(tempdir())

ModelFirth <- R6::R6Class(
  "ModelFirth",
  inherit = ModelInterface,
  public = list(
    #' @field models data.frame with modelling data and model.
    models = NULL,
    #' @field model_name name of model
    model_name = character(),
    #' @field subject_id e.g. protein_Id
    subject_id = character(),
    #' @field anova_df function to compute anova
    anova_df = NULL,
    #' @field p.adjust function to adjust p-values
    p.adjust = NULL,
    #' @description
    #' initialize
    #' @param models dataframe with modelling results
    #' @param model_name name of model
    #' @param subject_id subject column name
    #' @param p.adjust method to adjust p-values
    #'
    initialize = function(
      models,
      model_name = "modelFirth",
      subject_id = "protein_Id",
      p.adjust = prolfqua::adjust_p_values
    ) {
      self$models = models
      self$model_name = model_name
      self$subject_id = subject_id
      self$p.adjust = p.adjust
    },
    #' @description
    #' return model coefficient table
    get_coefficients = function() {
      lmermodel <- "linear_model"
      # Extract coefficients
      .coef_df <- function(x) {
        object <- x
        chi2 <- qchisq(1 - object$prob, 1)
        out <- cbind(
          object$coefficients,
          diag(object$var)^0.5,
          object$ci.lower,
          object$ci.upper,
          chi2,
          object$prob,
          ifelse(object$method.ci == "Wald", 1, ifelse(object$method.ci == "-", 3, 2))
        )
        dimnames(out) <- list(
          names(object$coefficients),
          c("Estimate", "se(coef)", paste(c("lower", "upper"), 1 - object$alpha), "Chisq", "p", "method")
        )
        out <- data.frame(factor = row.names(out), out)
        return(out)
      }

      models <- dplyr::bind_rows(self$models[[1]]$model_df, self$models[[2]]$model_df)
      res <- vector(mode = "list", nrow(models))
      model_coeff <- models |>
        dplyr::mutate(!!"Coeffs_model" := purrr::map(!!sym(lmermodel), .coef_df))
      model_coeff <- model_coeff |>
        dplyr::select(!!!syms(self$subject_id), !!sym("Coeffs_model"), isSingular, nr_coef)
      model_coeff <- tidyr::unnest(model_coeff, cols = c(Coeffs_model))
      if (!is.null(self$models$hkey)) {
        model_coeff <- model_coeff |> dplyr::filter(!grepl(self$models$hkey, factor))
      }
      return(model_coeff)
    },
    #' @description
    #' return anova table
    get_anova = function() {
      warning("method not implemented!")
      return(NULL)
    },

    #' @description
    #' histogram of model coefficient
    coef_histogram = function() {
      model_coeff <- self$get_coefficients()
      model_coeff <- tidyr::unite(model_coeff, "subject_id", self$subject_id)
      ## Coef_Histogram
      fname_histogram_coeff <- paste0("Coef_Histogram_", self$model_name, ".pdf")
      histogram_coeff <- ggplot(data = model_coeff, aes(x = p, group = factor)) +
        geom_histogram(breaks = seq(0, 1, by = 0.05)) +
        facet_wrap(~factor)
      return(list(plot = histogram_coeff, name = fname_histogram_coeff))
    },
    #' @description
    #' volcano plot of non intercept coefficients
    coef_volcano = function() {
      model_coeff <- self$get_coefficients()
      model_coeff <- tidyr::unite(model_coeff, "subject_id", self$subject_id)
      fname_volcano_plot <- paste0("Coef_volcano_plot_", self$model_name, ".pdf")
      volcano_plot <- model_coeff |>
        dplyr::filter(factor != "(Intercept)") |>
        prolfqua::multigroup_volcano(
          effect = "Estimate",
          significance = "p",
          contrast = "factor",
          label = "subject_id",
          xintercept = c(-1, 1),
          colour = "isSingular"
        )
      return(list(plot = volcano_plot, name = fname_volcano_plot))
    },
    #' @description
    #' pairs-plot of coefficients
    coef_pairs = function() {
      model_coeff <- self$get_coefficients()
      model_coeff <- tidyr::unite(model_coeff, "subject_id", self$subject_id)
      ## Coef_Pairsplot
      for_pairs <- model_coeff |>
        dplyr::select(all_of(c("subject_id", "factor", "Estimate"))) |>
        tidyr::pivot_wider(names_from = "factor", values_from = "Estimate")
      fname_pairsplot_coef <- paste0("Coef_Pairsplot_", self$model_name, ".pdf")
      return(list(plot = for_pairs, name = fname_pairsplot_coef))
    },
    #' @description
    #' histogram of ANOVA results
    #' @param what show either "Pr..F." or "FDR.Pr..F."
    anova_histogram = function(what = c("p.value", "FDR")) {
      warning("not implemented")
      return(NULL)
    },
    #' @description
    #' write figures related to ANOVA into pdf file
    #' @param path folder name
    #' @param width figure width
    #' @param height figure height
    #'
    write_anova_figures = function(path, width = 10, height = 10) {
      warning("not implemented")
      return(NULL)
    },
    #' @description
    #' write figures related to Coefficients into pdf file
    #' @param path folder name
    #' @param width figure width
    #' @param height figure height
    #'
    write_coef_figures = function(path, width = 10, height = 10) {
      private$write_fig(self$coef_histogram(), path, width, height)
      private$write_fig(self$coef_volcano(), path, width, height)
      private$write_fig(self$coef_pairs(), path, width, height)
    }
  ),
  private = list(
    write_fig = function(res, path, width = 10, height = 10) {
      fpath <- file.path(path, res$name)
      message("Writing figure into : ", fpath, "\n")
      pdf(fpath, width = width, height = height)
      print(res$plot)
      dev.off()
    }
  )
)
