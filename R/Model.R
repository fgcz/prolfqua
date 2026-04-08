# Model -----
#' R6 class representing modelling result
#'
#' @export
#' @family modelling
#' @examples
#'
#'
#'
#' istar <- prolfqua_data('data_ionstar')$normalized()
#' istar_data <- dplyr::filter(istar$data ,protein_Id %in% sample(protein_Id, 100))
#' model_name <- "f_condtion_r_peptide"
#' formula_randomPeptide <-
#'   strategy_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id)",
#'    model_name = model_name)
#' pepIntensity <- istar_data
#' config <- istar$config
#' config$hierarchy_keys_depth()
#' mod <- prolfqua::build_model(
#'  pepIntensity,
#'  formula_randomPeptide,
#'  model_name = model_name,
#'  subject_id = config$hierarchy_keys_depth())
#'
#' mod$model_df
#' aovtable  <- mod$get_anova()
#' mod$get_coefficients()
#' mod$coef_histogram()
#' mod$coef_volcano()
#' mod$coef_pairs()
#' mod$anova_histogram()
#'
Model <- R6::R6Class(
  "Model",
  inherit = ModelInterface,
  public = list(
    #' @field model_df data.frame with modelling data and model.
    model_df = NULL,
    #' @field model_name name of model
    model_name = character(),
    #' @field subject_id e.g. protein_Id
    subject_id = character(),
    #' @field model_strategy function to create the models
    model_strategy = NULL,
    #' @field anova_df function to compute anova
    anova_df = NULL,
    #' @field p.adjust function to adjust p-values
    p.adjust = NULL,
    #' @description
    #' initialize
    #' @param model_df dataframe with modelling results
    #' @param model_strategy model_strategy see \code{\link{strategy_lmer}}
    #' @param model_name name of model
    #' @param subject_id subject column name
    #' @param p.adjust method to adjust p-values
    #'
    initialize = function(
      model_df,
      model_strategy,
      model_name,
      subject_id = "protein_Id",
      p.adjust = prolfqua::adjust_p_values
    ) {
      self$model_df = model_df
      self$model_strategy = model_strategy
      self$model_name = model_name
      self$subject_id = subject_id
      self$p.adjust = p.adjust
    },
    #' @description
    #' return model coefficient table
    get_coefficients = function() {
      lmermodel <- "linear_model"
      complete_models <- get_complete_model_fit(self$model_df)
      # Extract coefficients
      .coef_df <- function(x) {
        x <- coef(summary(x))
        x <- data.frame(factor = row.names(x), x)
        return(x)
      }
      model_coeff <- complete_models |>
        dplyr::mutate(!!"Coeffs_model" := purrr::map(!!sym(lmermodel), .coef_df))
      model_coeff <- model_coeff |>
        dplyr::select(!!!syms(self$subject_id), !!sym("Coeffs_model"), isSingular, nr_coef)
      model_coeff <- tidyr::unnest(model_coeff, cols = c(Coeffs_model))
      return(model_coeff)
    },
    #' @description
    #' return anova table
    get_anova = function() {
      lmermodel <- "linear_model"
      complete_models <- get_complete_model_fit(self$model_df)

      model_anova <- complete_models |>
        dplyr::mutate(!!"Anova_model" := purrr::map(!!sym(lmermodel), self$model_strategy$anova_df$fun))

      model_anova <- model_anova |>
        dplyr::select(!!!syms(self$subject_id), !!sym("Anova_model"), isSingular, nr_coef)
      model_anova <- tidyr::unnest(model_anova, cols = c(Anova_model))

      model_anova <- model_anova |> dplyr::filter(factor != "Residuals")
      model_anova <- model_anova |> dplyr::filter(factor != "NULL")

      model_anova <- self$p.adjust(model_anova, column = self$model_strategy$anova_df$col_pval, group_by_col = "factor")

      model_anova <- model_anova |> dplyr::rename(p.value = !!sym(self$model_strategy$anova_df$col_pval))
      return(dplyr::ungroup(model_anova))
    },

    #' @description
    #' histogram of model coefficient
    coef_histogram = function() {
      model_coeff <- self$get_coefficients()
      model_coeff <- tidyr::unite(model_coeff, "subject_id", self$subject_id)
      ## Coef_Histogram
      fname_histogram_coeff <- paste0("Coef_Histogram_", self$model_name, ".pdf")
      histogram_coeff <- ggplot(data = model_coeff, aes(x = Pr...t.., group = factor)) +
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
          significance = "Pr...t..",
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
      ## Anova_p.values
      what <- match.arg(what)
      model_anova <- self$get_anova()
      fname_histogram_anova <- paste0("Anova_p.values_", self$model_name, ".pdf")
      histogram_anova <- model_anova |>
        ggplot(aes(x = !!sym(what), group = factor)) +
        geom_histogram(breaks = seq(0, 1, by = 0.05)) +
        facet_wrap(~factor)
      return(list(plot = histogram_anova, name = fname_histogram_anova))
    },
    #' @description
    #' write figures related to ANOVA into pdf file
    #' @param path folder name
    #' @param width figure width
    #' @param height figure height
    #'
    write_anova_figures = function(path, width = 10, height = 10) {
      private$write_fig(self$anova_histogram(), path, width, height)
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
