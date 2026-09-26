# Model -----
#' R6 class representing modelling result
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#'
#'
#'
#' istar <- sim_lfq_data_peptide_config(Nprot = 20)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#' model_name <- "f_condtion_r_peptide"
#' formula_randomPeptide <-
#'   strategy_lmer(paste0(lfqdata$response(), " ~ group_ + (1 | peptide_Id)"),
#'    model_name = model_name)
#' mod <- prolfqua::build_model(
#'  lfqdata,
#'  formula_randomPeptide,
#'  model_name = model_name)
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
      self$model_df <- model_df
      self$model_strategy <- model_strategy
      self$model_name <- model_name
      self$subject_id <- subject_id
      self$p.adjust <- p.adjust
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
      # Drop degenerate anova rows from failed/empty model fits. The real
      # failure signal is an NA factor; the previous `factor != "NULL"`
      # checked a string sentinel that the anova backends never emit.
      model_anova <- model_anova |> dplyr::filter(!is.na(factor), factor != "NULL")

      model_anova <- self$p.adjust(model_anova, column = self$model_strategy$anova_df$col_pval, group_by_col = "factor")

      model_anova <- model_anova |> dplyr::rename(p.value = !!sym(self$model_strategy$anova_df$col_pval))
      return(dplyr::ungroup(model_anova))
    }
  )
)
