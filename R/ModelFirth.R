# Model -----
#' R6 class representing modelling result
#'
#' @export
#' @family modelling
#' @examplesIf requireNamespace("logistf", quietly = TRUE)
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
#' @return An R6 class generator.

ModelFirth <- R6::R6Class(
  "ModelFirth",
  inherit = ModelInterface,
  public = list(
    #' @field models \code{\link{model_analyse}} result: \code{model_df}, one
    #'   fitted model per subject, and the \code{strategy}
    models = NULL,
    #' @field model_name name of model
    model_name = character(),
    #' @field subject_id e.g. protein_Id
    subject_id = character(),
    #' @field p.adjust function to adjust p-values
    p.adjust = NULL,
    #' @description
    #' initialize
    #' @param models \code{\link{model_analyse}} result with a \code{model_df}
    #'   data.frame and the \code{strategy}
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
      if (!is.list(models) || !is.data.frame(models$model_df)) {
        abort_bad_argument("models", "be a model_analyse() result with a model_df data.frame")
      }
      self$models <- models
      self$model_name <- model_name
      self$subject_id <- subject_id
      self$p.adjust <- p.adjust
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
          dplyr::case_when(object$method.ci == "Wald" ~ 1, object$method.ci == "-" ~ 3, .default = 2)
        )
        dimnames(out) <- list(
          names(object$coefficients),
          c("Estimate", "se(coef)", paste(c("lower", "upper"), 1 - object$alpha), "Chisq", "p", "method")
        )
        out <- data.frame(factor = row.names(out), out)
        return(out)
      }

      model_coeff <- self$models$model_df |>
        dplyr::mutate(!!"Coeffs_model" := purrr::map(!!sym(lmermodel), .coef_df))
      model_coeff <- model_coeff |>
        dplyr::select(!!!syms(self$subject_id), !!sym("Coeffs_model"), isSingular, nr_coef)
      model_coeff <- tidyr::unnest(model_coeff, cols = c(Coeffs_model))
      return(model_coeff)
    },
    #' @description
    #' return anova table
    get_anova = function() {
      warning("method not implemented!")
      return(NULL)
    },
    #' @description
    #' histogram of ANOVA results
    #' @param what show either "Pr..F." or "FDR.Pr..F."
    anova_histogram = function(what = c("p.value", "FDR")) {
      warning("not implemented")
      return(NULL)
    }
  ),
  private = list(
    coef_pvalue = "p"
  )
)
