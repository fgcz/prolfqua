# ContrastsModerated -----

#' Limma moderated contrasts
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#'
#' istar <- sim_lfq_data_protein_config(Nprot = 50)
#' protIntensity <- istar$data
#' config <- istar$config
#'
#'
#' lProt <- LFQData$new(protIntensity, config)
#' lProt$rename_response("transformedIntensity")
#' modelFunction <-
#'   strategy_lm("transformedIntensity  ~ group_")
#' mod <- build_model(
#'  lProt,
#'  modelFunction)
#'
#' Contr <- c("dil.b_vs_a" = "group_A - group_Ctrl")
#' contrast <- prolfqua::Contrasts$new(mod,
#'  Contr)
#' contrast <- ContrastsModerated$new(contrast)
#' bb <- contrast$get_contrasts()
#' csi <- ContrastsMissing$new(lProt, contrasts = Contr)
#'
#' contrast$get_contrasts() |> dim()
#' (xx <- csi$get_contrasts())   |> dim()
#' merged <- merge_contrasts_results(contrast, csi)
#'
#' merged$more$get_contrasts() |> dim()
#' stopifnot(all(dim(merged$merged$get_contrasts()) == c(50, 14)))
#' stopifnot(all(dim(merged$same$get_contrasts()) == c(49, 14)))
#'
#' cs <- contrast$get_contrast_sides()
#' cslf <- contrast$get_linfct()
#' ctr <- contrast$get_contrasts()
#' ctrwide <- contrast$to_wide()
#' cp <- contrast$get_Plotter()
#'
#' print(cp$histogram()$p.value, vp=NULL)
#' print(cp$histogram()$FDR, vp = NULL)
#'
#' cp$volcano()
#' cp$ma_plot()
#'
#'
#'
ContrastsModerated <- R6::R6Class(
  classname = "ContrastsModerated",
  inherit = ContrastsInterface,
  public = list(
    #' @field Contrast Class implementing the Contrast interface
    Contrast = NULL,
    #' @field model_name name of model
    model_name = character(),
    #' @field subject_id columns with subject_id (proteinID)
    subject_id = character(),
    #' @field p.adjust function to adjust p-values
    p.adjust = NULL,
    #' @field variance_floor optional lower bound for posterior variances
    variance_floor = NULL,
    #' @description
    #' initialize
    #' @param Contrast class implementing the ContrastInterface
    #' @param model_name name of the model
    #' @param p.adjust function to adjust p-values - default BH
    #' @param variance_floor optional positive lower bound for posterior variances
    initialize = function(
      Contrast,
      model_name = paste0(Contrast$model_name, "_moderated"),
      p.adjust = prolfqua::adjust_p_values,
      variance_floor = NULL
    ) {
      self$Contrast <- Contrast
      self$subject_id <- Contrast$subject_id
      self$model_name <- model_name
      self$p.adjust <- p.adjust
      self$variance_floor <- variance_floor
      self$config <- Contrast$get_config()
    },
    #' @description
    #' get both sides of contrasts
    get_contrast_sides = function() {
      self$Contrast$get_contrast_sides()
    },
    #' @description
    #' get linear functions from contrasts
    #' @param global logical TRUE - get the a linear functions for all models, FALSE - linear function for each model
    get_linfct = function(global = TRUE) {
      self$Contrast$get_linfct()
    },
    #' @description
    #' applies limma moderation
    #' @seealso \code{\link{moderated_p_limma_long}}
    #' @param all should all columns be returned (default FALSE)
    #' @param global use a global linear function (determined by get_linfct)
    get_contrasts = function(all = FALSE) {
      contrast_result <- self$Contrast$get_contrasts(all = FALSE)
      contrast_result <- moderated_p_limma_long(
        contrast_result,
        group_by_col = "contrast",
        estimate = "diff",
        variance_floor = self$variance_floor
      )
      contrast_result <- .finalize_moderated_columns(contrast_result, self$p.adjust, all)
      contrast_result <- dplyr::ungroup(contrast_result)
      # modelName and estimate_type pass through from the wrapped contrast
      # unchanged; the moderated-Wald-test wording belongs in the methods text,
      # not in the per-row model identity.
      stopifnot(all(super$column_description()$column_name %in% colnames(contrast_result)))

      return(contrast_result)
    },
    #' @description
    #' get \code{\link{ContrastsPlotter}}
    #' @param fc_threshold fold change threshold to show in plots
    #' @param fdr_threshold FDR threshold to show in plots
    #'
    get_Plotter = function(
      fc_threshold = 1,
      fdr_threshold = 0.1
    ) {
      contrast_result <- self$get_contrasts()
      res <- ContrastsPlotter$new(
        contrast_result,
        subject_id = self$subject_id,
        fcthresh = fc_threshold,
        volcano = list(list(score = "FDR", thresh = fdr_threshold)),
        histogram = list(list(score = "p.value", xlim = c(0, 1, 0.05)), list(score = "FDR", xlim = c(0, 1, 0.05))),
        score = list(list(score = "statistic", thresh = 5)),
        modelName = "modelName",
        diff = "diff",
        contrast = "contrast"
      )
      return(res)
    },
    #' @description
    #' convert to wide format
    #' @param columns value column default moderated.p.value
    #' @return data.frame
    to_wide = function(columns = c("p.value", "FDR", "statistic")) {
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_to_wide(
        contrast_minimal,
        subject_id = self$subject_id,
        columns = c("diff", columns),
        contrast = "contrast"
      )
      return(contrasts_wide)
    }
  )
)
