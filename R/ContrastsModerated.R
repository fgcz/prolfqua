# ContrastsModerated -----


#' Limma moderated contrasts
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
#'
#' merged <- merge_contrasts_results(contrast, csi)
#'
#' merged$more$get_contrasts() |> dim()
#' stopifnot(all(dim(merged$merged$get_contrasts() == c(50,13))))
#' stopifnot(all(dim(merged$same$get_contrasts()) == c(50,13)))
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
    #' @field modelName name of model
    modelName = character(),
    #' @field subject_Id columns with subject_Id (proteinID)
    subject_Id = character(),
    #' @field p.adjust function to adjust p-values
    p.adjust = NULL,
    #' @description
    #' initialize
    #' @param Contrast class implementing the ContrastInterface
    #' @param modelName name of the model
    #' @param p.adjust function to adjust p-values - default BH
    initialize = function(Contrast,
                          modelName = paste0(Contrast$modelName, "_moderated"),
                          p.adjust = prolfqua::adjust_p_values
    ){
      self$Contrast = Contrast
      self$subject_Id = Contrast$subject_Id
      self$modelName = modelName
      self$p.adjust = p.adjust
    },
    #' @description
    #' get both sides of contrasts
    get_contrast_sides = function(){
      self$Contrast$get_contrast_sides()
    },
    #' @description
    #' get linear functions from contrasts
    #' @param global logical TRUE - get the a linear functions for all models, FALSE - linear function for each model
    get_linfct = function(global = TRUE){
      self$Contrast$get_linfct()
    },
    #' @description
    #' applies limma moderation
    #' @seealso \code{\link{moderated_p_limma_long}}
    #' @param all should all columns be returned (default FALSE)
    #' @param global use a global linear function (determined by get_linfct)
    get_contrasts = function(all = FALSE){
      contrast_result <- self$Contrast$get_contrasts(all = FALSE)
      contrast_result <- moderated_p_limma_long(
        contrast_result ,
        group_by_col = "contrast",
        estimate = "diff")
      if (!all) {
        contrast_result <- contrast_result |> select(-c( "sigma","df",
                                                         "statistic", "p.value","conf.low","conf.high",
                                                         "FDR",  "moderated.df.prior" ,
                                                         "moderated.var.prior"))
        contrast_result <- contrast_result |> mutate(sigma = sqrt(moderated.var.post),.keep = "unused")
        contrast_result <- contrast_result |> rename(
          conf.low = "moderated.conf.low",
          conf.high = "moderated.conf.high",
          statistic = "moderated.statistic" ,
          df = "moderated.df.total",
          p.value = "moderated.p.value"
        )
        contrast_result <- self$p.adjust(contrast_result, column = "p.value",
                                         group_by_col = "contrast",
                                         newname = "FDR")
      }else{
        contrast_result <- self$p.adjust(contrast_result,
                                         column = "moderated.p.value",
                                         group_by_col = "contrast",
                                         newname = "FDR.moderated")
      }
      contrast_result <- dplyr::ungroup(contrast_result)
      if (class(contrast_result$modelName) == "factor") {
        mname <- factor(paste0(contrast_result$modelName,"_moderated"),
                        levels = paste0(levels(contrast_result$modelName), "_moderated"))
      }else{
        mname <- paste0(contrast_result$modelName,"_moderated")
      }
      contrast_result$modelName <- mname
      stopifnot(all(.requiredContrastColumns %in% colnames(contrast_result)))

      return(contrast_result)
    },
    #' @description
    #' get \code{\link{ContrastsPlotter}}
    #' @param FCthreshold fold change threshold to show in plots
    #' @param FDRthreshold FDR threshold to show in plots
    #'
    get_Plotter = function(
    FCthreshold = 1,
    FDRthreshold = 0.1
    ){
      contrast_result <- self$get_contrasts()
      res <- ContrastsPlotter$new(
        contrast_result,
        subject_Id = self$subject_Id,
        fcthresh = FCthreshold,
        volcano = list(list(score = "FDR", thresh = FDRthreshold)),
        histogram = list(list(score = "p.value", xlim = c(0,1,0.05)),
                         list(score = "FDR", xlim = c(0,1,0.05))),
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
    to_wide = function(columns = c("p.value", "FDR")){
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(contrast_minimal,
                                                     subject_Id = self$subject_Id,
                                                     columns = c("diff", columns),
                                                     contrast = 'contrast')
      return(contrasts_wide)
    }
  )
)


