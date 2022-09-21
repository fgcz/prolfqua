# ContrastsModerated -----


#' Limma moderated contrasts
#' @export
#' @family modelling
#' @examples
#'
#' istar <- old2new(prolfqua_data('data_ionstar')$normalized())
#' istar_data <- dplyr::filter(istar$data ,protein_Id %in% sample(protein_Id, 100))
#' pepIntensity <- istar_data
#' config <- istar$config$clone(deep = TRUE)
#'
#'
#' ld <- LFQData$new(pepIntensity, config)
#' lProt <- ld$get_Aggregator()$medpolish()
#' lProt$rename_response("transformedIntensity")
#' modelFunction <-
#'   strategy_lm("transformedIntensity  ~ dilution.")
#' mod <- build_model(
#'  lProt,
#'  modelFunction)
#'
#' Contr <- c("dil.b_vs_a" = "dilution.a - dilution.b")
#' contrast <- prolfqua::Contrasts$new(mod,
#'  Contr)
#' contrast <- ContrastsModerated$new(contrast)
#' bb <- contrast$get_contrasts()
#'
#' csi <- ContrastsSimpleImpute$new(lProt, contrasts = Contr)
#'
#' csi$get_contrasts()
#' contrast$get_contrasts()
#'
#' merged <- addContrastResults(contrast, csi)
#' merged$more$get_contrasts() |> dim()
#' merged$merged$get_contrasts() |> dim()
#' merged$same$get_contrasts() |> dim()
#'
#' cs <- contrast$get_contrast_sides()
#' cslf <- contrast$get_linfct()
#' ctr <- contrast$get_contrasts()
#' ctrwide <- contrast$to_wide()
#' cp <- contrast$get_Plotter()
#' cp$histogram()
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
      contrast_result <- mutate(contrast_result,modelName = self$modelName, .before  = 1)
      contrast_result <- dplyr::ungroup(contrast_result)
      stopifnot(all(.requiredContrastColumns %in% colnames(contrast_result)))

      return(contrast_result)
    },
    #' @description
    #' get \code{\link{Contrasts_Plotter}}
    #' @param FCthreshold fold change threshold to show in plots
    #' @param FDRthreshold FDR threshold to show in plots
    #'
    get_Plotter = function(
    FCthreshold = 1,
    FDRthreshold = 0.1
    ){
      contrast_result <- self$get_contrasts()
      res <- Contrasts_Plotter$new(
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


