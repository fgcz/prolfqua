# ContrastsROPECA -----

#'
#' ROPECA reproducibility-optimization method
#'
#' ROPECA optimizes the reproducibility of statistical testing
#'  on peptide-level and aggregates the peptide-level changes
#'   to determine differential protein-level expression.
#'
#' @export
#' @family modelling
#' @examples
#'
#' istar <- old2new(prolfqua_data('data_ionstar')$normalized())
#' istar_data <- dplyr::filter(istar$data ,protein_Id %in% sample(protein_Id, 100))
#' modelFunction <-
#'   strategy_lm("transformedIntensity  ~ dilution.")
#' pepIntensity <- istar_data
#' config <- istar$config$clone(deep = TRUE)
#' config$table$hierarchyDepth <- 2
#' config$table$hkeysDepth()
#'
#' mod <- build_model(
#'  pepIntensity,
#'  modelFunction,
#'  subject_Id = config$table$hkeysDepth())
#'
#'  Contr <- c("dil.b_vs_a" = "dilution.b - dilution.a")
#'
#'
#'  contr <- prolfqua::Contrasts$new(mod, Contr)
#'  dim(contr$get_contrasts())
#'  contrM <- prolfqua::ContrastsModerated$new(contr)
#'  dim(contrM$get_contrasts())
#'  contrast <- prolfqua::ContrastsROPECA$new(contrM)
#'  contrast$get_contrasts()
#'  #ContrastsROPECA$debug("to_wide")
#'  contrast <- prolfqua::ContrastsROPECA$new(contr)
#'  tmp <- contrast$get_contrasts()
#'  dim(tmp)
#'  pl <- contrast$get_Plotter()
#'  contrast$to_wide()
#'  pl$histogram()
#'  pl$ma_plot()
#'
ContrastsROPECA <- R6::R6Class(
  "ContrastsROPECA",
  inherit = ContrastsInterface,
  public = list(
    #' @field Contrast Contrast
    Contrast = NULL,
    #' @field contrast_result contrast result
    contrast_result = NULL,
    #' @field modelName model name
    modelName = character(),
    #' @field subject_Id columns with protein ID's
    subject_Id = character(),
    #' @field p.adjust method to use for p.value adjustment
    p.adjust = NULL,
    #' @description
    #' initialize
    #' @param Contrast e.g. instance of Contrasts class, or ContrastsModerated
    #' @param modelName default ROPECA
    #' @param p.adjust function to use for p.value adjustement
    initialize = function(Contrast,
                          modelName = "ROPECA",
                          p.adjust = prolfqua::adjust_p_values
    ){
      self$Contrast = Contrast
      stopifnot(length(Contrast$subject_Id) > 1)
      self$modelName = modelName
      self$subject_Id = Contrast$subject_Id
      self$p.adjust = p.adjust
    },
    #' @description
    #' show names of contrasts
    #' @return data.frame
    get_contrast_sides = function(){
      self$contrast$get_contrast_sides()
    },
    #' @description
    #' get linear function used to determine contrasts
    #' @return data.frame
    get_linfct = function(){
      self$contrast$get_linfct()
    },
    #' @description
    #' get contrasts
    #' @seealso \code{\link{summary_ROPECA_median_p.scaled}}
    #' @param all should all columns be returned (default FALSE)
    #' @param global use a global linear function (determined by get_linfct)
    #' @return data.frame
    get_contrasts = function(all = FALSE){
      if (is.null(self$contrast_result)) {
        contrast_result <- self$Contrast$get_contrasts(all = FALSE)
        contrast_result <- summary_ROPECA_median_p.scaled(
          contrast_result,
          contrast = "contrast",
          subject_Id = self$subject_Id[length(self$subject_Id) - 1],
          estimate = "diff",
          statistic = "statistic",
          p.value = "p.value",
          max.n = 10)
        contrast_result <- dplyr::rename(contrast_result, diff = "estimate")
        contrast_result <- self$p.adjust(contrast_result,
                                         column = "beta.based.significance",
                                         group_by_col = "contrast",
                                         newname = "FDR.beta.based.significance")
        contrast_result <- self$p.adjust(contrast_result,
                                         column = "median.p.value",
                                         group_by_col = "contrast",
                                         newname = "FDR.median.p.value")
        contrast_result <- mutate(contrast_result,modelName = self$modelName, .before  = 1)
        self$contrast_result <- contrast_result

      }

      if (!all) {
        res <- select(self$contrast_result ,
                      -all_of(c( "n_not_na", "mad.estimate",
                                 "n.beta", "isSingular",
                                 "median.p.scaled","median.p.value",
                                 "FDR.median.p.value")) )
      }else{
        res <- self$contrast_result
      }

      return(res)
    },
    #' @description
    #' get \code{\link{Contrasts_Plotter}}
    #' @return \code{\link{Contrasts_Plotter}}
    #' @param FDRthreshold FDR threshold
    #' @param FCthreshold FC threshold
    get_Plotter = function(FDRthreshold = 0.1,
                           FCthreshold = 2){
      contrast_result <- self$get_contrasts()
      res <- Contrasts_Plotter$new(
        contrast_result,
        subject_Id = self$subject_Id[1],
        fcthresh = FCthreshold,
        volcano = list(list(score = "FDR.beta.based.significance", thresh = FDRthreshold)),
        histogram = list(list(score = "beta.based.significance", xlim = c(0,1,0.05)),
                         list(score = "FDR.beta.based.significance", xlim = c(0,1,0.05))),
        modelName = "modelName",
        diff = "diff",
        contrast = "contrast")
      return(res)
    },
    #' @description convert to wide format
    #' @param columns value column default beta.based.significance
    #' @return data.frame
    to_wide = function(columns = c("beta.based.significance", "FDR.beta.based.significance")){
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(
        contrast_minimal,
        subject_Id = self$subject_Id[length(self$subject_Id) - 1],
        columns = c("diff", columns),
        contrast = 'contrast')
      return(contrasts_wide)
    }

  ))
