# ContrastsSaintExpress ----
#' ContrastsSaintExpress Wrapper to results produced by SaintExpress (list.txt file)
#'
#' @export
#' @family modelling
#'
ContrastsSaintExpress <- R6::R6Class(
  "ContrastsSaint",
  inherit = ContrastsInterface,
  public = list(
    #' @field contrast_result data.frame with the contrast computation results
    contrast_result = NULL,
    #' @field modelName model name
    modelName = character(),
    #' @field subject_Id subject id defualt 'Prey'
    subject_Id = character(),
    #' @description
    #' initialize
    #' @param contrastsdf return value of \code{\link{runSaint}}
    #' @param subject_Id default "Prey"
    #' @param modelName name of model
    initialize = function(contrastsdf,
                          subject_Id = "Prey",
                          modelName = "ContrastSaint"
    ){


      self$contrast_result = contrastsdf
      self$subject_Id = subject_Id
      self$modelName = modelName

      if ( "AvgIntensity" %in% colnames(contrastsdf)) {
        self$contrast_result <- contrastsdf |> mutate(log2_EFCs = log2(FoldChange),
                                                      avgAbd = log2(AvgIntensity),
                                                      modelName = modelName)

      }else{
        self$contrast_result <- contrastsdf |> mutate(log2_EFCs = log2(FoldChange),
                                                      avgAbd = log2(AvgSpec) ,
                                                      modelName = modelName)

      }

    },
    #' @description
    #' show contrasts
    #' @return data.frame
    get_contrast_sides = function(){
      # extract contrast sides
      tt <- self$contrasts[grep("-",self$contrasts)]
      tt <- tibble(contrast = names(tt) , rhs = tt)
      tt <- tt |> mutate(rhs = gsub("[` ]","",rhs)) |>
        tidyr::separate(rhs, c("group_1", "group_2"), sep = "-")
      return(tt)
    },
    #' @description
    #' no available for SaintExpress
    #'
    get_linfct = function(){
      NULL
    },

    #' @description
    #' get contrasts
    #' @seealso \code{\link{summary_ROPECA_median_p.scaled}}
    #' @param all should all columns be returned (default FALSE)
    #' @param global use a global linear function (determined by get_linfct)
    get_contrasts = function(all = FALSE){
      res <- self$contrast_result |> select(
        all_of(c(self$subject_Id,
                 "modelName",
                 "Bait",
                 "avgAbd",
                 "log2_EFCs",
                 "SaintScore",
                 "BFDR"
        )))
      res
    },
    #' @description get \code{\link{Contrasts_Plotter}}
    #' @param FCthreshold fold change threshold to show
    #' @param SaintScore SaintScore threshold to show in the heatmap.
    #' @param BFDRthreshold BDRF threshold
    #' @return \code{\link{Contrasts_Plotter}}
    get_Plotter = function(FCthreshold = 1, SaintScore = 0.75, BFDRthreshold = 0.1){
      res <- Contrasts_Plotter$new(
        self$contrast_result,
        subject_Id = self$subject_Id,
        fcthresh = FCthreshold,
        volcano = list(list(score = "BFDR", thresh = BFDRthreshold)),
        histogram = list(list(score = "BFDR", xlim = c(0,1,0.05)), list(score = "SaintScore", xlim = c(0,1,0.05))),
        score = list(list(score = "SaintScore", thresh = SaintScore )),
        modelName = "modelName",
        diff = "log2_EFCs",
        contrast = "Bait")
      return(res)
    },
    #' @description convert to wide format
    #' @param columns value column default SaintScore, BFDR
    #' @return data.frame
    to_wide = function(columns = c("SaintScore", "BFDR")){
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(contrast_minimal,
                                                     subject_Id = self$subject_Id,
                                                     columns = c("log2_EFCs", columns),
                                                     contrast = 'Bait')
      return(contrasts_wide)
    }
  ))
