
# ContrastsMissing----
#' Compute contrasts with group mean imputation
#'
#' If there are no observations in one of the groups for some of the proteins,
#' the group mean cannot be estimated. Therefore, assuming that the observation
#' is missing because the protein abundance is below the detection limit,
#' we substitute the unobserved group with the median of protein abundances
#'  observed only in one sample of the group.
#' The variance of a protein is estimated using the pooled variance
#' of all observations of all groups.
#'
#' @family modelling
#' @export
#' @examples
#' Nprot <- 120
#' istar <- prolfqua::sim_lfq_data_protein_config(Nprot = Nprot,weight_missing = .4)
#' istar$data$abundance |> is.na() |> sum()
#' protIntensity <- istar$data
#' config <- istar$config
#'
#'
#' lProt <- LFQData$new(protIntensity, config)
#' lProt$rename_response("transformedIntensity")
#'
#' Contr <- c("dil.b_vs_a" = "group_A - group_Ctrl")
#' csi <- ContrastsMissing$new(lProt, contrasts = Contr)
#' csi$get_contrast_sides()
#'
#' res <- csi$get_contrasts()
#'
#' stopifnot(nrow(res) ==  (protIntensity$protein_Id |> unique() |> length()))
#' res$contrast |> table()
#' stopifnot((res$p.value |> is.na() |> sum()) == 0)
#' plot(res$diff, -log10(res$p.value), pch = ".")
#' csi$column_description()
#' x<- csi$get_Plotter()
#' p <- x$volcano()
#' pdf(file = NULL)
#' print(p)
#' dev.off()
#'
#' dd <- prolfqua::sim_lfq_data_2Factor_config(Nprot = 100,weight_missing = 0.1)
#'
#' Contrasts <- c("c1" = "TreatmentA - TreatmentB",
#'                "C2" = "BackgroundX- BackgroundZ",
#'                "c3" = "`TreatmentA:BackgroundX` - `TreatmentA:BackgroundZ`",
#'                "c4" = "`TreatmentB:BackgroundX` - `TreatmentB:BackgroundZ`"
#'                )
#' lProt <- LFQData$new(dd$data, dd$config)
#' lProt$rename_response("transformedIntensity")
#'
#' csi <- ContrastsMissing$new(lProt, contrasts = Contrasts)
#' res <- csi$get_contrasts()
#' pl <- csi$get_Plotter()
#' pdf(file = NULL)
#' pl$volcano()
#' dev.off()
ContrastsMissing <- R6::R6Class(
  "ContrastsMissing",
  inherit = ContrastsInterface,
  private = list(
    method = "V1"
  ),
  public = list(
    #' @field subject_Id subject_id e.g. protein_ID column
    subject_Id = character(),
    #' @field contrasts array with contrasts (see example)
    contrasts = character(),
    #' @field modelName model name
    modelName = character(),
    #' @field contrast_result data frame with results of contrast computation
    contrast_result = NULL,
    #' @field lfqdata data frame
    lfqdata = NULL,
    #' @field confint confidence interval
    confint = 0.95,
    #' @field p.adjust function to adjust p-values
    p.adjust = NULL,
    #' @field global Take global or local values for imputation
    global = logical(),
    #' @field present default 1, presence in interaction to infer limit of detection.
    present = 1,
    #' @field minsd default 1, if standard deviation can not be estimated, what is the prior minimum sd, default = 1s
    minsd = 1,
    #' @description
    #' initialize
    #' @param lfqdata LFQData
    #' @param contrasts array of contrasts (see example)
    #' @param confint confidence interval
    #' @param p.adjust method for p-value adjustment - default Benjamini Hochberg
    #' @param modelName default "groupAverage"
    initialize = function(lfqdata,
                          contrasts,
                          confint = 0.95,
                          p.adjust = prolfqua::adjust_p_values,
                          modelName = "groupAverage"
                          ){
      self$subject_Id = lfqdata$config$table$hierarchy_keys_depth()
      self$contrasts = contrasts
      self$modelName = modelName
      self$lfqdata = lfqdata
      self$confint = confint
      self$p.adjust = p.adjust
    },
    #' @description
    #' get contrasts sides
    #'
    get_contrast_sides = function(){
      # extract contrast sides
      tt <- self$contrasts[grep("-",self$contrasts)]
      tt <- tibble(contrast = names(tt) , rhs = tt)
      tt <- tt |> dplyr::mutate(rhs = gsub("[` ]","",rhs)) |>
        tidyr::separate(rhs, c("group_1", "group_2"), sep = "-")
      return(tt)
    },
    #' @description
    #' table with results of contrast computation
    #' @param all FALSE, do not show all columns (default)
    get_contrasts = function(all = FALSE){
      if (is.null(self$contrast_result)) {
        if (self$lfqdata$config$table$hierarchyDepth < length(self$lfqdata$config$table$hierarchy_keys())) {
          stop("hierarchy depth < hierarchy_keys(). Please aggregate first.")
        } else {
          mh1 <- prolfqua::MissingHelpers$new(self$lfqdata$data, self$lfqdata$config, prob = 0.5, weighted = TRUE)
          result <- mh1$get_contrasts(Contrasts = self$contrasts, confint = self$confint, all = all)
          result <- self$p.adjust(result, column = "p.value", group_by_col = "contrast", newname = "FDR")
        }
        result <- result |> rename(diff = estimate, sigma = sd, std.error = sdT )
        result <- mutate(result, modelName = self$modelName, .before = 1)
        self$contrast_result <- ungroup(result)
      }
      res <- self$contrast_result
      stopifnot(all(names(super$column_description()$column_name) %in% colnames(res)))
      invisible(res)
    },
    #' @description
    #' get ContrastsPlotter
    #' @return Contrast_Plotter
    get_Plotter = function(){
      res <- ContrastsPlotter$new(
        self$get_contrasts(),
        subject_Id = self$subject_Id,
        volcano = list(list(score = "p.value", thresh = 0.1),list(score = "FDR", thresh = 0.1)),
        histogram = list(list(score = "p.value", xlim = c(0,1,0.05)),
                         list(score = "FDR", xlim = c(0,1,0.05))),
        modelName = "modelName",
        diff = "diff",
        contrast = "contrast")
      return(res)
    },
    #' @description
    #' convert contrast results to wide format
    #' @param columns value column default p.value
    #' @return data.frame
    to_wide = function(columns = c("p.value", "FDR","statistic")){
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(contrast_minimal,
                                                     subject_Id = self$subject_Id,
                                                     columns = c("diff", columns),
                                                     contrast = 'contrast')
      return(contrasts_wide)
    }

  )
)
