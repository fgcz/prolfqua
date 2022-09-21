# ContrastsSimpleImpute----
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
#'
#'
#'
ContrastsSimpleImpute <- R6::R6Class(
  "ContrastsSimpleImpute",
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
    #' @description
    #' initialize
    #' @param lfqdata LFQData
    #' @param contrasts array of contrasts (see example)
    #' @param confint confidence interval
    #' @param p.adjust method for p-value adjustment - default Benjamini Hochberg
    #' @param modelName default "groupAverage"
    #' @param method internal default V1
    #' @param global default TRUE use all or per condition data to impute from
    #' @param present in at most how many samples the protein should be observed
    initialize = function(lfqdata,
                          contrasts,
                          confint = 0.95,
                          p.adjust = prolfqua::adjust_p_values,
                          modelName = "groupAverage",
                          method = "V1",
                          global = TRUE,
                          present = 1){
      self$subject_Id = lfqdata$config$table$hkeysDepth()
      self$contrasts = contrasts
      self$modelName = modelName
      self$lfqdata = lfqdata
      self$confint = confint
      self$p.adjust = p.adjust
      private$method = method
      self$global  = global
      self$present = present
    },
    #' @description
    #' get contrasts sides
    #'
    get_contrast_sides = function(){
      # extract contrast sides
      tt <- self$contrasts[grep("-",self$contrasts)]
      tt <- tibble(contrast = names(tt) , rhs = tt)
      tt <- tt |> mutate(rhs = gsub("[` ]","",rhs)) |>
        tidyr::separate(rhs, c("group_1", "group_2"), sep = "-")
      return(tt)
    },
    #' @description
    #' table with results of contrast computation
    #' @param all FALSE, do not show all columns (default)
    get_contrasts = function(all = FALSE){
      if (is.null(self$contrast_result)) {
        if (self$lfqdata$config$table$hierarchyDepth < length(self$lfqdata$config$table$hierarchyKeys())) {
          stop("hierarchy depth < hierarchyKeys(). Please aggregate first.")
        } else {

          result = get_imputed_contrasts(
            self$lfqdata$data,
            self$lfqdata$config,
            self$contrasts,
            present = self$present,
            global = self$global)

          # compute statistics using pooled variance
          result$isSingular <- TRUE
          result <- select(result , -all_of(c("n","estimate_mad")))

          var = summarize_stats(self$lfqdata$data, self$lfqdata$config)

          pooled <- poolvar(var, self$lfqdata$config, method = self$method)
          pooled <- dplyr::select(pooled ,-all_of(c(self$lfqdata$config$table$factor_keys_depth()[1],"var")))

          result <- dplyr::inner_join(result, pooled, by = self$lfqdata$config$table$hkeysDepth())

          resultNA <- result[result$n == 0, ]
          resultnotNa <- result[result$n != 0,]
          meandf <- resultnotNa |> summarize(n = 1, df = 1, sd = mean(sd),sdT = mean(sdT))
          resultNA$n <- 0
          resultNA$df <- 1
          resultNA$sd <- meandf$sd
          resultNA$sdT <- meandf$sdT
          result <- bind_rows(resultNA, resultnotNa)

          result <- dplyr::mutate(result, statistic = .data$estimate_median / .data$sdT,
                                  p.value = 2*pt(abs(.data$statistic), df = .data$df, lower.tail = FALSE))

          prqt <- -qt((1 - self$confint)/2, df = result$df)
          result$conf.low <- result$estimate_median  - prqt * (result$sdT)
          result$conf.high <- result$estimate_median + prqt * (result$sdT)
          result <- self$p.adjust(result, column = "p.value", group_by_col = "contrast", newname = "FDR")

          if (!all) {
            result <- select(result, -all_of( c("isSingular", "not_na" , "mean" ,"n.groups", "n", "meanAll") ) )
          }

        }

        result <- result |> rename(diff = estimate_median, sigma = sd, std.error = sdT )
        result <- mutate(result, modelName = self$modelName, .before = 1)
        self$contrast_result <- ungroup(result)
      }
      res <- self$contrast_result
      stopifnot(all(.requiredContrastColumns %in% colnames(res)))
      invisible(res)
    },
    #' @description
    #' get ContrastsPlotter
    #' @return Contrast_Plotter
    get_Plotter = function(){
      res <- ContrastsPlotter$new(
        self$get_contrasts(),
        subject_Id = self$subject_Id,
        volcano = list(list(score = "FDR", thresh = 0.1)),
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
