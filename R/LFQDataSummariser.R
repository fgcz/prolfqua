
# LFQDataSummariser ----
#' Summarize LFQData
#'
#' @export
#' @family LFQData
#' @examples
#'
#' istar <- prolfqua_data('data_ionstar')
#' istar$config <- old2new(istar$config)
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(data, istar$config)
#' sum <- lfqdata$get_Summariser()
#' sum
#' sum$hierarchy_counts()
#' sum$hierarchy_counts_sample("wide")
#' sum$hierarchy_counts_sample("long")
#' sum$plot_hierarchy_counts_sample()
#' sum$interaction_missing_stats()
#' sum$missingness_per_group()
#' sum$missingness_per_group_cumsum()
#' sum$plot_missingness_per_group()
#' sum$plot_missingness_per_group_cumsum()
#'
LFQDataSummariser <- R6::R6Class(
  "LFQDataSummariser",
  public = list(
    #' @field lfq LFQData
    lfq = NULL,
    #' @description
    #' initialize
    #' @param lfqdata LFQData
    initialize = function(lfqdata ) {
      self$lfq <- lfqdata
    },
    #' @description
    #' summarize hierarchy
    hierarchy_counts = function(){
      self$lfq$hierarchy_counts()
    },
    #' @description
    #' number of elements at each level in every sample
    #' @param value wide - wide format, long - long format, plot - ggplot
    hierarchy_counts_sample = function(value=c("wide","long")){
      value <- match.arg(value)
      fun <- prolfqua::hierarchy_counts_sample(self$lfq$data, self$lfq$config)
      return(fun(value))
    },
    #' @description
    #' barplot showing number of elements at each level in every sample
    #' @param value wide - wide format, long - long format, plot - ggplot
    plot_hierarchy_counts_sample = function(){
      fun <- prolfqua::hierarchy_counts_sample(self$lfq$data, self$lfq$config)
      return(fun("plot"))
    },
    #' @description
    #' missing per condition and protein
    interaction_missing_stats = function(){
      prolfqua::interaction_missing_stats(self$lfq$data, self$lfq$config)
    },
    #' @description
    #' missing stats per condition
    missingness_per_group = function(){
      prolfqua::missingness_per_condition(self$lfq$data, self$lfq$config)$data
    },
    #' @description
    #' missing stats per condition as cumulative sum
    missingness_per_group_cumsum = function(){
      prolfqua::missingness_per_condition_cumsum(self$lfq$data, self$lfq$config)$data
    },
    #' @description
    #' barplot with number of features with 1,2, etc missing in condition
    #' @return ggplot
    plot_missingness_per_group = function(){
      prolfqua::missingness_per_condition(self$lfq$data, self$lfq$config)$figure
    },
    #' @description
    #' barplot with cumulative sum of features with 1,2, etc missing in condition
    #' @return ggplot
    plot_missingness_per_group_cumsum = function(){
      prolfqua::missingness_per_condition_cumsum(self$lfq$data, self$lfq$config)$figure
    }
  )
)
