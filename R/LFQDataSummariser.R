
# LFQDataSummariser ----
#' Summarize LFQData
#'
#' @export
#' @family LFQData
#' @examples
#' library(prolfqua)
#' istar <- prolfqua::sim_lfq_data_peptide_config()
#'
#' data <- istar$data
#' lfqdata <- LFQData$new(data, istar$config)
#' sum <- lfqdata$get_Summariser()
#' sum
#' sum$hierarchy_counts()
#' sum$hierarchy_counts_sample("wide")
#' sum$hierarchy_counts_sample("long")
#' sum$plot_hierarchy_counts_sample()
#' sum$plot_hierarchy_counts_sample()
#' tmp <- sum$interaction_missing_stats()
#'
#' sum$missingness_per_group()
#' sum$missingness_per_group_cumsum()
#' sum$plot_missingness_per_group()
#' sum$plot_missingness_per_group_cumsum()
#' sum$upset_interaction_missing_stats()
#' sum$percentage_abundance()
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
    #' @param nr_children get summary for 1,2 or more number of children
    hierarchy_counts_sample = function(value=c("wide","long"), nr_children = 1){
      value <- match.arg(value)
      fun <- prolfqua::hierarchy_counts_sample(self$lfq$data, self$lfq$config,
                                               nr_children = nr_children)
      return(fun(value))
    },
    #' @description
    #' barplot showing number of elements at each level in every sample
    #' @param value wide - wide format, long - long format, plot - ggplot
    #' @param nr_children get summary for 1,2 or more number of children
    plot_hierarchy_counts_sample = function(nr_children = 1){
      fun <- prolfqua::hierarchy_counts_sample(self$lfq$data, self$lfq$config,
                                               nr_children = nr_children)
      return(fun("plot"))
    },
    #' @description
    #' missing per condition and protein
    interaction_missing_stats = function(){
      prolfqua::interaction_missing_stats(self$lfq$data, self$lfq$config)
    },
    #' @description
    #' upset plot with missing information per protein and condition
    #' @param tr if less than tr observations in condition then missing
    upset_interaction_missing_stats = function(tr = 2){
      pups <- UpSet_interaction_missing_stats(self$lfq$data,  self$lfq$config, tr = tr)
      res <- UpSetR::upset(pups$data, order.by = "freq", nsets = pups$nsets)
      return(res)
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
    },
    #' @description
    #' Does roll up to highest hierarchy and
    #' Computes the percent abundance of proteins overall and within each group
    #' @param N default 1000
    #' @return data frame
    percentage_abundance = function(N = 1000){
      # roll up to protein intensities

      if (self$lfq$is_transformed()) {
        warning("The abundances are transformed.\n Since this function sums up
                protein abundances,\n it is best to use untransformed data.")
      }

      ag <- try(self$lfq$get_Aggregator())
      if (class(ag)[1] == "try-error") {
        bb <- self$lfq$get_copy()
        bb$rename_response("totalIntensity")
      } else{
        bb <- ag$sum_topN(N = N)
        bb$rename_response("totalIntensity")
      }

      # compute protein level summaries
      dall <- interaction_missing_stats(bb$data, bb$config, factors = NULL)
      dfac <- interaction_missing_stats(bb$data, bb$config)
      xd <- setdiff(colnames(dfac$data), colnames(dall$data))
      for (i in xd) {
        dall$data[[i]] <- "ALL"
      }
      all <- dplyr::bind_rows(dfac$data, dall$data)
      nested <- all |> dplyr::group_by(!!!rlang::syms(self$lfq$config$table$factor_keys_depth())) |> tidyr::nest()
      for (i in seq_len(nrow(nested))) {
        nested$data[[i]] <- nested$data[[i]] |>
          dplyr::arrange(.data$meanAbundance) |>
          dplyr::mutate(id = dplyr::row_number()) |>
          dplyr::mutate(abundance_percent = meanAbundance/sum(meanAbundance, na.rm = TRUE)*100 ) |>
          dplyr::mutate(abundance_percent_cumulative = cumsum(ifelse(is.na(abundance_percent), 0, abundance_percent)) + abundance_percent*0) |>
          dplyr::mutate(percent_prot = id / max(id) * 100)
      }
      res <- tidyr::unnest(nested, cols = "data")
      return(res)
    }
  )
)
