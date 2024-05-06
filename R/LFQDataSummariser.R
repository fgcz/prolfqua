
# LFQDataSummariser ----
#' Summarize LFQData
#'
#' @export
#' @family LFQData
#' @examples
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
      x1 <- prolfqua::interaction_missing_stats(self$lfq$data, self$lfq$config)
      x2 <- prolfqua::summarize_stats_factors(self$lfq$data, self$lfq$config)
      return(x2)
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
    #' @return data frame
    percentage_abundance = function(){
      # roll up to protein intensities
      # compute protein level summaries

      dall <- prolfqua::summarize_stats_all(self$lfq$data, self$lfq$config)
      dfac <- prolfqua::summarize_stats_factors(self$lfq$data, self$lfq$config)

      all <- dplyr::bind_rows(dfac, dall)
      nested <- all |> dplyr::group_by(!!sym("interaction")) |> tidyr::nest()
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
