# LFQDataAggregator ----
#'
#' Decorates LFQData with methods to aggregate protein intensities
#' aggregates intensities
#'
#' @export
#' @family LFQData
#'
#' @examples
#' istar <-  prolfqua::sim_lfq_data_peptide_config()
#' istar$config <- istar$config
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(data, istar$config)
#'
#' lfqTrans <- lfqdata$clone()$get_Transformer()
#' lfqTrans$log2()
#' lfqTrans <- lfqTrans$robscale()$lfq
#' lfqAggregator <- LFQDataAggregator$new(lfqTrans, "protein")
#'
#' lfqAggregator$medpolish()
#' pmed <- lfqAggregator$plot()
#' pmed$plots[[1]]
#' lfqAggregator$lmrob()
#' prob <- lfqAggregator$plot()
#' prob$plots[[1]]
#'
#' lfqCopy <- lfqdata$clone()
#' lfqCopy$is_transformed()
#' lfqAggregator <- LFQDataAggregator$new(lfqCopy, "protein")
#' lfqAggregator$sum_topN()
#' pSum <- lfqAggregator$plot()
#' pSum$plots[[1]]
#'
#' lfqAggregator$mean_topN()
#' pMean <- lfqAggregator$plot()
#' pMean$plots[[1]]
#' protPlotter <- lfqAggregator$lfq_agg$get_Plotter()
#' protPlotter$heatmap()
#' \dontrun{
#' lfqAggregator$write_plots(tempdir())
#' }
#'
LFQDataAggregator <- R6::R6Class(
  "LFQDataAggregator",
  public = list(
    #' @field lfq LFQData
    lfq = NULL,
    #' @field lfq_agg aggregation result
    lfq_agg = NULL,
    #' @field prefix to use for aggregation results e.g. protein
    prefix = character(),
    ## @field filepath
    ## filepath = character(),
    #' @description
    #' initialize
    #' @param lfq LFQData
    #' @param prefix default protein
    initialize = function(lfq, prefix = "protein"){
      if ( length(lfq$config$table$hierarchy_keys()) == 1 ) {
        stop("no hierarchies to aggregate from: ",  lfq$config$table$hierarchy_keys())
      }
      if (length(lfq$config$table$hierarchy_keys()) == lfq$config$table$hierarchyDepth) {
        stop("no hierarchies to aggregate from: ",
             lfq$config$table$hierarchy_keys(),
             ", hierarchyDepth :",
             lfq$config$table$hierarchyDepth)
      }
      self$lfq = lfq$clone(deep = TRUE)
      self$prefix = prefix
    },
    #' @description
    #' aggregate using median polish
    #' @param N top N by intensity
    #' @return LFQData
    medpolish = function(){
      if (!self$lfq$is_transformed()) {
        warning("You did not transform the intensities.",
                "medpolish works best with already variance stabilized intensities.",
                "Use LFQData$get_Transformer to transform the data :",
                self$lfq$config$table$workIntensity)
      }
      res <- estimate_intensity(self$lfq$data, self$lfq$config, .func = medpolish_estimate_dfconfig)
      self$lfq_agg <- LFQData$new(res$data, res$config, prefix = self$prefix)
      invisible(self$lfq_agg)
    },
    #' @description
    #' aggregate using robust regression
    #' @param N top N by intensity
    #' @return LFQData
    lmrob = function(){
      if (!self$lfq$is_transformed()) {
        warning("You did not transform the intensities.",
                "Robust regression works best with already variance stabilized intensities.",
                "Use LFQData$get_Transformer to transform the data.",
                self$lfq$config$table$workIntensity,)
      }

      res <- estimate_intensity(self$lfq$data, self$lfq$config, .func = rlm_estimate_dfconfig)
      res <-
      self$lfq_agg <- LFQData$new(res$data, res$config, prefix = self$prefix)
      invisible(self$lfq_agg)
    },
    #' @description
    #' aggregate topN using mean
    #' @param N top N by intensity
    #' @return LFQData
    mean_topN = function(N = 3){
      mean_f <- function(x, name = FALSE){
        if (name) {return("mean")}
        return(mean(x, na.rm = TRUE))
      }
      private$.topN(N = N, .func = mean_f)
    },
    #' @description
    #' aggregate topN using sum
    #' @param N top N by intensity
    #' @return LFQData
    sum_topN = function(N = 3){
      sum_f <- function(x, name = FALSE){
        if (name) { return("sum") }
        sum(x, na.rm = TRUE)}

      private$.topN(N = N, .func = sum_f)
    },
    #' @description
    #' creates aggregation plots
    #' @param subset create plots for a subset of the data only, e.g. proteins with more then 2 peptides.
    #' @param show.legend default FALSE
    #' @return data.frame
    #'
    plot = function(subset = NULL, show.legend = FALSE){
      if (is.null(self$lfq_agg)) {
        stop("please aggregate the data first")
      }

      if (!is.null(subset)) {
        lfqagg <- self$lfq_agg$get_subset(subset)
      }else {
        lfqagg <- self$lfq_agg
      }
      df <- prolfqua::plot_estimate(
        self$lfq$data,
        self$lfq$config,
        lfqagg$data,
        lfqagg$config,
        show.legend = show.legend)
      invisible(df)
    },
    #' @description
    #' writes plots to folder
    #'
    #' @param qcpath qcpath
    #' @param subset write plots only for some
    #' @param show.legend legend
    #' @param width figure width
    #' @param height figure height
    #' @return file path
    write_plots = function(qcpath, subset = NULL, show.legend = FALSE, width = 6, height = 6){
      pl <- self$plot(subset)
      pb <- progress::progress_bar$new(total = nrow(pl))
      filepath <- file.path(qcpath, paste0(self$prefix, "_aggregation_plot.pdf"))
      pdf(filepath , width = width, height = height)
      for (i in seq_len(nrow(pl))) {
        print(pl$plots[[i]])
        pb$tick()
      }
      dev.off()
      invisible(filepath)
    }
  ),
  private = list(
    .topN = function(.func, N = 3){
      if (self$lfq$is_transformed()) {
        warning("You did transform the intensities.",
                "top N works with raw data.",
                self$lfq$config$table$workIntensity)
      }
      ranked <- rank_peptide_by_intensity(self$lfq$data , self$lfq$config)
      resTOPN <- aggregate_intensity_topN(ranked, self$lfq$config, .func = .func, N = N)
      self$lfq_agg <- LFQData$new(resTOPN$data, resTOPN$config, prefix = self$prefix)
      invisible(self$lfq_agg)
    }
  )
)
