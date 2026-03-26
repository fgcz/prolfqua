# LFQDataAggregator ----

# Shared helper: validate hierarchy is aggregatable
.check_aggregatable <- function(lfq) {
  if (length(lfq$config$hierarchy_keys()) == 1) {
    stop("no hierarchies to aggregate from: ", lfq$config$hierarchy_keys())
  }
  if (length(lfq$config$hierarchy_keys()) == lfq$config$hierarchyDepth) {
    stop(
      "no hierarchies to aggregate from: ",
      lfq$config$hierarchy_keys(),
      ", hierarchyDepth :",
      lfq$config$hierarchyDepth
    )
  }
}

# Shared helper: plot aggregation result
.aggregator_plot <- function(lfq, lfq_agg, subset = NULL, show.legend = FALSE) {
  if (is.null(lfq_agg)) {
    stop("please aggregate the data first")
  }
  if (!is.null(subset)) {
    lfqagg <- lfq_agg$get_subset(subset)
  } else {
    lfqagg <- lfq_agg
  }
  df <- prolfqua::plot_estimate(
    lfq$data,
    lfq$config,
    lfqagg$data,
    lfqagg$config,
    show.legend = show.legend
  )
  invisible(df)
}

# Shared helper: write aggregation plots to PDF
.aggregator_write_plots <- function(self, qcpath, subset = NULL, show.legend = FALSE, width = 6, height = 6) {
  pl <- self$plot(subset)
  pb <- progress::progress_bar$new(total = nrow(pl))
  filepath <- file.path(qcpath, paste0(self$prefix, "_aggregation_plot.pdf"))
  pdf(filepath, width = width, height = height)
  for (i in seq_len(nrow(pl))) {
    print(pl$plots[[i]])
    pb$tick()
  }
  dev.off()
  invisible(filepath)
}


#' AggregateMedpolish
#'
#' Aggregates peptide intensities to protein level using median polish.
#' Works best with variance-stabilized (log-transformed) intensities.
#'
#' @export
#' @family LFQData
#'
#' @examples
#' istar <- prolfqua::sim_lfq_data_peptide_config()
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(data, istar$config)
#' lfqTrans <- lfqdata$clone()$get_Transformer()$log2()$robscale()$lfq
#'
#' agg <- AggregateMedpolish$new(lfqTrans, "protein")
#' agg$aggregate()
#' p <- agg$plot()
#' p$plots[[1]]
#'
AggregateMedpolish <- R6::R6Class(
  "AggregateMedpolish",
  public = list(
    #' @field lfq LFQData
    lfq = NULL,
    #' @field lfq_agg aggregation result
    lfq_agg = NULL,
    #' @field prefix to use for aggregation results e.g. protein
    prefix = character(),
    #' @description
    #' initialize
    #' @param lfq LFQData
    #' @param prefix default protein
    initialize = function(lfq, prefix = "protein") {
      .check_aggregatable(lfq)
      if (!lfq$is_transformed()) {
        warning(
          "You did not transform the intensities. ",
          "medpolish works best with already variance stabilized intensities. ",
          "Use LFQData$get_Transformer to transform the data: ",
          lfq$config$workIntensity
        )
      }
      self$lfq <- lfq$clone(deep = TRUE)
      self$prefix <- prefix
    },
    #' @description
    #' run median polish aggregation
    #' @return LFQData
    aggregate = function() {
      res <- estimate_intensity(self$lfq$data, self$lfq$config, .func = medpolish_estimate_dfconfig)
      self$lfq_agg <- LFQData$new(res$data, res$config, prefix = self$prefix)
      invisible(self$lfq_agg)
    },
    #' @description
    #' creates aggregation plots
    #' @param subset create plots for a subset of the data only
    #' @param show.legend default FALSE
    #' @return data.frame
    plot = function(subset = NULL, show.legend = FALSE) {
      .aggregator_plot(self$lfq, self$lfq_agg, subset = subset, show.legend = show.legend)
    },
    #' @description
    #' writes plots to folder
    #' @param qcpath qcpath
    #' @param subset write plots only for some
    #' @param show.legend legend
    #' @param width figure width
    #' @param height figure height
    #' @return file path
    write_plots = function(qcpath, subset = NULL, show.legend = FALSE, width = 6, height = 6) {
      .aggregator_write_plots(self, qcpath, subset, show.legend, width, height)
    }
  )
)


#' AggregateRlm
#'
#' Aggregates peptide intensities to protein level using robust regression (rlm).
#' Works best with variance-stabilized (log-transformed) intensities.
#'
#' @export
#' @family LFQData
#'
#' @examples
#' istar <- prolfqua::sim_lfq_data_peptide_config()
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(data, istar$config)
#' lfqTrans <- lfqdata$clone()$get_Transformer()$log2()$robscale()$lfq
#'
#' agg <- AggregateRlm$new(lfqTrans, "protein")
#' agg$aggregate()
#' p <- agg$plot()
#' p$plots[[1]]
#'
AggregateRlm <- R6::R6Class(
  "AggregateRlm",
  public = list(
    #' @field lfq LFQData
    lfq = NULL,
    #' @field lfq_agg aggregation result
    lfq_agg = NULL,
    #' @field prefix to use for aggregation results e.g. protein
    prefix = character(),
    #' @description
    #' initialize
    #' @param lfq LFQData
    #' @param prefix default protein
    initialize = function(lfq, prefix = "protein") {
      .check_aggregatable(lfq)
      if (!lfq$is_transformed()) {
        warning(
          "You did not transform the intensities. ",
          "Robust regression works best with already variance stabilized intensities. ",
          "Use LFQData$get_Transformer to transform the data. ",
          lfq$config$workIntensity
        )
      }
      self$lfq <- lfq$clone(deep = TRUE)
      self$prefix <- prefix
    },
    #' @description
    #' run robust regression aggregation
    #' @return LFQData
    aggregate = function() {
      res <- estimate_intensity(self$lfq$data, self$lfq$config, .func = rlm_estimate_dfconfig)
      self$lfq_agg <- LFQData$new(res$data, res$config, prefix = self$prefix)
      invisible(self$lfq_agg)
    },
    #' @description
    #' creates aggregation plots
    #' @param subset create plots for a subset of the data only
    #' @param show.legend default FALSE
    #' @return data.frame
    plot = function(subset = NULL, show.legend = FALSE) {
      .aggregator_plot(self$lfq, self$lfq_agg, subset = subset, show.legend = show.legend)
    },
    #' @description
    #' writes plots to folder
    #' @param qcpath qcpath
    #' @param subset write plots only for some
    #' @param show.legend legend
    #' @param width figure width
    #' @param height figure height
    #' @return file path
    write_plots = function(qcpath, subset = NULL, show.legend = FALSE, width = 6, height = 6) {
      .aggregator_write_plots(self, qcpath, subset, show.legend, width, height)
    }
  )
)


#' AggregateTopN
#'
#' Aggregates peptide intensities to protein level using top N peptides.
#' Works with raw (untransformed) intensities.
#'
#' @export
#' @family LFQData
#'
#' @examples
#' istar <- prolfqua::sim_lfq_data_peptide_config()
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(data, istar$config)
#'
#' agg <- AggregateTopN$new(lfqdata, "protein", N = 3, func = "sum")
#' agg$aggregate()
#' p <- agg$plot()
#' p$plots[[1]]
#'
#' agg_mean <- AggregateTopN$new(lfqdata, "protein", N = 3, func = "mean")
#' agg_mean$aggregate()
#' protPlotter <- agg_mean$lfq_agg$get_Plotter()
#' protPlotter$heatmap()
#' agg_mean$write_plots(tempdir())
#'
AggregateTopN <- R6::R6Class(
  "AggregateTopN",
  public = list(
    #' @field lfq LFQData
    lfq = NULL,
    #' @field lfq_agg aggregation result
    lfq_agg = NULL,
    #' @field prefix to use for aggregation results e.g. protein
    prefix = character(),
    #' @field N top N peptides by intensity
    N = 3L,
    #' @field func aggregation function name: "sum" or "mean"
    func = "sum",
    #' @description
    #' initialize
    #' @param lfq LFQData
    #' @param prefix default protein
    #' @param N top N peptides (default 3)
    #' @param func "sum" or "mean" (default "sum")
    initialize = function(lfq, prefix = "protein", N = 3, func = "sum") {
      .check_aggregatable(lfq)
      if (lfq$is_transformed()) {
        warning("You did transform the intensities. top N works with raw data. ", lfq$config$workIntensity)
      }
      match.arg(func, c("sum", "mean"))
      self$lfq <- lfq$clone(deep = TRUE)
      self$prefix <- prefix
      self$N <- N
      self$func <- func
    },
    #' @description
    #' run top N aggregation
    #' @return LFQData
    aggregate = function() {
      .func <- if (self$func == "sum") {
        function(x, name = FALSE) {
          if (name) {
            return("sum")
          }
          sum(x, na.rm = TRUE)
        }
      } else {
        function(x, name = FALSE) {
          if (name) {
            return("mean")
          }
          mean(x, na.rm = TRUE)
        }
      }
      ranked <- rank_peptide_by_intensity(self$lfq$data, self$lfq$config)
      resTOPN <- aggregate_intensity_topN(ranked, self$lfq$config, .func = .func, N = self$N)
      self$lfq_agg <- LFQData$new(resTOPN$data, resTOPN$config, prefix = self$prefix)
      invisible(self$lfq_agg)
    },
    #' @description
    #' creates aggregation plots
    #' @param subset create plots for a subset of the data only
    #' @param show.legend default FALSE
    #' @return data.frame
    plot = function(subset = NULL, show.legend = FALSE) {
      .aggregator_plot(self$lfq, self$lfq_agg, subset = subset, show.legend = show.legend)
    },
    #' @description
    #' writes plots to folder
    #' @param qcpath qcpath
    #' @param subset write plots only for some
    #' @param show.legend legend
    #' @param width figure width
    #' @param height figure height
    #' @return file path
    write_plots = function(qcpath, subset = NULL, show.legend = FALSE, width = 6, height = 6) {
      .aggregator_write_plots(self, qcpath, subset, show.legend, width, height)
    }
  )
)
