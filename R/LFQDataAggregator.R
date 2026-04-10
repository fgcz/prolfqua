# LFQDataAggregator ----

# Shared helper: validate hierarchy is aggregatable
.check_aggregatable <- function(lfq) {
  if (length(lfq$hierarchy_keys()) == 1) {
    stop("no hierarchies to aggregate from: ", lfq$hierarchy_keys())
  }
  if (length(lfq$hierarchy_keys()) == length(lfq$relevant_hierarchy_keys())) {
    stop(
      "no hierarchies to aggregate from: ",
      lfq$hierarchy_keys(),
      ", hierarchyDepth :",
      length(lfq$relevant_hierarchy_keys())
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
    lfq,
    lfqagg,
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
          lfq$config$work_intensity
        )
      }
      self$lfq <- lfq$clone(deep = TRUE)
      self$prefix <- prefix
    },
    #' @description
    #' run median polish aggregation
    #' @return LFQData
    aggregate = function() {
      res <- estimate_intensity(self$lfq, .func = medpolish_estimate_dfconfig)
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
          lfq$config$work_intensity
        )
      }
      self$lfq <- lfq$clone(deep = TRUE)
      self$prefix <- prefix
    },
    #' @description
    #' run robust regression aggregation
    #' @return LFQData
    aggregate = function() {
      res <- estimate_intensity(self$lfq, .func = rlm_estimate_dfconfig)
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
        warning("You did transform the intensities. top N works with raw data. ", lfq$config$work_intensity)
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
      ranked <- rank_peptide_by_intensity(
        self$lfq$get_data(),
        self$lfq$response(),
        self$lfq$hierarchy_keys()
      )
      res_topn <- aggregate_intensity_top_n(ranked, self$lfq, .func = .func, N = self$N)
      self$lfq_agg <- LFQData$new(res_topn$data, res_topn$config, prefix = self$prefix)
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


# Helper: convert limpa EList back to long-format LFQData
.elist_to_lfqdata <- function(elist, wide, config, prefix, impute_only) {
  annotation <- wide$annotation
  sample_col <- config$sample_name

  intensity_name <- "limpa"
  se_name <- "limpa_se"

  # Row IDs from the EList (united hierarchy keys with ~lfq~ separator from to_wide)
  row_ids_raw <- rownames(elist$E)
  sample_names <- colnames(elist$E)

  if (impute_only) {
    hierarchy_keys <- config$hierarchy_keys()
  } else {
    hierarchy_keys <- config$hierarchy_keys_depth()
  }

  # Parse row IDs back to hierarchy columns
  sep <- "~lfq~"
  all_keys <- c(hierarchy_keys, config$isotope_label)
  if (length(all_keys) == 1) {
    row_ids <- data.frame(V1 = row_ids_raw, stringsAsFactors = FALSE)
    colnames(row_ids) <- all_keys
  } else {
    row_ids <- as.data.frame(
      do.call(rbind, strsplit(row_ids_raw, sep, fixed = TRUE)),
      stringsAsFactors = FALSE
    )
    if (ncol(row_ids) == length(all_keys)) {
      colnames(row_ids) <- all_keys
    } else {
      colnames(row_ids) <- hierarchy_keys
    }
  }

  # Pivot three matrices to long format
  nr <- nrow(elist$E)
  nc <- ncol(elist$E)

  row_idx <- rep(seq_len(nr), each = nc)
  col_idx <- rep(seq_len(nc), times = nr)

  long_data <- row_ids[row_idx, , drop = FALSE]
  rownames(long_data) <- NULL
  long_data[[sample_col]] <- sample_names[col_idx]
  long_data[[intensity_name]] <- as.vector(t(elist$E))
  long_data[[se_name]] <- as.vector(t(elist$other$standard.error))

  n_obs_vec <- as.vector(t(elist$other$n.observations))

  # Join annotation (factors, fileName, etc.)
  # Exclude columns already parsed from row IDs to avoid .x/.y suffixes
  anno_cols <- intersect(
    colnames(annotation),
    c(config$file_name, sample_col, config$factor_keys(), config$isotope_label, config$norm_value)
  )
  anno_cols <- setdiff(anno_cols, c(colnames(long_data), sample_col))
  anno_cols <- c(sample_col, anno_cols)
  long_data <- dplyr::left_join(long_data, annotation[, anno_cols, drop = FALSE], by = sample_col)

  # Build new config
  if (impute_only) {
    newconfig <- config$clone(deep = TRUE)
    newconfig$work_intensity <- intensity_name
  } else {
    newconfig <- make_reduced_hierarchy_config(
      config,
      work_intensity = intensity_name,
      hierarchy = config$hierarchy_keys_depth(names = FALSE)
    )
  }

  # nr_children from n_observations
  nr_children_name <- paste0("nr_children_", paste(hierarchy_keys, collapse = "_"))
  long_data[[nr_children_name]] <- n_obs_vec
  newconfig$nr_children <- nr_children_name

  # opt_se for standard errors
  newconfig$opt_se <- se_name

  LFQData$new(long_data, newconfig, prefix = prefix)
}


#' AggregateLimpa
#'
#' Aggregates peptide/precursor intensities to protein level using limpa's
#' probabilistic DPC-based quantification (\code{limpa::dpcQuant}).
#' With \code{impute_only = TRUE}, performs same-level imputation without
#' aggregation using \code{limpa::dpcQuantByRow}.
#'
#' The output LFQData contains three value columns:
#' \itemize{
#'   \item intensity (\code{"limpa"}) — protein-level log-expression (complete, no NAs)
#'   \item standard error (in \code{config$opt_se}) — per-protein, per-sample posterior SEs
#'   \item observation count (in \code{config$nr_children}) — number of observed precursors
#' }
#'
#' @export
#' @family LFQData
#'
#' @examples
#' \dontrun{
#' istar <- prolfqua::sim_lfq_data_peptide_config()
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#'
#' agg <- AggregateLimpa$new(lfqdata, "protein")
#' agg$aggregate()
#' agg$lfq_agg$to_wide()
#' }
AggregateLimpa <- R6::R6Class(
  "AggregateLimpa",
  public = list(
    #' @field lfq LFQData (deep cloned input)
    lfq = NULL,
    #' @field lfq_agg aggregation result
    lfq_agg = NULL,
    #' @field prefix to use for aggregation results e.g. protein
    prefix = character(),
    #' @field dpc_result estimated DPC object from limpa::dpc
    dpc_result = NULL,
    #' @field dpc_slope DPC slope parameter (default 0.8)
    dpc_slope = 0.8,
    #' @field impute_only if TRUE use dpcQuantByRow (no aggregation)
    impute_only = FALSE,

    #' @description
    #' initialize
    #' @param lfq LFQData with log2-transformed intensities
    #' @param prefix default "protein"
    #' @param dpc_slope DPC slope parameter passed to limpa (default 0.8)
    #' @param impute_only if TRUE, use dpcQuantByRow instead of dpcQuant
    initialize = function(lfq, prefix = "protein", dpc_slope = 0.8, impute_only = FALSE) {
      if (!requireNamespace("limpa", quietly = TRUE)) {
        stop(
          "Package 'limpa' is required for AggregateLimpa. ",
          "Install it from Bioconductor: BiocManager::install('limpa')"
        )
      }
      if (!impute_only) {
        .check_aggregatable(lfq)
      }
      self$lfq <- lfq$clone(deep = TRUE)
      self$prefix <- prefix
      self$dpc_slope <- dpc_slope
      self$impute_only <- impute_only
    },

    #' @description
    #' run limpa DPC-based aggregation (or imputation)
    #' @return LFQData
    aggregate = function() {
      wide <- self$lfq$to_wide(as.matrix = TRUE)
      expr_matrix <- wide$data

      # Estimate DPC
      self$dpc_result <- limpa::dpc(expr_matrix, b1.upper = 1)

      if (self$impute_only) {
        elist <- limpa::dpcQuantByRow(
          expr_matrix,
          dpc = self$dpc_result,
          dpc.slope = self$dpc_slope,
          verbose = FALSE
        )
        # Restore rownames from input (dpcQuantByRow uses integer IDs internally)
        rownames(elist$E) <- rownames(expr_matrix)
        if (!is.null(elist$other$standard.error)) {
          rownames(elist$other$standard.error) <- rownames(expr_matrix)
        }
        if (!is.null(elist$other$n.observations)) {
          rownames(elist$other$n.observations) <- rownames(expr_matrix)
        }
      } else {
        # Build protein ID vector from hierarchy_keys_depth in rowdata
        rowdata <- wide$rowdata
        hierarchy_keys_depth <- self$lfq$relevant_hierarchy_keys()
        if (length(hierarchy_keys_depth) == 1) {
          protein_ids <- rowdata[[hierarchy_keys_depth]]
        } else {
          protein_ids <- do.call(
            paste,
            c(rowdata[, hierarchy_keys_depth, drop = FALSE], sep = "~lfq~")
          )
        }
        elist <- limpa::dpcQuant(
          expr_matrix,
          protein.id = protein_ids,
          dpc = self$dpc_result,
          dpc.slope = self$dpc_slope,
          verbose = FALSE
        )
      }

      self$lfq_agg <- .elist_to_lfqdata(
        elist,
        wide,
        self$lfq$config,
        self$prefix,
        self$impute_only
      )
      invisible(self$lfq_agg)
    },

    #' @description
    #' creates aggregation plots (only for aggregation mode, not impute_only)
    #' @param subset create plots for a subset of the data only
    #' @param show.legend default FALSE
    #' @return data.frame
    plot = function(subset = NULL, show.legend = FALSE) {
      if (self$impute_only) {
        message("Aggregation plots not available in impute_only mode")
        return(invisible(NULL))
      }
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
      if (self$impute_only) {
        message("Aggregation plots not available in impute_only mode")
        return(invisible(NULL))
      }
      .aggregator_write_plots(self, qcpath, subset, show.legend, width, height)
    }
  )
)
