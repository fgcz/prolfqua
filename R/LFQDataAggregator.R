# LFQDataAggregator ----

#' Shared helper: validate hierarchy is aggregatable
#' @param lfq LFQData to inspect
#' @noRd
.check_aggregatable <- function(lfq) {
  if (length(lfq$hierarchy_keys()) == 1) {
    abort_bad_argument(
      "lfq",
      "have more than one hierarchy level to aggregate from",
      not = toString(lfq$hierarchy_keys())
    )
  }
  if (length(lfq$hierarchy_keys()) == length(lfq$relevant_hierarchy_keys())) {
    abort_bad_argument(
      "lfq",
      "have a hierarchy_depth below its number of hierarchy levels to aggregate",
      not = paste0(
        "hierarchy_keys = ",
        toString(lfq$hierarchy_keys()),
        "; hierarchy_depth = ",
        length(lfq$relevant_hierarchy_keys())
      )
    )
  }
}

#' Base class of the peptide to protein aggregators
#'
#' Holds the fields, construction, and plotting shared by \code{\link{AggregateMedpolish}},
#' \code{\link{AggregateRlm}}, \code{\link{AggregateTopN}} and \code{\link{AggregateLimpa}}.
#' Subclasses implement \code{aggregate()}.
#'
#' @return An R6 class generator.
#' @family LFQData
#' @keywords internal
AggregatorBase <- R6::R6Class(
  "AggregatorBase",
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
      .assert_model_input(lfq)
      private$.check_transformed(lfq)
      self$lfq <- lfq$clone(deep = TRUE)
      self$prefix <- prefix
    },
    #' @description
    #' creates aggregation plots
    #' @param subset create plots for a subset of the data only
    #' @param show.legend default FALSE
    #' @return data.frame
    plot = function(subset = NULL, show.legend = FALSE) {
      invisible(.plot_estimate(self$lfq, private$.plotted(subset), show.legend = show.legend))
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
      # one plot at a time: memory does not grow with the number of proteins
      lfq_agg <- private$.plotted(subset)
      plot_data <- .estimate_plot_data(self$lfq, lfq_agg)
      estimate_aes <- .estimate_aes(self$lfq$sample_name(), lfq_agg$response())
      pb <- .make_progress(nrow(plot_data))
      filepath <- file.path(qcpath, paste0(self$prefix, "_aggregation_plot.pdf"))
      withr::local_pdf(filepath, width = width, height = height)
      for (i in seq_len(nrow(plot_data))) {
        .render_plot_to_device(.estimate_plot(plot_data, i, self$lfq, estimate_aes, show.legend))
        pb$tick()
      }
      invisible(filepath)
    }
  ),
  private = list(
    .plotted = function(subset) {
      if (is.null(self$lfq_agg)) {
        abort_bad_argument("lfq_agg", "be an aggregated LFQData (call the aggregator first)", not = "NULL")
      }
      if (is.null(subset)) self$lfq_agg else self$lfq_agg$get_subset(subset)
    },
    .check_transformed = function(lfq) NULL,
    .set_result = function(lfq_agg) {
      self$lfq_agg <- lfq_agg
      self$lfq_agg$complete_cases()
      invisible(self$lfq_agg)
    }
  )
)


#' Aggregate peptide intensities by median polish
#'
#' Aggregates peptide intensities to protein level using median polish.
#' Works best with variance-stabilized (log-transformed) intensities.
#'
#' The aggregated data carry the peptide consistency of each protein in each sample, from the median polish
#' residuals of its peptides: \code{consistency_sigma2} (the protein's residual variance, shrunk over proteins with
#' \code{limma::squeezeVar()}), \code{consistency_factor} (the inverse of the sample's residual variance relative to
#' the protein's, shrunk over all protein and sample cells; cells with one peptide get the prior) and
#' \code{consistency_factor_rel} (the factor divided by its mean over the protein's samples). They also carry weight
#' variants built from them (see \code{AnalysisConfiguration$weight_columns}). The facades' default weights multiply
#' in \code{consistency_factor_rel} (see \code{AnalysisConfiguration$use_peptide_consistency}).
#'
#' @return An R6 class generator.
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
  inherit = AggregatorBase,
  public = list(
    #' @field control convergence controls of \code{stats::medpolish}: \code{maxiter} (default 10) and
    #'   \code{eps} (default 0.01; the polish stops when the sum of absolute residuals changes by less than
    #'   \code{eps} times that sum)
    control = NULL,
    #' @field consistency_priors prior df and variance of the peptide consistency, set by \code{aggregate()}: one
    #'   row with \code{protein_df_prior}, \code{protein_var_prior} (protein level) and \code{cell_df_prior},
    #'   \code{cell_var_prior} (protein and sample level); NA when fewer than three proteins or cells have residual
    #'   df
    consistency_priors = NULL,
    #' @description
    #' initialize
    #' @param lfq LFQData
    #' @param prefix default protein
    #' @param control named list merged over the defaults; convergence controls of \code{stats::medpolish}:
    #'   \code{maxiter} (default 10) and
    #'   \code{eps} (default 0.01; the polish stops when the sum of absolute residuals changes by less than
    #'   \code{eps} times that sum)
    initialize = function(lfq, prefix = "protein", control = list()) {
      self$control <- .aggregation_control("medpolish", control)
      super$initialize(lfq, prefix)
    },
    #' @description
    #' run median polish aggregation
    #' @return LFQData
    aggregate = function() {
      res <- .estimate_intensity(self$lfq, method = "medpolish", control = self$control)
      self$consistency_priors <- res$priors
      private$.set_result(.aggregated_with_failures(res, self$prefix, "medpolish"))
    }
  ),
  private = list(
    .check_transformed = function(lfq) {
      if (!lfq$is_transformed()) {
        warning(
          "You did not transform the intensities. ",
          "medpolish works best with already variance stabilized intensities. ",
          "Use LFQData$get_Transformer to transform the data: ",
          lfq$get_config()$work_intensity
        )
      }
    }
  )
)


#' Aggregate peptide intensities by robust linear regression
#'
#' Aggregates peptide intensities to protein level using robust regression (rlm).
#' Works best with variance-stabilized (log-transformed) intensities.
#'
#' The aggregated data carry the peptide consistency columns and weight variants of
#' \code{\link{AggregateMedpolish}}, computed from the rlm residuals.
#'
#' @return An R6 class generator.
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
  inherit = AggregatorBase,
  public = list(
    #' @field control convergence controls of \code{MASS::rlm}: \code{maxit} (default 20) and \code{acc}
    #'   (default 1e-4, relative change of the IRLS residuals)
    control = NULL,
    #' @field consistency_priors prior df and variance of the peptide consistency, set by \code{aggregate()}; see
    #'   \code{\link{AggregateMedpolish}}
    consistency_priors = NULL,
    #' @description
    #' initialize
    #' @param lfq LFQData
    #' @param prefix default protein
    #' @param control named list merged over the defaults; convergence controls of \code{MASS::rlm}: \code{maxit}
    #'   (default 20) and \code{acc}
    #'   (default 1e-4, relative change of the IRLS residuals)
    initialize = function(lfq, prefix = "protein", control = list()) {
      self$control <- .aggregation_control("rlm", control)
      super$initialize(lfq, prefix)
    },
    #' @description
    #' run robust regression aggregation
    #' @return LFQData
    aggregate = function() {
      res <- .estimate_intensity(self$lfq, method = "rlm", control = self$control)
      self$consistency_priors <- res$priors
      private$.set_result(.aggregated_with_failures(res, self$prefix, "rlm"))
    }
  ),
  private = list(
    .check_transformed = function(lfq) {
      if (!lfq$is_transformed()) {
        warning(
          "You did not transform the intensities. ",
          "Robust regression works best with already variance stabilized intensities. ",
          "Use LFQData$get_Transformer to transform the data. ",
          lfq$get_config()$work_intensity
        )
      }
    }
  )
)


#' Aggregate peptide intensities from the top N peptides
#'
#' Aggregates peptide intensities to protein level using top N peptides.
#' Works with raw (untransformed) intensities.
#'
#' @return An R6 class generator.
#' @export
#' @family LFQData
#'
#' @examplesIf requireNamespace("ComplexHeatmap", quietly = TRUE) && requireNamespace("circlize", quietly = TRUE)
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
  inherit = AggregatorBase,
  public = list(
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
      super$initialize(lfq, prefix)
      self$N <- N
      self$func <- match.arg(func, c("sum", "mean"))
    },
    #' @description
    #' run top N aggregation
    #' @return LFQData
    aggregate = function() {
      fun <- match.fun(self$func)
      .func <- function(x, name = FALSE) {
        if (name) {
          return(self$func)
        }
        fun(x, na.rm = TRUE)
      }
      ranked <- .rank_peptide_by_intensity(
        self$lfq$data_long(),
        self$lfq$response(),
        self$lfq$hierarchy_keys()
      )
      res_topn <- .aggregate_intensity_top_n(ranked, self$lfq, .func = .func, N = self$N)
      private$.set_result(LFQData$new(res_topn$data, res_topn$config, prefix = self$prefix))
    }
  ),
  private = list(
    .check_transformed = function(lfq) {
      if (lfq$is_transformed()) {
        warning("You did transform the intensities. top N works with raw data. ", lfq$get_config()$work_intensity)
      }
    }
  )
)


#' Helper: convert limpa EList back to long-format LFQData
#' @param elist limpa EList with expression and uncertainty
#' @param wide wide-format data, row identifiers and sample annotations
#' @param config AnalysisConfiguration defining the column roles
#' @param data long-format input data; the identification confidence of an output feature in a sample is the
#'   maximum over its observed input features
#' @param prefix prefix identifying the output level
#' @param impute_only return imputed input without protein rollup
#' @noRd
.elist_to_lfqdata <- function(elist, wide, config, data, prefix, impute_only) {
  if (impute_only) {
    newconfig <- config$clone(deep = TRUE)
    newconfig$work_intensity <- "limpa"
    newconfig$weight_columns <- character()
    newconfig$consistency_factor <- NULL
  } else {
    newconfig <- .make_reduced_hierarchy_config(config, "limpa", config$hierarchy_keys_depth(names = FALSE))
  }
  # rows of E are the united hierarchy keys; pivot_longer keeps row-major order,
  # which the transposed standard error and observation count matrices match
  long_data <- .response_matrix_as_tibble(elist$E, "limpa", newconfig)
  long_data$limpa_se <- as.vector(t(elist$other$standard.error))
  long_data <- dplyr::left_join(
    long_data,
    wide$annotation,
    by = config$sample_name,
    relationship = "many-to-one"
  )
  newconfig$nr_children <- paste0("nr_children_", paste(newconfig$hierarchy_keys(), collapse = "_"))
  long_data[[newconfig$nr_children]] <- as.vector(t(elist$other$n.observations))
  newconfig$opt_se <- "limpa_se"
  keys <- newconfig$hierarchy_keys()
  confidence <- .aggregate_confidence(
    data,
    config$get_response(),
    c(keys, config$sample_name),
    config$identification_confidence,
    newconfig$identification_confidence
  )
  long_data <- dplyr::left_join(long_data, confidence, by = c(keys, config$sample_name), relationship = "one-to-one")
  long_data <- .add_feature_columns(long_data, keys, newconfig)
  LFQData$new(long_data, newconfig, prefix = prefix)
}

#' Call the limpa detection probability curve with its slope bound
#' @param expr_matrix numeric feature-by-sample expression matrix
#' @param b1.upper upper bound for the detection probability slope
#' @param dpc_fun limpa detection-curve function
#' @noRd
.call_limpa_dpc <- function(expr_matrix, b1.upper = 1, dpc_fun = limpa::dpc) {
  dpc_args <- list(y = expr_matrix)
  if ("b1.upper" %in% names(formals(dpc_fun))) {
    dpc_args$b1.upper <- b1.upper
  }
  do.call(dpc_fun, dpc_args)
}


#' Aggregate or impute intensities with limpa's DPC quantification
#'
#' Aggregates peptide/precursor intensities to protein level using limpa's
#' probabilistic DPC-based quantification (\code{limpa::dpcQuant}).
#' With \code{impute_only = TRUE}, performs same-level imputation without
#' aggregation using \code{limpa::dpcQuantByRow}.
#'
#' The output LFQData contains three value columns:
#' \itemize{
#'   \item intensity (\code{"limpa"}) - protein-level log-expression (complete, no NAs)
#'   \item standard error (in \code{config$opt_se}) - per-protein, per-sample posterior SEs
#'   \item observation count (in \code{config$nr_children}) - number of observed precursors
#' }
#'
#' limpa's iterative steps run with limpa's defaults: \code{limpa::dpc()} estimates the detection probability curve
#' in 2 iterations on a subset of 2000 rows, and \code{limpa::dpcQuant()} optimises each protein with its own
#' settings. Only the slope bound \code{dpc_slope} is exposed. limpa returns no convergence status for these fits, so
#' a changed setting could not be checked; the fit status of limpa analyses reports \code{not_converged} as NA.
#'
#' @return An R6 class generator.
#' @export
#' @family LFQData
#'
#' @examplesIf requireNamespace("limpa", quietly = TRUE)
#' istar <- prolfqua::sim_lfq_data_peptide_config()
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#'
#' agg <- AggregateLimpa$new(lfqdata, "protein")
#' agg$aggregate()
#' agg$lfq_agg$data_wide()
AggregateLimpa <- R6::R6Class(
  "AggregateLimpa",
  inherit = AggregatorBase,
  public = list(
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
      rlang::check_installed("limpa", reason = "for AggregateLimpa.")
      if (!impute_only) {
        .check_aggregatable(lfq)
      }
      .assert_model_input(lfq)
      self$lfq <- lfq$clone(deep = TRUE)
      self$prefix <- prefix
      self$dpc_slope <- dpc_slope
      self$impute_only <- impute_only
    },

    #' @description
    #' run limpa DPC-based aggregation (or imputation)
    #' @return LFQData
    aggregate = function() {
      wide <- self$lfq$data_wide(as.matrix = TRUE)
      expr_matrix <- wide$data

      # Estimate DPC. Older limpa versions do not expose the b1.upper argument.
      self$dpc_result <- .call_limpa_dpc(expr_matrix, b1.upper = 1)

      if (self$impute_only) {
        elist <- limpa::dpcQuantByRow(
          expr_matrix,
          dpc = self$dpc_result,
          dpc.slope = self$dpc_slope,
          verbose = FALSE
        )
        # Restore rownames from input (dpcQuantByRow uses integer IDs internally)
        rownames(elist$E) <- rownames(expr_matrix)
      } else {
        # protein ID: hierarchy keys at depth, united like the rownames of data_wide
        id_cols <- self$lfq$relevant_hierarchy_keys()
        protein_ids <- do.call(paste, c(wide$rowdata[, id_cols, drop = FALSE], sep = "~lfq~"))
        elist <- limpa::dpcQuant(
          expr_matrix,
          protein.id = protein_ids,
          dpc = self$dpc_result,
          dpc.slope = self$dpc_slope,
          verbose = FALSE
        )
      }

      private$.set_result(.elist_to_lfqdata(
        elist,
        wide,
        self$lfq$get_config(),
        self$lfq$data_long(),
        self$prefix,
        self$impute_only
      ))
    },

    #' @description
    #' creates aggregation plots (only for aggregation mode, not impute_only)
    #' @param subset create plots for a subset of the data only
    #' @param show.legend default FALSE
    #' @return data.frame
    plot = function(subset = NULL, show.legend = FALSE) {
      if (self$impute_only) {
        inform_prolfqua("AggregateLimpa$plot: aggregation plots not available in impute_only mode")
        return(invisible(NULL))
      }
      super$plot(subset = subset, show.legend = show.legend)
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
        inform_prolfqua("AggregateLimpa$write_plots: aggregation plots not available in impute_only mode")
        return(invisible(NULL))
      }
      super$write_plots(qcpath, subset, show.legend, width, height)
    }
  )
)


#' Protein-level LFQData of an estimate_intensity() result, with its estimation warnings in a new record
#' @param res result of \code{estimate_intensity()}
#' @param prefix output prefix
#' @param method aggregation method
#' @noRd
.aggregated_with_failures <- function(res, prefix, method) {
  result <- LFQData$new(res$data, res$config, prefix = prefix)
  log <- result$get_provenance_log()
  for (reason in unique(res$failures$warning)) {
    log$record(
      paste0("aggregate (", method, ")"),
      "subjects with estimation warning",
      res$failures$id[res$failures$warning == reason],
      reason,
      failure = TRUE,
      response = result$response()
    )
  }
  log$report()
  result
}
