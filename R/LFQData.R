#LFQData ----
#'
#' LFQData R6 class
#'
#' @section Missing Data Assumptions:
#' The filtering and imputation methods in this package assume that
#' missing values are Missing Completely At Random (MCAR) or Missing At
#' Random (MAR). Abundance-dependent missingness (MNAR), which is common
#' in DDA proteomics, is not modelled. Users should be aware that MNAR
#' can bias fold-change estimates and inflate false discovery rates.
#'
#' @return An R6 class generator.
#' @export
#' @family LFQData
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$filter_proteins_by_peptide_count()
#' tmp <- lfqdata$data_wide()
#' testthat::expect_equal(nrow(tmp$data) , nrow(tmp$rowdata))
#' testthat::expect_equal(ncol(tmp$data) , nrow(tmp$annotation) + ncol(tmp$rowdata))
#'
#' stopifnot("data.frame" %in% class(tmp$data))
#' tmp <- lfqdata$data_wide(as.matrix = TRUE)
#' stopifnot("matrix" %in% class(tmp$data))
#' stopifnot(lfqdata$is_transformed()==FALSE)
#' lfqdata$summarize_hierarchy()
#'
#' lfqdata$response()
#' lfqdata$rename_response("peptide.intensity")
#' lfqdata$response()
#' stopifnot("LFQData" %in% class(lfqdata$get_copy()))
#' stopifnot("LFQDataTransformer" %in% class(lfqdata$get_Transformer()))
#' stopifnot("LFQDataStats" %in% class(lfqdata$get_Stats()))
#' stopifnot("LFQDataSummariser" %in% class(lfqdata$get_Summariser()))
#' stopifnot("LFQDataPlotter" %in% class(lfqdata$get_Plotter()))
#' stopifnot("AggregateMedpolish" %in% class(lfqdata$get_Aggregator("medpolish")))
#'
#' tmp <-lfqdata$get_sample(5, seed = 4)
#' stopifnot(nrow(tmp$hierarchy()) == 5)
#'
LFQData <- R6::R6Class(
  "LFQData",
  private = list(
    .data = NULL,
    .config = NULL,
    # proportion of distinct modelling-level keys (subject_id) whose top-level
    # hierarchy id (protein_Id) is flagged by `detector` (is_decoy / is_contaminant).
    .prefix_proportion = function(detector, pattern) {
      top <- self$hierarchy_keys()[1]
      sid <- self$subject_id()
      keys <- self$data_long() |>
        dplyr::select(dplyr::all_of(unique(c(top, sid)))) |>
        dplyr::distinct()
      n_total <- nrow(dplyr::distinct(dplyr::select(keys, dplyr::all_of(sid))))
      if (n_total == 0) {
        return(0)
      }
      flagged <- detector(keys[[top]], pattern)
      n_flag <- nrow(dplyr::distinct(
        dplyr::select(keys[flagged, , drop = FALSE], dplyr::all_of(sid))
      ))
      n_flag / n_total
    },
    # Verify that the internal data still carries the columns implied by the
    # current configuration (structural join / grouping keys). Guards the
    # mutators, which would otherwise accept data that silently breaks
    # downstream operations.
    .validate_state = function() {
      cfg <- private$.config
      must_have <- unique(c(cfg$file_name, cfg$sample_name, cfg$factor_keys(), cfg$hierarchy_keys()))
      missing <- setdiff(must_have, names(private$.data))
      if (length(missing) > 0) {
        abort_missing_columns(missing, data_nm = "LFQData$data")
      }
      invisible(NULL)
    }
  ),
  public = list(
    #' @field prefix e.g. "peptide_", "protein_", "compound_"
    prefix = "",
    #' @description
    #' initialize
    #' @param data data.frame
    #' @param config configuration
    #' @param is_pep todo
    #' @param prefix will be use as output prefix
    #' @param setup is data setup needed, default = FALSE, if TRUE, calls \code{\link{setup_analysis}} on data first.
    initialize = function(data, config, prefix = "ms_", setup = FALSE) {
      private$.data <- if (setup) {
        setup_analysis(data, config)
      } else {
        data
      }
      private$.config <- config$clone(deep = TRUE)
      self$prefix <- prefix
    },
    #' @description
    #' set data (replaces the internal data frame)
    #' @param new_data data.frame
    #' @return self (invisible)
    set_data = function(new_data) {
      if (!is.data.frame(new_data)) {
        abort_bad_argument("new_data", "be a data frame", not = paste(class(new_data), collapse = "/"))
      }
      private$.data <- new_data
      private$.validate_state()
      invisible(self)
    },
    #' @description
    #' return the AnalysisConfiguration object
    #' @return AnalysisConfiguration
    get_config = function() {
      private$.config
    },
    #' @description
    #' set a config field value
    #' @param field character — field name
    #' @param value the value to set
    set_config_value = function(field, value) {
      private$.config[[field]] <- value
      private$.validate_state()
      invisible(self)
    },
    #' @description
    #' get deep copy
    get_copy = function() {
      return(self$clone(deep = TRUE))
    },
    #' @description
    #' samples subset of data
    #' @param size size of subset default 100
    #' @param seed set seed
    get_sample = function(size = 100, seed = NULL) {
      sample_data <- function() {
        prolfqua::sample_subset(size = size, self$data_long(), self$relevant_hierarchy_keys())
      }
      subset <- if (is.null(seed)) {
        sample_data()
      } else {
        withr::with_seed(seed, sample_data())
      }
      return(LFQData$new(subset, self$get_config()))
    },
    #' @description
    #' get subset of data
    #' @param x data frame with columns containing subject_id
    get_subset = function(x) {
      x <- select(x, any_of(self$subject_id())) |> distinct()
      subset <- inner_join(x, self$data_long())
      return(LFQData$new(subset, self$get_config()))
    },
    #' @description
    #' get subject ID columns
    subject_id = function() {
      return(private$.config$hierarchy_keys_depth())
    },
    #' @description
    #' is data transformed
    #' @param is_transformed logical
    #' @return logical
    is_transformed = function(is_transformed) {
      if (missing(is_transformed)) {
        return(private$.config$is_response_transformed)
      } else {
        private$.config$is_response_transformed <- is_transformed
      }
    },
    #' @description
    #' some software is reporting NA's as 0, you must remove it from your data
    #' @param threshold default 4.
    #' @return self
    remove_small_intensities = function(threshold = 4) {
      private$.data <- dplyr::filter(self$data_long(), !!sym(self$response()) >= threshold)
      private$.data <- prolfqua::complete_cases(self)
      invisible(self)
    },
    #' @description
    #' remove proteins with less than X peptides
    #' @return self
    filter_proteins_by_peptide_count = function() {
      cfg <- self$get_config()
      message("removing proteins with less than: ", cfg$min_peptides_protein, " peptpides")
      level_a <- cfg$hierarchy_keys_depth()
      level_b <- cfg$hierarchy_keys()[length(level_a) + 1]
      if (is.na(level_b)) {
        warning("here is no B in A")
        return(invisible(self))
      }
      c_name <- paste0("nr_", level_b, "_IN_", paste(level_a, collapse = "_"))
      counts <- dplyr::distinct(dplyr::select(self$data_long(), dplyr::all_of(c(level_a, level_b)))) |>
        dplyr::group_by(dplyr::across(dplyr::all_of(level_a))) |>
        dplyr::summarize(!!c_name := dplyr::n())
      data <- dplyr::inner_join(self$data_long(), counts, by = level_a)
      message("Column added : ", c_name)
      private$.data <- dplyr::filter(data, !!sym(c_name) >= cfg$min_peptides_protein)
      invisible(self)
    },
    #' @description
    #' remove decoy / reverse-database entries (rows whose top-level hierarchy id,
    #' e.g. protein_Id, is a decoy). Detection uses the configured
    #' `pattern_decoys` unioned with the built-in defaults (see
    #' \code{\link{is_decoy}}); decoys are always detectable even when no pattern
    #' is configured. Returns a new decoy-free LFQData (self is not modified).
    #' @return LFQData without decoys
    remove_decoys = function() {
      top <- self$hierarchy_keys()[1]
      pdata <- self$data_long()
      keep <- !prolfqua::is_decoy(pdata[[top]], self$get_config()$pattern_decoys)
      return(LFQData$new(pdata[keep, , drop = FALSE], self$get_config()))
    },
    #' @description
    #' proportion of modelling-level keys (subject_id, e.g. protein or peptide by
    #' `hierarchy_depth`) that are decoys — an empirical-FDR / target-decoy QC
    #' signal. Decoy status derives from the top-level hierarchy id (protein_Id).
    #' @return numeric in [0, 1]
    decoy_proportion = function() {
      private$.prefix_proportion(prolfqua::is_decoy, self$get_config()$pattern_decoys)
    },
    #'
    #' @description
    #' some software is reporting NA's as 0, you must remove it from your data
    #' @param threshold default 4.
    #' @return self
    complete_cases = function() {
      private$.data <- prolfqua::complete_cases(self)
      invisible(self)
    },
    #' @description
    #' converts the data to wide
    #' @param as.matrix return as data.frame or matrix
    #' @param value see possible lfqdata$get_config()$value_vars()
    #' @return list with data, annotation, and configuration
    #'
    data_wide = function(as.matrix = FALSE, value = NULL) {
      cfg <- self$get_config()
      if (!is.null(value)) {
        stopifnot(value %in% cfg$value_vars())
      }
      wide <- prolfqua::tidy_to_wide_config(self, as.matrix = as.matrix, value = value %||% self$response())
      wide$config <- cfg$clone(deep = TRUE)
      return(wide)
    },
    #' @description
    #' Annotation table
    #' @return data.frame
    factors = function() {
      prolfqua::table_factors(self$data_long(), self$file_name(), self$sample_name(), self$factor_keys())
    },
    #' @description
    #' Hierarchy table
    hierarchy = function() {
      hk <- self$relevant_hierarchy_keys()
      hkdf <- self$data_long() |> select(all_of(hk)) |> distinct()
      return(hkdf)
    },
    #' @description
    #' name of response variable
    #' @return data.frame
    response = function() {
      private$.config$get_response()
    },
    #' @description
    #' return all hierarchy column names
    hierarchy_keys = function() {
      private$.config$hierarchy_keys()
    },
    #' @description
    #' return hierarchy column names at current depth (alias for subject_id)
    relevant_hierarchy_keys = function() {
      private$.config$hierarchy_keys_depth()
    },
    #' @description
    #' return all factor column names
    factor_keys = function() {
      private$.config$factor_keys()
    },
    #' @description
    #' return factor column names at current depth
    relevant_factor_keys = function() {
      private$.config$factor_keys_depth()
    },
    #' @description
    #' return sample name column
    sample_name = function() {
      private$.config$sample_name
    },
    #' @description
    #' return file name column
    file_name = function() {
      private$.config$file_name
    },
    #' @description
    #' return name of nr_children column
    nr_children_col = function() {
      private$.config$nr_children
    },
    #' @description
    #' return isotope label column name
    isotope_label = function() {
      private$.config$isotope_label
    },
    #' @description
    #' return the tidy (long-format) data frame
    #' @param na.omit if TRUE, remove rows with NA in response column
    data_long = function(na.omit = FALSE) {
      if (na.omit) {
        return(stats::na.omit(private$.data))
      }
      private$.data
    },
    #' @description
    #' new name of response variable
    #' @param newname default Intensity
    rename_response = function(newname = "Intensity") {
      if ((newname %in% colnames(self$data_long()))) {
        msg <- paste(newname, " already in data :", paste(colnames(self$data_long()), collapse = " "), ".")
        message(msg)
      } else {
        cfg <- self$get_config()
        old <- cfg$pop_response()
        cfg$set_response(newname)
        private$.data <- self$data_long() |> dplyr::rename(!!newname := !!sym(old))
      }
    },
    #' @description
    #' number of elements at each level
    hierarchy_counts = function() {
      prolfqua::hierarchy_counts(self$data_long(), self$hierarchy_keys(), self$isotope_label())
    },
    #' @description
    #' e.g. number of peptides per protein etc
    #' @return data.frame
    summarize_hierarchy = function() {
      prolfqua::summarize_hierarchy(
        self$data_long(),
        self$hierarchy_keys(),
        self$isotope_label(),
        hierarchy = self$relevant_hierarchy_keys()
      )
    },
    #' @description
    #' get Plotter
    #' @return LFQDataPlotter
    get_Plotter = function() {
      return(LFQDataPlotter$new(self, self$prefix))
    },
    #' @description
    #' get Summariser
    #' @return LFQDataSummarizer
    get_Summariser = function() {
      return(LFQDataSummariser$new(self))
    },
    #' @description
    #' Get \code{\link{LFQDataStats}}. For more details see \code{\link{LFQDataStats}}.
    #' @param stats default interaction, computes statistics within interaction.
    #' @return LFQDataStats
    get_Stats = function(stats = c("everything", "interaction", "all")) {
      stats <- match.arg(stats)
      return(LFQDataStats$new(self, stats = stats))
    },
    #' @description
    #' get Stats
    #' @return LFQDataTransformer
    get_Transformer = function() {
      return(LFQDataTransformer$new(self))
    },
    #' @description
    #' get Aggregator
    #' @param method aggregation method: "medpolish", "rlm", or "topN"
    #' @param ... passed to aggregator constructor (e.g. prefix, N, func)
    #' @return AggregateMedpolish, AggregateRlm, or AggregateTopN
    get_Aggregator = function(method = "medpolish", ...) {
      switch(
        method,
        "medpolish" = AggregateMedpolish$new(self, ...),
        "rlm" = AggregateRlm$new(self, ...),
        "topN" = AggregateTopN$new(self, ...),
        abort_bad_argument("method", 'be one of: "medpolish", "rlm", "topN"', not = method)
      )
    }
  )
)

#' converts LFQData object to SummarizedExperiment
#'
#' For compatibility with Bioconductor
#' @param lfqdata LFQData object
#' @return SummarizedExperiment (bioconductor)
#' @family LFQData
#' @export
#' @examples
#'
#' istar <- prolfqua::sim_lfq_data_peptide_config()
#' data <- istar$data
#' lfqdata <- LFQData$new(data, istar$config)
#' lfqdata$data_wide()
#' if(require("SummarizedExperiment")){
#'    tmp <- LFQDataToSummarizedExperiment(lfqdata)
#' }
#'
LFQDataToSummarizedExperiment <- function(lfqdata) {
  if (requireNamespace("SummarizedExperiment")) {
    wide <- lfqdata$data_wide(as.matrix = TRUE)
    ann <- data.frame(wide$annotation)
    rownames(ann) <- wide$annotation[[lfqdata$sample_name()]]
    assays <- S4Vectors::SimpleList(LFQ = wide$data)
    if ("nr_children" %in% lfqdata$get_config()$value_vars()) {
      nr_children <- lfqdata$data_wide(as.matrix = TRUE, value = "nr_children")
      assays[["nr_children"]] <- nr_children$data
    }
    se <- SummarizedExperiment::SummarizedExperiment(
      assays,
      colData = ann,
      rowData = wide$rowdata
    )
    return(se)
  }
}
