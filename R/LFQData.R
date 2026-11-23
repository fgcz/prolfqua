#LFQData ----
#'
#' LFQData R6 class
#'
#' The data container of prolfqua: one table in long format, one row per feature (for example peptide) and sample,
#' plus an [AnalysisConfiguration] that names the hierarchy, factor, sample and response columns. Its decorators
#' (`get_Transformer()`, `get_Aggregator()`, `get_Stats()`, `get_Plotter()`, `get_Summariser()`) transform, aggregate,
#' summarise and plot the data, and [build_contrast_analysis()] fits models to it.
#'
#' @section Missing Data Assumptions:
#' The filtering and imputation methods in this package assume that
#' missing values are Missing Completely At Random (MCAR) or Missing At
#' Random (MAR). Abundance-dependent missingness (MNAR), which is common
#' in DDA proteomics, is not modelled. When low intensities are missing,
#' fold changes of down-regulated proteins shrink toward zero; LOD imputation
#' (the \code{_impute} methods) caps estimates at the LOD and does not remove
#' this bias. MNAR can also inflate false discovery rates. How \code{NA},
#' zero and non-finite intensities are treated is described in
#' \code{\link{prolfqua-package}}, section "Missing and non-finite values".
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
#' stopifnot(nrow(tmp$data) == nrow(tmp$rowdata))
#' stopifnot(ncol(tmp$data) == nrow(tmp$annotation) + ncol(tmp$rowdata))
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
#' @srrstats {G2.8} LFQData wraps validated numeric measurements and their AnalysisConfiguration for all facades.
#' @srrstats {RE4.13} LFQData exposes categorical predictor values, factor keys and configuration.
#' @srrstats {RE7.2} Tests preserve subject IDs and sample alignment, including wide-to-long roundtrips.
LFQData <- R6::R6Class(
  "LFQData",
  private = list(
    .data = NULL,
    .config = NULL,
    .provenance = NULL,
    # new LFQData without the dropped rows, continuing this record with the removed top-level ids
    .remove_rows = function(pdata, keep, step, reason) {
      top <- self$hierarchy_keys()[1]
      result <- LFQData$new(pdata[keep, , drop = FALSE], self$get_config(), provenance = private$.provenance)
      log <- result$get_provenance_log()
      log$record(step, paste0(top, " removed"), pdata[[top]][!keep], reason, response = self$response())
      log$report()
      result
    },
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
    #' @param data long-format data.frame, tibble or data.table; stored as a tibble
    #' @param config AnalysisConfiguration naming the column roles
    #' @param prefix one character string, used as output prefix
    #' @param setup one TRUE/FALSE, default FALSE; TRUE calls \code{\link{setup_analysis}} on data first
    #' @param provenance internal: record of an LFQData this one continues (the same data transformed or
    #'   subset); a copy is kept. NULL starts a new record.
    initialize = function(data, config, prefix = "ms_", setup = FALSE, provenance = NULL) {
      .assert_response_data(data, config)
      .assert_string(prefix)
      .assert_flag(setup)
      # one class for every input (data.frame, tibble, data.table), so column extraction behaves the same
      data <- tibble::as_tibble(data)
      private$.data <- if (setup) {
        setup_analysis(data, config)
      } else {
        keys <- c(config$file_name, config$hierarchy_keys())
        missing_keys <- setdiff(keys, names(data))
        if (length(missing_keys) > 0L) {
          abort_missing_columns(missing_keys)
        }
        .assert_unique_keys(data, keys)
        data
      }
      private$.data <- .drop_response_attributes(private$.data, config$get_response())
      private$.config <- config$clone(deep = TRUE)
      self$prefix <- prefix
      private$.provenance <- if (is.null(provenance)) Provenance$new() else provenance$clone(deep = TRUE)
    },
    #' @description
    #' problems and dropped items recorded for this data: values set to NA, rows, features or proteins removed
    #' @return tibble with columns step, event, n, ids (list), reason, failure, response
    provenance = function() {
      private$.provenance$entries()
    },
    #' @description
    #' internal: the record this data continues, passed to \code{LFQData$new()} for derived data
    #' @return Provenance object
    get_provenance_log = function() {
      private$.provenance
    },
    #' @description
    #' print samples, response, hierarchy counts, samples per factor level, missing responses and provenance
    #' @param digits significant digits of printed numbers (default 3)
    #' @param ... unused
    #' @return self (invisible)
    print = function(digits = 3, ...) {
      data <- self$data_long()
      response <- self$response()
      samples <- self$factors()
      keys <- self$hierarchy_keys()
      levels <- vapply(seq_along(keys), function(i) dplyr::n_distinct(data[keys[seq_len(i)]]), integer(1))
      n_missing <- sum(is.na(data[[response]]))
      range <- range(data[[response]], na.rm = TRUE)
      cat("LFQData:", nrow(samples), "samples\n")
      cat("Response:", response, if (self$is_transformed()) "(transformed)" else "(not transformed)", "\n")
      cat("Hierarchy:", paste(keys, levels, collapse = ", "), "\n")
      for (key in self$relevant_factor_keys()) {
        counts <- table(samples[[key]])
        cat("Samples per ", key, ": ", paste(names(counts), counts, collapse = ", "), "\n", sep = "")
      }
      cat(
        "Response range: ",
        paste(trimws(formatC(range, digits = digits, format = "g")), collapse = " to "),
        "; missing ",
        n_missing,
        " of ",
        nrow(data),
        " (",
        trimws(formatC(100 * n_missing / nrow(data), digits = digits, format = "g")),
        "%)\n",
        sep = ""
      )
      cat(.provenance_line(self$provenance()), "\n")
      invisible(self)
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
    #' @return the AnalysisConfiguration of this object itself, not a copy: changing a field of it changes this
    #'   `LFQData`. Use `get_config()$clone(deep = TRUE)` to modify a configuration for another object.
    get_config = function() {
      private$.config
    },
    #' @description
    #' set a config field value
    #' @param field character - field name
    #' @param value the value to set
    set_config_value = function(field, value) {
      private$.config[[field]] <- value
      private$.validate_state()
      invisible(self)
    },
    #' @description
    #' get deep copy
    get_copy = function() {
      self$clone(deep = TRUE)
    },
    #' @description
    #' can the data be aggregated, i.e. is there a hierarchy level below the
    #' modelled one (`hierarchy_depth`)
    #' @return logical
    is_aggregatable = function() {
      n_levels <- length(self$hierarchy_keys())
      n_levels > 1 && length(self$relevant_hierarchy_keys()) < n_levels
    },
    #' @description
    #' samples subset of data
    #' @param size size of subset default 100
    #' @param seed set seed
    get_sample = function(size = 100, seed = NULL) {
      sample_data <- function() {
        .sample_subset(size = size, self$data_long(), self$relevant_hierarchy_keys())
      }
      subset <- if (is.null(seed)) {
        sample_data()
      } else {
        withr::with_seed(seed, sample_data())
      }
      LFQData$new(subset, self$get_config(), provenance = private$.provenance)
    },
    #' @description
    #' get subset of data
    #' @param x data frame with columns containing subject_id
    get_subset = function(x) {
      x <- select(x, any_of(self$subject_id())) |> distinct()
      if (ncol(x) == 0L) {
        abort_missing_columns(self$subject_id(), data_nm = "x")
      }
      subset <- inner_join(x, self$data_long(), by = colnames(x), relationship = "one-to-many")
      LFQData$new(subset, self$get_config(), provenance = private$.provenance)
    },
    #' @description
    #' get subject ID columns
    subject_id = function() {
      private$.config$hierarchy_keys_depth()
    },
    #' @description
    #' is data transformed
    #' @param is_transformed logical
    #' @return logical when called without argument, otherwise self (invisible)
    is_transformed = function(is_transformed) {
      if (missing(is_transformed)) {
        return(private$.config$is_response_transformed)
      }
      private$.config$is_response_transformed <- is_transformed
      invisible(self)
    },
    #' @description
    #' some software is reporting NA's as 0, you must remove it from your data
    #' @param threshold default 4.
    #' @return self
    remove_small_intensities = function(threshold = 4) {
      response <- self$response()
      before <- self$data_long()
      small <- !is.na(before[[response]]) & before[[response]] < threshold
      private$.data <- dplyr::filter(before, .data[[response]] >= .env$threshold)
      features <- self$hierarchy_keys()
      removed <- setdiff(.key_ids(before, features), .key_ids(private$.data, features))
      private$.data <- .complete_cases(self)
      step <- "remove_small_intensities"
      reason <- paste0("intensity below ", threshold)
      private$.provenance$record(
        step,
        "values set to NA",
        .key_ids(before[small, ], features),
        reason,
        n = sum(small),
        response = response
      )
      private$.provenance$record(step, "features removed", removed, reason, response = response)
      private$.provenance$report()
      invisible(self)
    },
    #' @description
    #' remove parents (e.g. proteins) with fewer than `threshold` distinct
    #' children (e.g. peptides) and add the count as a column. By default the
    #' second hierarchy level is counted within the first, independent of
    #' `hierarchy_depth`; a hierarchy with one level is an error.
    #' @param threshold minimum number of children, default `config$min_peptides_protein`
    #' @param parent parent key column, default the first hierarchy key
    #' @param child child key column, default the second hierarchy key
    #' @return self (invisible)
    filter_proteins_by_peptide_count = function(threshold = NULL, parent = NULL, child = NULL) {
      cfg <- self$get_config()
      keys <- cfg$hierarchy_keys()
      parent <- parent %||% keys[1]
      child <- child %||% keys[2]
      threshold <- threshold %||% cfg$min_peptides_protein
      before <- self$data_long()
      filtered <- filter_proteins_by_peptide_count(
        before,
        cfg,
        parent = parent,
        child = child,
        threshold = threshold
      )
      inform_prolfqua("filter_proteins_by_peptide_count: column added: ", filtered$name)
      private$.data <- filtered$data
      private$.provenance$record(
        "filter_proteins_by_peptide_count",
        paste0(parent, " removed"),
        setdiff(unique(before[[parent]]), unique(filtered$data[[parent]])),
        paste0("fewer than ", threshold, " ", child),
        response = self$response()
      )
      private$.provenance$report()
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
      private$.remove_rows(pdata, keep, "remove_decoys", "decoy")
    },
    #' @description
    #' remove contaminant entries (rows whose top-level hierarchy id, e.g.
    #' protein_Id, is a contaminant). Detection uses the configured
    #' `pattern_contaminants` unioned with the built-in defaults (see
    #' \code{\link{is_contaminant}}). Returns a new contaminant-free LFQData
    #' (self is not modified).
    #' @return LFQData without contaminants
    remove_contaminants = function() {
      top <- self$hierarchy_keys()[1]
      pdata <- self$data_long()
      keep <- !prolfqua::is_contaminant(pdata[[top]], self$get_config()$pattern_contaminants)
      private$.remove_rows(pdata, keep, "remove_contaminants", "contaminant")
    },
    #' @description
    #' proportion of modelling-level keys (subject_id, e.g. protein or peptide by
    #' `hierarchy_depth`) that are decoys - an empirical-FDR / target-decoy QC
    #' signal. Decoy status derives from the top-level hierarchy id (protein_Id).
    #' @return numeric in \[0, 1\]
    decoy_proportion = function() {
      private$.prefix_proportion(prolfqua::is_decoy, self$get_config()$pattern_decoys)
    },
    #'
    #' @description
    #' some software is reporting NA's as 0, you must remove it from your data
    #' @param threshold default 4.
    #' @return self
    complete_cases = function() {
      private$.data <- .complete_cases(self)
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
        if (!all(value %in% cfg$value_vars())) {
          abort_bad_argument("value", paste0("be one of ", toString(cfg$value_vars())), not = value)
        }
      }
      wide <- .tidy_to_wide_config(self, as.matrix = as.matrix, value = value %||% self$response())
      wide$config <- cfg$clone(deep = TRUE)
      wide
    },
    #' @description
    #' Annotation table
    #' @return data.frame
    factors = function() {
      .table_factors(self$data_long(), self$file_name(), self$sample_name(), self$factor_keys())
    },
    #' @description
    #' Hierarchy table
    hierarchy = function() {
      hk <- self$relevant_hierarchy_keys()
      hkdf <- self$data_long() |> select(all_of(hk)) |> distinct()
      hkdf
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
    #' return the tidy (long-format) data frame
    #' @param na.omit if TRUE, remove rows with NA in response column
    data_long = function(na.omit = FALSE) {
      if (na.omit) {
        return(stats::na.omit(private$.data))
      }
      private$.data
    },
    #' @description
    #' new name of response variable. Renaming to the current response name
    #' changes nothing; renaming to another existing column is an error.
    #' @param newname default Intensity
    #' @return self (invisible)
    rename_response = function(newname = "Intensity") {
      if (identical(newname, self$response())) {
        return(invisible(self))
      }
      if (newname %in% colnames(self$data_long())) {
        abort_bad_argument("newname", "not name an existing column", not = newname)
      }
      cfg <- self$get_config()
      old <- cfg$pop_response()
      cfg$set_response(newname)
      private$.data <- self$data_long() |> dplyr::rename(!!newname := dplyr::all_of(old))
      invisible(self)
    },
    #' @description
    #' number of elements at each level
    hierarchy_counts = function() {
      .hierarchy_counts(self$data_long(), self$hierarchy_keys())
    },
    #' @description
    #' e.g. number of peptides per protein etc
    #' @return data.frame
    summarize_hierarchy = function() {
      .summarize_hierarchy(
        self$data_long(),
        self$hierarchy_keys(),
        hierarchy = self$relevant_hierarchy_keys()
      )
    },
    #' @description
    #' get Plotter
    #' @return LFQDataPlotter
    get_Plotter = function() {
      LFQDataPlotter$new(self, self$prefix)
    },
    #' @description
    #' get Summariser
    #' @return LFQDataSummariser
    get_Summariser = function() {
      LFQDataSummariser$new(self)
    },
    #' @description
    #' Get \code{\link{LFQDataStats}}. For more details see \code{\link{LFQDataStats}}.
    #' @param stats default interaction, computes statistics within interaction.
    #' @return LFQDataStats
    get_Stats = function(stats = c("everything", "interaction", "all")) {
      stats <- match.arg(stats)
      LFQDataStats$new(self, stats = stats)
    },
    #' @description
    #' get Stats
    #' @return LFQDataTransformer
    get_Transformer = function() {
      LFQDataTransformer$new(self)
    },
    #' @description
    #' get Aggregator
    #' @param method aggregation method: "medpolish", "rlm", "topN", or "limpa"
    #' @param ... passed to aggregator constructor (e.g. prefix, N, func, dpc_slope)
    #' @return AggregateMedpolish, AggregateRlm, AggregateTopN, or AggregateLimpa
    get_Aggregator = function(method = "medpolish", ...) {
      method <- .arg_match_prolfqua(method, c("medpolish", "rlm", "topN", "limpa"), "method")
      switch(
        method,
        "medpolish" = AggregateMedpolish$new(self, ...),
        "rlm" = AggregateRlm$new(self, ...),
        "topN" = AggregateTopN$new(self, ...),
        "limpa" = AggregateLimpa$new(self, ...)
      )
    }
  )
)

#' Store the response as a plain double: drop extra classes and attributes (AsIs, units, labels), convert integers
#' @param data long-format data
#' @param response response column
#' @noRd
#' @srrstats {EA2.6} A response with extra classes or attributes (e.g. `I()`, units, labelled) becomes a plain
#'   numeric vector at `LFQData$new()`, with a message, so transformations and models see the values only.
#' @srrstats {G2.4b} The response is converted with `as.numeric()` to double, also when it is an integer count, in
#'   `LFQData$new()` and `setup_analysis()` (test-input-types.R).
.drop_response_attributes <- function(data, response) {
  values <- data[[response]]
  if (is.null(attributes(values))) {
    if (is.integer(values)) {
      data[[response]] <- as.numeric(values)
    }
    return(data)
  }
  inform_prolfqua(
    "Response `",
    response,
    "`: dropped attributes ",
    toString(names(attributes(values))),
    "; prolfqua uses the numeric values only."
  )
  data[[response]] <- as.numeric(values)
  data
}

#' Plot an LFQData object
#'
#' Default plot: distribution of the response per sample (violin plot), on a log scale when the response is not
#' transformed. Other plots are available from `x$get_Plotter()`.
#' @param x LFQData object
#' @param ... unused
#' @return ggplot
#' @family LFQData
#' @exportS3Method graphics::plot
#' @srrstats {EA4.2} LFQData has default `print()` (samples, response, hierarchy, factor levels, missing values,
#'   provenance) and `plot()` (response distribution per sample); fitted analyses have `print()`, `summary()` and
#'   `plot()`.
#' @examples
#' istar <- sim_lfq_data_peptide_config()
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata
#' plot(lfqdata)
plot.LFQData <- function(x, ...) {
  x$get_Plotter()$intensity_distribution_violin()
}

#' Convert an LFQData object to a SummarizedExperiment
#'
#' Builds a `SummarizedExperiment` with the response as the assay `LFQ` (and the child counts as `nr_children`
#' when present; features in rows, samples in columns), the hierarchy keys as row data and the sample annotation as
#' column data, for use with Bioconductor packages.
#'
#' For compatibility with Bioconductor
#' @param lfqdata [LFQData] object with the data and its configuration
#' @return SummarizedExperiment (bioconductor)
#' @family LFQData
#' @export
#' @examples
#'
#' istar <- prolfqua::sim_lfq_data_peptide_config()
#' data <- istar$data
#' lfqdata <- LFQData$new(data, istar$config)
#' lfqdata$data_wide()
#' if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'    tmp <- lfqdata_to_summarized_experiment(lfqdata)
#' }
#'
lfqdata_to_summarized_experiment <- function(lfqdata) {
  rlang::check_installed(c("SummarizedExperiment", "S4Vectors"), reason = "to build a SummarizedExperiment")
  wide <- lfqdata$data_wide(as.matrix = TRUE)
  ann <- data.frame(wide$annotation)
  rownames(ann) <- wide$annotation[[lfqdata$sample_name()]]
  assays <- S4Vectors::SimpleList(LFQ = wide$data)
  if ("nr_children" %in% lfqdata$get_config()$value_vars()) {
    nr_children <- lfqdata$data_wide(as.matrix = TRUE, value = "nr_children")
    assays[["nr_children"]] <- nr_children$data
  }
  SummarizedExperiment::SummarizedExperiment(
    assays,
    colData = ann,
    rowData = wide$rowdata
  )
}
