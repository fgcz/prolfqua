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
#' # filter for missing values
#'
#' f1 <- lfqdata$omit_na(nr_na = 0)
#' stopifnot(f1$hierarchy_counts() <= lfqdata$hierarchy_counts())
#'
#' f2 <- lfqdata$omit_na(factor_depth = 0)
#' stopifnot(f2$hierarchy_counts() <= lfqdata$hierarchy_counts())
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
#' lfqdata2 <- lfqdata$get_copy()
#' lfqdata2$set_data(lfqdata2$data_long()[1:100, ])
#' res <- lfqdata$filter_difference(lfqdata2)
#' stopifnot(nrow(res$data_long()) == nrow(lfqdata$data_long()) - 100)
#'
#' tmp <- lfqdata$get_sample(5, seed = 4)
#' stopifnot(nrow(tmp$hierarchy()) == 5)
#'
LFQData <- R6::R6Class(
  "LFQData",
  private = list(
    .data = NULL,
    .config = NULL
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
      private$.data <- new_data
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
      if (!is.null(seed)) {
        set.seed(seed)
      }
      subset <- prolfqua::sample_subset(size = size, self$data_long(), self$relevant_hierarchy_keys())
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
        private$.config$is_response_transformed = is_transformed
      }
    },
    #' @description
    #' some software is reporting NA's as 0, you must remove it from your data
    #' @param threshold default 4.
    #' @return self
    remove_small_intensities = function(threshold = 4) {
      private$.data <- prolfqua::remove_small_intensities(
        self$data_long(),
        self$response(),
        threshold = threshold
      )
      invisible(self)
    },
    #' @description
    #' remove proteins with less than X peptides
    #' @return self
    filter_proteins_by_peptide_count = function() {
      message("removing proteins with less than: ", self$get_config()$min_peptides_protein, " peptpides")
      private$.data <- prolfqua::filter_proteins_by_peptide_count(self$data_long(), self$get_config())$data
      invisible(self)
    },
    #' @description
    #' Omit NA from intensities per hierarchy (e.g. protein or peptide), idea is to use it for normalization
    #' For instance if a peptide has a missing value in more then nrNA of the samples within a condition
    #' it will be removed
    #' @param nr_na number of NA values
    #' @param factor_depth control whether `nr_na` is applied per condition or more
    #'   globally, e.g. `factor_depth = 0` means per experiment
    #' @return LFQData with NA omitted.
    #'
    omit_na = function(nr_na = 0, factor_depth = NULL) {
      if (is.null(factor_depth)) {
        missing <- prolfqua::summarize_stats_factors(self)
      } else {
        if (factor_depth >= 1) {
          lfq_copy <- self$get_copy()
          lfq_copy$set_config_value("factor_depth", factor_depth)
          missing <- prolfqua::summarize_stats_factors(lfq_copy)
        } else {
          missing <- prolfqua::summarize_stats_all(self)
        }
      }
      not_na <- missing |> dplyr::filter(nrNAs <= nr_na)
      sum_n <- not_na |> group_by(across(all_of(self$hierarchy_keys()))) |> summarise(n = n())
      not_na <- sum_n |> dplyr::filter(n == max(n))

      not_na <- not_na |> dplyr::select(dplyr::all_of(self$hierarchy_keys()))
      not_na_data <- dplyr::inner_join(not_na, self$data_long()) |> ungroup()
      return(LFQData$new(not_na_data, self$get_config()))
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
      if (is.null(value)) {
        wide <- prolfqua::tidy_to_wide_config(self, as.matrix = as.matrix)
      } else {
        stopifnot(value %in% cfg$value_vars())
        wide <- prolfqua::tidy_to_wide_config(
          self,
          as.matrix = as.matrix,
          value = value
        )
      }
      wide$config <- cfg$clone(deep = TRUE)
      return(wide)
    },
    #' @description
    #' deprecated — use data_wide() instead
    #' @param as.matrix return as data.frame or matrix
    #' @param value see possible lfqdata$get_config()$value_vars()
    to_wide = function(as.matrix = FALSE, value = NULL) {
      .Deprecated("data_wide")
      self$data_wide(as.matrix = as.matrix, value = value)
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
    #' deprecated — use data_long() instead
    get_data = function() {
      .Deprecated("data_long")
      self$data_long()
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
        stop("Unknown aggregation method: ", method)
      )
    },
    #' @description
    #' get difference of self with other if other is subset of self
    #'
    #' @details
    #' Use to compare filtering results obtained from self, e.g. which proteins and peptides were removed (other)
    #'
    #' @param other a filtered LFQData set
    #' @return LFQData
    #'
    filter_difference = function(other) {
      diffdata <- prolfqua::filter_difference(
        self$data_long(),
        other$data_long(),
        self$get_config()
      )
      return(LFQData$new(diffdata, self$get_config()))
    }
  )
)

# Direct intensity manipulation ----

#' Remove rows when intensity lower then threshold
#' @param pdata data.frame
#' @param response character — name of the intensity column
#' @param threshold numeric — minimum intensity to keep (default 1)
#' @return data.frame
#' @export
#' @keywords internal
#' @family filtering
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config(Nprot = 20)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' res1 <- remove_small_intensities(lfqdata$data_long(), lfqdata$response(), threshold = 1)
#' res1000 <- remove_small_intensities(lfqdata$data_long(), lfqdata$response(), threshold = 1000)
#' stopifnot(nrow(res1) > nrow(res1000))
#'
remove_small_intensities <- function(pdata, response, threshold = 1) {
  res_data <- pdata |> dplyr::filter(!!sym(response) >= threshold)
  return(res_data)
}

# Hierarchy counting helpers ----

.make_name_AinB <- function(level_a, level_b, prefix = "nr_") {
  c_name <- paste(prefix, level_b, "_IN_", level_a, sep = "")
  return(c_name)
}

.nr_B_in_A <- function(data, level_a, level_b, merge = TRUE) {
  nam_a <- paste(level_a, collapse = "_")
  nam_b <- paste(level_b, collapse = "_")
  c_name <- .make_name_AinB(nam_a, nam_b)
  if (!c_name %in% colnames(data)) {
    data$c_name <- NULL
  }
  tmp <- data |>
    dplyr::select(all_of(c(level_a, level_b))) |>
    dplyr::distinct() |>
    dplyr::group_by(across(all_of(level_a))) |>
    dplyr::summarize(!!c_name := n())

  if (!merge) {
    return(tmp)
  }
  data <- dplyr::inner_join(data, tmp, by = level_a)
  message("Column added : ", c_name)
  return(list(data = data, name = c_name))
}


#' Compute nr of B per A
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @export
#' @keywords internal
#' @examples
#'
#' bb <- sim_lfq_data_peptide_config(Nprot = 100)
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#' hierarchy <- config$hierarchy_keys()
#' res <- nr_B_in_A(data, config)
#'
#' res$data |>
#'   dplyr::select(all_of(c(config$hierarchy_keys_depth(),  res$name))) |>
#'   dplyr::distinct() |>
#'   dplyr::pull() |> table()
#'
#'
#' bb <- prolfqua::prolfqua_data('data_skylineSRM_HL_A')
#' config <- bb$config_f()
#' data <- bb$data
#' data$Area[data$Area == 0] <- NA
#' analysis <- setup_analysis(data, config)
#'
#' resDataStart <- bb$analysis(bb$data, config)
#'
#' nr_B_in_A(resDataStart, config)
#' nr_B_in_A(resDataStart, config, merge = FALSE)
#' config$hierarchy_depth <- 2
#' nr_B_in_A(resDataStart, config, merge = FALSE)
#'
nr_B_in_A <- function(pdata, config, merge = TRUE) {
  level_a <- config$hierarchy_keys_depth()
  level_b <- config$hierarchy_keys()[length(level_a) + 1]
  if (is.na(level_b)) {
    warning("here is no B in A")
    return(NULL)
  } else {
    .nr_B_in_A(pdata, level_a, level_b, merge = merge)
  }
}

# Filtering helpers -----

#' Keep only those proteins with 2 IDENTIFIED peptides
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @return list with data.frame (data) and name of new column (name)
#' @export
#' @keywords internal
#' @family preprocessing
#' @examples
#'
#' istar <- prolfqua::sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' filterPep <- prolfqua::filter_proteins_by_peptide_count(istar$data, istar$config)
#' x <- prolfqua::summarize_hierarchy(filterPep$data,
#'   lfq$hierarchy_keys(), lfq$isotope_label())
#' stopifnot(x$peptide_Id_n >= istar$config$min_peptides_protein)
#'
filter_proteins_by_peptide_count <-
  function(pdata, config) {
    # remove single hit wonders
    tmp <- prolfqua::nr_B_in_A(pdata, config)
    if (!is.null(tmp)) {
      res <- dplyr::filter(tmp$data, !!sym(tmp$name) >= config$min_peptides_protein)
      name <- tmp$name
    } else {
      res <- pdata
      name <- NULL
    }
    return(list(data = res, name = name))
  }


#' get the difference of two dataset where one is a subset of the other.
#'
#' @param x data.frame
#' @param y data.frame
#' @param config AnlysisConfiguration
#' @return data.frame
#' @export
#' @keywords internal
#'
#' @examples
#'
#'
#' istar <- prolfqua::sim_lfq_data_peptide_config()
#' istar$config <- istar$config
#' istar_data <- istar$data
#' filterPep <- prolfqua:::filter_proteins_by_peptide_count( istar_data ,  istar$config )
#' tmp <- filter_difference(istar_data, filterPep$data, istar$config)
#' stopifnot(nrow(istar_data )  - nrow(filterPep$data) == nrow(tmp))
#' tmp <- filter_difference(filterPep$data, istar_data , istar$config)
#' stopifnot(nrow(istar_data )  - nrow(filterPep$data) == nrow(tmp))
#'
filter_difference <- function(x, y, config) {
  if (nrow(y) > nrow(x)) {
    dplyr::anti_join(y, x, by = config$id_vars())
  } else {
    dplyr::anti_join(x, y, by = config$id_vars())
  }
}


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
