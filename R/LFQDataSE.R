# LFQDataSE ----
#'
#' LFQDataSE R6 class - SummarizedExperiment-backed alternative to LFQData
#'
#' This class provides the same interface as LFQData but uses SummarizedExperiment
#' as the internal storage backend. This enables better Bioconductor interoperability
#' and improved performance for matrix operations.
#'
#' @export
#' @family LFQData
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' lfqdata <- LFQData$new(istar$data, istar$config)
#'
#' # Convert to LFQDataSE
#' lfq_se <- LFQDataSE$new(data = istar$data, config = istar$config)
#'
#' # Or from existing LFQData
#' lfq_se2 <- as_LFQDataSE(lfqdata)
#'
#' # Access underlying SummarizedExperiment
#' se <- lfq_se$se()
#' SummarizedExperiment::assay(se)[1:5, 1:3]
#'
#' # Use like LFQData (same interface)
#' lfq_se$to_wide()
#' lfq_se$hierarchy_counts()
#'
#' # Convert back to LFQData
#' lfq_back <- as_LFQData(lfq_se)
#'
LFQDataSE <- R6::R6Class(
  "LFQDataSE",
  public = list(
    #' @field config AnalysisConfiguration
    config = NULL,
    #' @field .se SummarizedExperiment object (internal storage)
    .se = NULL,
    #' @field prefix e.g. "peptide_", "protein_", "compound_"
    prefix = "",
    #' @field .data_cache cached long-format data
    .data_cache = NULL,
    #' @field .data_dirty flag indicating if cache needs refresh
    .data_dirty = TRUE,

    #' @description
    #' initialize LFQDataSE
    #' @param se SummarizedExperiment object (optional, provide either se or data+config)
    #' @param data data.frame (optional, used if se is NULL)
    #' @param config AnalysisConfiguration (required if using data)
    #' @param prefix output prefix
    #' @param setup is data setup needed
    initialize = function(se = NULL, data = NULL, config = NULL, prefix = "ms_", setup = FALSE) {
      if (!is.null(se)) {
        # Initialize from SummarizedExperiment
        self$.se <- se
        if (!is.null(config)) {
          self$config <- config$clone(deep = TRUE)
        } else {
          # Try to get config from SE metadata
          self$config <- S4Vectors::metadata(se)$config
        }
        self$prefix <- if (!is.null(S4Vectors::metadata(se)$prefix)) {
          S4Vectors::metadata(se)$prefix
        } else {
          prefix
        }
      } else if (!is.null(data) && !is.null(config)) {
        # Initialize from data + config (like LFQData)
        if (setup) {
          data <- setup_analysis(data, config)
        }
        self$config <- config$clone(deep = TRUE)
        self$prefix <- prefix
        # Convert to SE
        self$.se <- private$long_to_se(data)
      } else {
        stop("Must provide either 'se' or both 'data' and 'config'")
      }
      self$.data_dirty <- TRUE
    },

    #' @description
    #' Access underlying SummarizedExperiment
    #' @return SummarizedExperiment object
    se = function() {
      return(self$.se)
    },

    #' @description
    #' Access assay data from SE
    #' @param name assay name (default: first assay)
    #' @return matrix
    assay = function(name = NULL) {
      if (is.null(name)) {
        return(SummarizedExperiment::assay(self$.se))
      } else {
        return(SummarizedExperiment::assay(self$.se, name))
      }
    },

    #' @description
    #' Access all assays from SE
    #' @return SimpleList of assays
    assays = function() {
      return(SummarizedExperiment::assays(self$.se))
    },

    #' @description
    #' Access column (sample) annotations from SE
    #' @return DataFrame
    colData = function() {
      return(SummarizedExperiment::colData(self$.se))
    },

    #' @description
    #' Access row (feature) annotations from SE
    #' @return DataFrame
    rowData = function() {
      return(SummarizedExperiment::rowData(self$.se))
    },

    #' @description
    #' get deep copy
    get_copy = function() {
      new_se <- self$.se
      return(LFQDataSE$new(se = new_se, config = self$config, prefix = self$prefix))
    },

    #' @description
    #' samples subset of data
    #' @param size size of subset default 100
    #' @param seed set seed
    get_sample = function(size = 100, seed = NULL) {
      if (!is.null(seed)) { set.seed(seed) }
      subset <- prolfqua::sample_subset(size = size, self$data, self$config)
      return(LFQDataSE$new(data = subset, config = self$config$clone(deep = TRUE), prefix = self$prefix))
    },

    #' @description
    #' get subset of data
    #' @param x data frame with columns containing subject_Id
    get_subset = function(x) {
      x <- dplyr::select(x, dplyr::any_of(self$subject_Id())) |> dplyr::distinct()
      subset <- dplyr::inner_join(x, self$data)
      return(LFQDataSE$new(data = subset, config = self$config$clone(deep = TRUE), prefix = self$prefix))
    },

    #' @description
    #' get subject ID columns
    subject_Id = function() {
      return(self$config$table$hierarchy_keys_depth())
    },

    #' @description
    #' is data transformed
    #' @param is_transformed logical
    #' @return logical
    is_transformed = function(is_transformed) {
      if (missing(is_transformed)) {
        return(self$config$table$is_response_transformed)
      } else {
        self$config$table$is_response_transformed <- is_transformed
      }
    },

    #' @description
    #' some software is reporting NA's as 0, you must remove it from your data
    #' @param threshold default 4.
    #' @return self
    remove_small_intensities = function(threshold = 4) {
      data <- prolfqua::remove_small_intensities(self$data, self$config, threshold = threshold)
      self$.se <- private$long_to_se(data)
      self$.data_dirty <- TRUE
      invisible(self)
    },

    #' @description
    #' remove proteins with less than X peptides
    #' @return self
    filter_proteins_by_peptide_count = function() {
      message("removing proteins with less than: ",
              self$config$parameter$min_peptides_protein,
              " peptides")
      data <- prolfqua::filter_proteins_by_peptide_count(self$data, self$config)$data
      self$.se <- private$long_to_se(data)
      self$.data_dirty <- TRUE
      invisible(self)
    },

    #' @description
    #' Omit NA from intensities per hierarchy
    #' @param nrNA number of NA values
    #' @param factorDepth factor depth for grouping
    #' @return LFQDataSE with NA omitted
    omit_NA = function(nrNA = 0, factorDepth = NULL) {
      if (is.null(factorDepth)) {
        missing <- prolfqua::summarize_stats_factors(self$data, self$config)
      } else {
        if (factorDepth >= 1) {
          cfg <- self$config$clone(deep = TRUE)
          cfg$table$factorDepth <- factorDepth
          missing <- prolfqua::summarize_stats_factors(self$data, cfg)
        } else {
          missing <- prolfqua::summarize_stats_all(self$data, self$config)
        }
      }
      notNA <- missing |> dplyr::filter(nrNAs <= nrNA)
      sumN <- notNA |> dplyr::group_by_at(self$config$table$hierarchy_keys()) |>
        dplyr::summarise(n = dplyr::n())
      notNA <- sumN |> dplyr::filter(n == max(n))
      notNA <- notNA |> dplyr::select(self$config$table$hierarchy_keys())
      notNAdata <- dplyr::inner_join(notNA, self$data) |> dplyr::ungroup()
      return(LFQDataSE$new(data = notNAdata, config = self$config$clone(deep = TRUE), prefix = self$prefix))
    },

    #' @description
    #' complete cases only
    #' @return self
    complete_cases = function() {
      data <- prolfqua::complete_cases(self$data, self$config)
      self$.se <- private$long_to_se(data)
      self$.data_dirty <- TRUE
      invisible(self)
    },

    #' @description
    #' converts the data to wide format
    #' @param as.matrix return as matrix instead of data.frame
    #' @param value value column to use
    #' @return list with data, annotation, and configuration
    to_wide = function(as.matrix = FALSE, value = NULL) {
      # Use SE directly for performance
      if (is.null(value)) {
        mat <- SummarizedExperiment::assay(self$.se)
      } else {
        stopifnot(value %in% SummarizedExperiment::assayNames(self$.se))
        mat <- SummarizedExperiment::assay(self$.se, value)
      }

      annotation <- as.data.frame(SummarizedExperiment::colData(self$.se))
      annotation[[self$config$table$sampleName]] <- rownames(annotation)
      rowdata <- as.data.frame(SummarizedExperiment::rowData(self$.se))

      if (!as.matrix) {
        mat <- as.data.frame(mat)
        mat <- cbind(rowdata, mat)
      }

      result <- list(
        data = mat,
        annotation = annotation,
        rowdata = rowdata,
        config = self$config$clone(deep = TRUE)
      )
      return(result)
    },

    #' @description
    #' Annotation table
    #' @return data.frame
    factors = function() {
      prolfqua::table_factors(self$data, self$config)
    },

    #' @description
    #' Hierarchy table
    hierarchy = function() {
      hk <- self$config$table$hierarchy_keys_depth()
      hkdf <- self$data |> dplyr::select(dplyr::all_of(hk)) |> dplyr::distinct()
      return(hkdf)
    },

    #' @description
    #' name of response variable
    #' @return character
    response = function() {
      self$config$table$get_response()
    },

    #' @description
    #' rename response variable
    #' @param newname new name for response variable
    rename_response = function(newname = "Intensity") {
      if (newname %in% colnames(self$data)) {
        msg <- paste(newname, " already in data :", paste(colnames(self$data), collapse = " "), ".")
        message(msg)
      } else {
        old <- self$config$table$pop_response()
        self$config$table$set_response(newname)
        # Update SE assay name
        names(SummarizedExperiment::assays(self$.se))[1] <- newname
        self$.data_dirty <- TRUE
      }
    },

    #' @description
    #' number of elements at each level
    hierarchy_counts = function() {
      # Can compute directly from SE dimensions
      prolfqua::hierarchy_counts(self$data, self$config)
    },

    #' @description
    #' e.g. number of peptides per protein etc
    #' @return data.frame
    summarize_hierarchy = function() {
      prolfqua::summarize_hierarchy(self$data, self$config)
    },

    #' @description
    #' get Plotter
    #' @return LFQDataPlotter
    get_Plotter = function() {
      return(LFQDataPlotter$new(self, self$prefix))
    },

    #' @description
    #' get Summariser
    #' @return LFQDataSummariser
    get_Summariser = function() {
      return(LFQDataSummariser$new(self))
    },

    #' @description
    #' Get LFQDataStats
    #' @param stats default interaction
    #' @return LFQDataStats
    get_Stats = function(stats = c("everything", "interaction", "all")) {
      stats <- match.arg(stats)
      return(LFQDataStats$new(self, stats = stats))
    },

    #' @description
    #' get Transformer
    #' @return LFQDataTransformer
    get_Transformer = function() {
      return(LFQDataTransformer$new(self))
    },

    #' @description
    #' get Aggregator
    #' @return LFQDataAggregator
    get_Aggregator = function() {
      return(LFQDataAggregator$new(self))
    },

    #' @description
    #' get difference of self with other
    #' @param other a filtered LFQDataSE set
    #' @return LFQDataSE
    filter_difference = function(other) {
      diffdata <- prolfqua::filter_difference(self$data, other$data, self$config)
      res <- LFQDataSE$new(data = diffdata, config = self$config$clone(deep = TRUE), prefix = self$prefix)
      return(res)
    }
  ),

  active = list(
    #' @field data access data in long format (computed from SE on demand)
    data = function(value) {
      if (missing(value)) {
        if (self$.data_dirty || is.null(self$.data_cache)) {
          self$.data_cache <- private$se_to_long()
          self$.data_dirty <- FALSE
        }
        return(self$.data_cache)
      } else {
        # Setting data: update SE
        self$.se <- private$long_to_se(value)
        self$.data_dirty <- TRUE
      }
    }
  ),

  private = list(
    #' Convert SummarizedExperiment to long format
    se_to_long = function() {
      # Get row and column annotations
      row_data <- as.data.frame(SummarizedExperiment::rowData(self$.se))
      col_data <- as.data.frame(SummarizedExperiment::colData(self$.se))
      sample_name <- self$config$table$sampleName

      # Ensure sample name column exists in col_data
      if (!sample_name %in% colnames(col_data)) {
        col_data[[sample_name]] <- rownames(col_data)
      }

      # Get all assay names (value variables)
      assay_names <- SummarizedExperiment::assayNames(self$.se)
      if (is.null(assay_names) || length(assay_names) == 0) {
        assay_names <- "value"
      }

      # Convert first assay to long format (primary response)
      mat <- SummarizedExperiment::assay(self$.se, 1)
      mat_df <- as.data.frame(mat)
      hierarchy_keys <- self$config$table$hierarchy_keys()

      # Add hierarchy columns from rowData
      for (hk in hierarchy_keys) {
        if (hk %in% colnames(row_data)) {
          mat_df[[hk]] <- row_data[[hk]]
        }
      }

      # If hierarchy key is rownames
      if (length(hierarchy_keys) == 1 && !hierarchy_keys[1] %in% colnames(mat_df)) {
        mat_df[[hierarchy_keys[1]]] <- rownames(mat)
      }

      # Pivot primary response to long format
      response_name <- assay_names[1]
      long_df <- tidyr::pivot_longer(
        mat_df,
        cols = colnames(mat),
        names_to = sample_name,
        values_to = response_name
      )

      # Add other assays (value variables) if present
      if (length(assay_names) > 1) {
        for (i in 2:length(assay_names)) {
          assay_mat <- SummarizedExperiment::assay(self$.se, i)
          assay_df <- as.data.frame(assay_mat)

          # Add hierarchy columns
          for (hk in hierarchy_keys) {
            if (hk %in% colnames(row_data)) {
              assay_df[[hk]] <- row_data[[hk]]
            }
          }
          if (length(hierarchy_keys) == 1 && !hierarchy_keys[1] %in% colnames(assay_df)) {
            assay_df[[hierarchy_keys[1]]] <- rownames(assay_mat)
          }

          # Pivot to long and join
          assay_long <- tidyr::pivot_longer(
            assay_df,
            cols = colnames(assay_mat),
            names_to = sample_name,
            values_to = assay_names[i]
          )

          # Join with main long_df
          join_cols <- c(hierarchy_keys, sample_name)
          long_df <- dplyr::left_join(long_df, assay_long, by = join_cols)
        }
      }

      # Add factor columns from colData
      long_df <- dplyr::left_join(long_df, col_data, by = sample_name)

      return(long_df)
    },

    #' Convert long format to SummarizedExperiment
    long_to_se = function(data) {
      # Get all value variables from config
      value_vars <- self$config$table$value_vars()
      response_name <- self$config$table$get_response()

      # Create assay list - one per value variable
      assay_list <- list()
      sample_col <- self$config$table$sampleName
      hierarchy_keys <- self$config$table$hierarchy_keys()

      for (vv in value_vars) {
        if (vv %in% colnames(data)) {
          # Create config copy for this value variable
          cfg_copy <- self$config$clone(deep = TRUE)
          cfg_copy$table$set_response(vv)

          wide <- prolfqua::tidy_to_wide_config(data, cfg_copy, as.matrix = TRUE)
          assay_list[[vv]] <- wide$data
        }
      }

      # If no value vars found, use the response
      if (length(assay_list) == 0) {
        wide <- prolfqua::tidy_to_wide_config(data, self$config, as.matrix = TRUE)
        assay_list[[response_name]] <- wide$data
      } else {
        # Use the first conversion for row/col data
        cfg_copy <- self$config$clone(deep = TRUE)
        cfg_copy$table$set_response(value_vars[1])
        wide <- prolfqua::tidy_to_wide_config(data, cfg_copy, as.matrix = TRUE)
      }

      # Create column annotation
      ann <- data.frame(wide$annotation)
      if (sample_col %in% colnames(ann)) {
        rownames(ann) <- ann[[sample_col]]
      }

      # Create SE with all assays
      se <- SummarizedExperiment::SummarizedExperiment(
        assays = S4Vectors::SimpleList(assay_list),
        colData = ann,
        rowData = wide$rowdata,
        metadata = list(config = self$config, prefix = self$prefix)
      )

      return(se)
    }
  )
)


#' Convert LFQData to LFQDataSE
#'
#' @param lfqdata LFQData object
#' @return LFQDataSE object
#' @export
#' @family LFQData
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfq_se <- as_LFQDataSE(lfqdata)
#' class(lfq_se)
#'
as_LFQDataSE <- function(lfqdata) {
  LFQDataSE$new(data = lfqdata$data, config = lfqdata$config, prefix = lfqdata$prefix)
}


#' Convert LFQDataSE to LFQData
#'
#' @param lfqdatase LFQDataSE object
#' @return LFQData object
#' @export
#' @family LFQData
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' lfq_se <- LFQDataSE$new(data = istar$data, config = istar$config)
#' lfqdata <- as_LFQData(lfq_se)
#' class(lfqdata)
#'
as_LFQData <- function(lfqdatase) {
  LFQData$new(lfqdatase$data, lfqdatase$config, prefix = lfqdatase$prefix)
}


#' Create LFQData or LFQDataSE based on global option
#'
#' This function creates either an LFQData or LFQDataSE object depending on
#' the value of `getOption("prolfqua.default_class")`.
#'
#' @param data data.frame in long format
#' @param config AnalysisConfiguration object
#' @param prefix output prefix (default "ms_")
#' @param setup if TRUE, calls setup_analysis on data first (default FALSE)
#' @return LFQData or LFQDataSE depending on getOption("prolfqua.default_class")
#' @export
#' @family LFQData
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' lfqdata <- make_LFQData(istar$data, istar$config)
#' class(lfqdata)  # LFQData by default
#'
#' # Switch to LFQDataSE
#' options(prolfqua.default_class = "LFQDataSE")
#' lfq_se <- make_LFQData(istar$data, istar$config)
#' class(lfq_se)  # LFQDataSE
#' options(prolfqua.default_class = NULL)  # reset
#'
make_LFQData <- function(data, config, prefix = "ms_", setup = FALSE) {
  class_type <- getOption("prolfqua.default_class", default = "LFQData")
  if (class_type == "LFQDataSE") {
    LFQDataSE$new(data = data, config = config, prefix = prefix, setup = setup)
  } else {
    LFQData$new(data = data, config = config, prefix = prefix, setup = setup)
  }
}


#' Ensure object is LFQData or LFQDataSE
#'
#' Converts a list with data/config to LFQData or LFQDataSE, or returns
#' an existing LFQData/LFQDataSE object unchanged.
#'
#' @param x list with data and config elements, or LFQData/LFQDataSE object
#' @return LFQData or LFQDataSE (respects global option for new objects)
#' @export
#' @family LFQData
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' lfqdata <- as_lfq(istar)
#' class(lfqdata)
#'
#' # Already an LFQData - returns unchanged
#' lfqdata2 <- as_lfq(lfqdata)
#' identical(lfqdata, lfqdata2)  # TRUE
#'
as_lfq <- function(x) {
  if (inherits(x, c("LFQData", "LFQDataSE"))) return(x)
  make_LFQData(x$data, x$config)
}
