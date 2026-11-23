# LFQDataTransformer ----

#' Decorate LFQData with Methods for transforming Intensities
#'
#' @return An R6 class generator.
#' @export
#'
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config(n_proteins = 50)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#'
#' lfqcopy <- lfqdata$get_copy()
#' lfqTrans <- lfqcopy$get_Transformer()
#'
#' x <- lfqTrans$intensity_array(log2)
#'
#' x$lfq$is_transformed()
#' x <- x$intensity_matrix()
#' plotter <- x$lfq$get_Plotter()
#' plotter$intensity_distribution_density()
#'
#' # transform by asinh root and scale
#'
#' lfqcopy <- lfqdata$get_copy()
#' lfqTrans <- lfqcopy$get_Transformer()
#' x <- lfqTrans$intensity_array(asinh)
#' mads1 <- mean(x$get_scales()$mads)
#' x <- lfqTrans$intensity_matrix(force = TRUE)
#' mads2 <- mean(x$get_scales()$mads)
#' stopifnot(abs(mads1 - mads2) < 1e-8)
#'
#'
#' stopifnot(abs(mean(x$get_scales()$medians)) < 1e-8)
#' lfqcopy <- lfqdata$get_copy()
#' lfqTrans <- lfqcopy$get_Transformer()
#' lfqTrans$log2()
#' before <- lfqTrans$get_scales()
#' lfqTrans$robscale()
#' after <- lfqTrans$get_scales()
#' stopifnot(abs(mean(before$medians) - mean(after$medians)) < 1e-8)
#' stopifnot(abs(mean(before$mads) - mean(after$mads)) < 1e-8)
#'
#' # normalize data using vsn
#' lfqcopy <- lfqdata$get_copy()
#' lfqTrans <- lfqcopy$get_Transformer()
#' lfqTransCheck <- lfqcopy$get_Transformer()
#'
#' lfqTransCheck$log2()
#' lfqTransCheck$get_scales()
#' lfqTransCheck$lfq$get_Plotter()$intensity_distribution_density()
#'
#' if (requireNamespace("vsn", quietly = TRUE)) {
#'  res <- lfqTrans$intensity_matrix( .func = vsn::justvsn)
#'  res$lfq$get_Plotter()$intensity_distribution_density()
#'  res$get_scales()
#' }
#' if (requireNamespace("preprocessCore", quietly = TRUE)) {
#' quant <- function(y){
#'  ynorm <- preprocessCore::normalize.quantiles(y)
#'  rownames(ynorm) <- rownames(y)
#'  colnames(ynorm) <- colnames(y)
#'  return(ynorm)
#' }
#'  res <- lfqTrans$intensity_matrix( .func = quant)
#'  res$lfq$get_Plotter()$intensity_distribution_density()
#' }
#'
#'
#' #subset scaling
#'
#' istar2 <- sim_lfq_data_peptide_config()
#' lfqdata2 <- LFQData$new(istar2$data, istar2$config)
#' lfqdata2 <- lfqdata2$get_Transformer()$intensity_array(log2)$lfq
#' head(lfqdata2$hierarchy())
#' internal <- lfqdata2$get_subset(head(lfqdata2$hierarchy()))
#' internal$hierarchy()
#' tr <- lfqdata2$get_Transformer()
#' tr$center_to_reference(internal)
#' pl <- tr$lfq$get_Plotter()
#' pl$intensity_distribution_density()
#' lfqdata2$get_Plotter()$intensity_distribution_density()
#' robscale <- lfqdata2$get_Transformer()$robscale_subset(internal)$lfq
#' robscale$get_Plotter()$intensity_distribution_density()
#'
#' # Apply a one-unit log2 offset to the first sample
#' tr <- lfqdata$get_Transformer()$log2()
#' samples <- unique(tr$lfq$data_long()[[tr$lfq$sample_name()]])
#' offsets <- setNames(numeric(length(samples)), samples)
#' offsets[1] <- 1
#' tr$subtract_sample_offsets(offsets)
#'
#' @srrstats {G2.16} Log transforms turn zeros into NA and warn on negatives; Inf, -Inf and NaN responses raise a
#'   classed error before every fit and aggregation. Limit: no option to ignore or remove non-finite values, by design.
#' @section Centring and contrasts:
#' Centring is an explicit transformation step, never applied by a model. `robscale()` and `robscale_subset()`
#' centre each sample at its median (and scale by its MAD); `preserve_mean = TRUE` adds the overall mean back, so
#' abundances stay on the original level, `FALSE` centres them at zero. `center_to_reference()` and
#' `subtract_sample_offsets()` shift each sample by one offset. A per-sample shift moves every feature of that
#' sample by the same amount: a contrast between groups changes by the difference of the mean offsets of the groups,
#' and a shift common to all samples (such as `preserve_mean`) changes the average abundance (`avgAbd`) but not the
#' contrast estimates, standard errors or p-values.
#' @srrstats {RE2.3} Centring is available through explicit transformer parameters (`preserve_mean`, reference
#'   subsets, sample offsets); its effect on contrasts is documented in the section "Centring and contrasts" and
#'   tested in test-centring.R.
LFQDataTransformer <- R6::R6Class(
  "LFQDataTransformer",
  public = list(
    #' @field lfq LFQData object
    lfq = NULL,
    #' @description
    #' initialize
    #' @param lfqdata
    #' LFQData object to transform
    initialize = function(lfqdata) {
      self$lfq <- lfqdata$clone(deep = TRUE)
    },
    #' @description
    #' log2 transform data
    #' @param force if `FALSE`, data already log2 transformed will not be
    #'   transformed a second time. `TRUE` forces log transformation
    #' @return LFQDataTransformer
    log2 = function(force = FALSE) {
      self$intensity_array(log2, force = force)
    },
    #' @description
    #' get mean and variance and standard deviation in each sample
    #' @return list with means and mads
    get_scales = function() {
      .lfqdata_robscales(self$lfq)
    },
    #' @description
    #' robust scale data
    #' @param colname new name of transformed column
    #' @param preserve_mean should original mean value be preserved TRUE, if FALSE then center at zero; changes the
    #'   average abundance, not the contrasts (see section "Centring and contrasts")
    #' @return LFQDataTransformer (self)
    robscale = function(preserve_mean = TRUE, colname = "transformedIntensity") {
      res <- self$robscale_subset(self$lfq, preserve_mean = preserve_mean, colname = colname)
      invisible(res)
    },
    #' @description
    #' log2 transform and robust scale data based on subset
    #' @param lfqsubset LFQData subset to use for normalization
    #' @param preserve_mean should original mean value be preserved TRUE, if FALSE then center at zero
    #' @param colname how to name the transformed intensities, default transformedIntensity
    #' @return LFQDataTransformer (self)
    #'
    robscale_subset = function(lfqsubset, preserve_mean = TRUE, colname = "transformed_abundance") {
      if (self$lfq$is_transformed() != lfqsubset$is_transformed()) {
        abort_bad_argument("lfqsubset", "have the same transformation state (is_transformed) as self")
      }
      cfg <- self$lfq$get_config()$clone(deep = TRUE)
      scaled <- .scale_with_subset(self$lfq, lfqsubset, preserve_mean = preserve_mean, colname = colname)
      cfg$set_response(scaled$colname)
      self$lfq <- LFQData$new(scaled$data, cfg, provenance = self$lfq$get_provenance_log())
      invisible(self)
    },
    #' @description
    #' center data to a reference subset
    #' @param lfqsubset LFQData subset to use as reference
    #' @return LFQDataTransformer (self)
    #'
    center_to_reference = function(lfqsubset) {
      if (self$lfq$is_transformed() != lfqsubset$is_transformed()) {
        abort_bad_argument("lfqsubset", "have the same transformation state (is_transformed) as self")
      }
      if (!self$lfq$is_transformed()) {
        warning("The data should be log2-transformed.")
      }
      self$lfq <- center_to_reference_cfg(self$lfq, lfqsubset)
      invisible(self)
    },

    #' @description
    #' Subtract caller-supplied offsets from the active log2 response. The caller supplies
    #' log2 data and chooses the global level through the offsets. Missing intensities
    #' remain missing. The result becomes the active response in a new column named
    #' \code{<response>_offset_adjusted}, with a numeric suffix if the name already exists.
    #' @param offsets Named numeric vector with one finite log2 offset per sample.
    #'   Names must match the configured sample IDs exactly and may be in any order.
    #' @return LFQDataTransformer (self, invisibly)
    subtract_sample_offsets = function(offsets) {
      data <- self$lfq$data_long()
      sample_ids <- as.character(data[[self$lfq$sample_name()]])
      matched_offsets <- .match_sample_offsets(offsets, sample_ids)
      response <- self$lfq$response()
      colname <- tail(make.unique(c(names(data), paste0(response, "_offset_adjusted"))), 1)
      data[[colname]] <- data[[response]] - matched_offsets
      cfg <- self$lfq$get_config()$clone(deep = TRUE)
      cfg$set_response(colname)
      self$lfq <- LFQData$new(data, cfg, prefix = self$lfq$prefix, provenance = self$lfq$get_provenance_log())
      invisible(self)
    },

    #' @description
    #' Transforms intensities
    #' @param .func transformation function working with arrays e.g. log2, log10, asinh etc.
    #' @param force transformation on already transformed data.
    #' @param transformed whether the result counts as transformed; set
    #'   \code{FALSE} for a back-transformation such as \code{exp}
    #'
    #' @return LFQDataTransformer (self)
    #'
    intensity_array = function(.func = log2, force = FALSE, transformed = TRUE) {
      .funcname <- deparse(match.call()$.func)
      private$.set_transformed(
        force,
        .transform_work_intensity(self$lfq$data_long(), self$lfq$response(), .func, .funcname),
        transformed = transformed
      )
    },
    #' @description
    #' pass a function which works with matrices, e.g., vsn::justvsn
    #' @param .func any function taking a matrix and returning a matrix
    #'   (columns sample, rows feature, e.g. `base::scale()`); default: robust scaling, which centres each
    #'   sample on its median and scales it by its MAD, as `$robscale()`
    #' @param force transformation on data already transformed
    #' @return LFQDataTransformer (self)
    #'
    intensity_matrix = function(.func = .robust_scale, force = FALSE) {
      .funcname <- deparse(match.call()$.func)
      private$.set_transformed(force, .apply_to_response_matrix(self$lfq, .func, .funcname))
    }
  ),
  private = list(
    # `res` (list with data and colname) is a promise: the transformation only runs if it is applied
    .set_transformed = function(force, res, transformed = TRUE) {
      if (self$lfq$is_transformed() && !force) {
        warning("The data are already transformed; set `force = TRUE` to log2-transform them again.")
        return(invisible(self))
      }
      cfg <- self$lfq$get_config()$clone(deep = TRUE)
      cfg$set_response(res$colname)
      cfg$is_response_transformed <- transformed
      response <- self$lfq$response()
      result <- LFQData$new(res$data, cfg, provenance = self$lfq$get_provenance_log())
      zeros <- res$zeros %||% rep(FALSE, nrow(res$data))
      log <- result$get_provenance_log()
      log$record(
        "log transform",
        "values set to NA",
        .key_ids(res$data[zeros, ], self$lfq$hierarchy_keys()),
        "zero intensity (not detected)",
        n = sum(zeros),
        response = response
      )
      log$report()
      self$lfq <- result
      invisible(self)
    }
  )
)

# Intensity transformation helpers ----

#' Validate log2 offsets and align them to the sample ID of each measurement
#' @param offsets Named numeric vector of finite log2 offsets, one per sample.
#' @param sample_ids Character vector of sample IDs in measurement order.
#' @return Unnamed numeric vector of offsets in measurement order.
#' @noRd
.match_sample_offsets <- function(offsets, sample_ids) {
  if (!is.numeric(offsets) || !is.null(dim(offsets)) || any(!is.finite(offsets))) {
    abort_bad_argument("offsets", "be a numeric vector of finite log2 offsets")
  }
  offset_ids <- names(offsets)
  valid_ids <- all(
    !is.null(offset_ids),
    !anyNA(offset_ids),
    all(nzchar(offset_ids)),
    anyDuplicated(offset_ids) == 0L
  )
  if (!valid_ids) {
    abort_bad_argument("offsets", "have unique, nonempty sample IDs as names")
  }
  if (!setequal(offset_ids, sample_ids)) {
    abort_bad_argument("offsets", "name exactly the samples in the data")
  }
  unname(offsets[sample_ids])
}

#' Transform intensity column by applying a function
#' @param pdata data.frame
#' @param response character, name of the response column to transform
#' @param .func function to transform intensities e.g. log2
#' @param .funcname name of function (used for creating new column name)
#' @return list with `data` (data.frame), `colname` (new column name) and `zeros` (logical, rows whose zero
#'   intensity became `NA`)
#' @details For `log`, `log2` and `log10`, zero intensities mean "not detected" and become `NA`; LFQData records
#'   them in its provenance. Negative intensities give `NaN` with a warning, which fitting and aggregation then
#'   reject.
#' @noRd
.transform_work_intensity <- function(pdata, response, .func, .funcname = NULL) {
  newcol <- paste(.funcname %||% deparse(match.call()$.func), response, sep = "_")
  vals <- pdata[[response]]
  zeros <- rep(FALSE, nrow(pdata))
  if (identical(.func, log2) || identical(.func, log) || identical(.func, log10)) {
    n_neg <- sum(vals < 0, na.rm = TRUE)
    if (n_neg > 0) {
      warning(
        "Log transform: ",
        n_neg,
        " negative values in `",
        response,
        "` will produce NaN. Consider filtering or using asinh."
      )
    }
    # 0 means not detected: a missing observation, not -Inf.
    zeros <- !is.na(vals) & vals == 0
    pdata[[response]][zeros] <- NA
  }

  pdata <- pdata |>
    dplyr::mutate(
      !!newcol := .func(.data[[response]])
    )
  inform_prolfqua("Column added : ", newcol)
  list(data = pdata, colname = newcol, zeros = zeros)
}

#' Takes matrix of responses and converts into tibble
#'
#' @param pdata matrix with rownames encoding hierarchy keys
#' @param value name of column to store values in
#' @param config AnalysisConfiguration (needed for column name mapping)
#' @param data optional data.frame to join back to
#' @param sep separator used to unite the hierarchy keys
#' @return a long-format tibble with the hierarchy key, sample name and `value` columns; if `data`
#'   is given, `data` joined with these values
#' @noRd
.response_matrix_as_tibble <- function(pdata, value, config, data = NULL, sep = "~lfq~") {
  pdata <- dplyr::bind_cols(
    tibble::tibble("row.names" := rownames(pdata)),
    tibble::as_tibble(pdata)
  )
  pdata <- tidyr::pivot_longer(
    pdata,
    cols = -1,
    names_to = config$sample_name,
    values_to = value
  )
  pdata <- tidyr::separate_wider_delim(
    pdata,
    "row.names",
    delim = sep,
    names = config$hierarchy_keys()
  )
  if (!is.null(data)) {
    pdata <- dplyr::inner_join(
      data,
      pdata,
      by = unique(c(config$hierarchy_keys(), config$sample_name)),
      relationship = "one-to-one"
    )
  }
  pdata
}

#' Compute medians and MADs of a matrix along one margin
#' @param data numeric matrix, samples in columns
#' @param dim margin to summarise: 2 for columns (samples), 1 for rows
#' @keywords internal
#'
.get_robscales <- function(data, dim = 2) {
  medians <- apply(data, dim, median, na.rm = TRUE)
  data <- sweep(data, dim, medians, "-")
  mads <- apply(data, dim, mad, na.rm = TRUE)
  list(medians = medians, mads = mads)
}

#' Compute the median and MAD of each sample
#' @param lfqdata [LFQData] object with the data and its configuration
#' @return a list with the numeric vectors `medians` and `mads`, with one element per sample
#' @noRd
.lfqdata_robscales <- function(lfqdata) {
  data <- .tidy_to_wide_config(lfqdata, as.matrix = TRUE)$data
  scales <- .get_robscales(data)
  scales
}


#' Robustly scale a matrix by medians and MADs
#' @param data numeric matrix
#' @param dim margin to scale: 2 scales each column (sample), 1 each row
#' @param preserve_mean if TRUE, add the mean of the medians back
#' @return a matrix like `data`, centred by the median of each column (`dim = 2`) or row (`dim = 1`) and divided
#'   by its MAD relative to the mean MAD; only centred, with a warning, if a MAD is zero; with
#'   `preserve_mean = TRUE` the mean of the medians is added back
#' @srrstats {G5.8c} Tests cover zero MAD, constant heatmaps, missing groups, and proteins with all-NA or constant
#'   responses in model fits (test-missing-nonfinite.R).
#' @noRd
.robust_scale <- function(data, dim = 2, preserve_mean = FALSE) {
  .apply_robscales(
    data,
    .get_robscales(data, dim = dim),
    preserve_mean,
    dim,
    "Skipping the scaling step in `robust_scale()`: one or more MAD values are zero."
  )
}

#' center by the medians and scale by the relative mads (see .get_robscales); warns with `zero_mad_warning`
#'
#' and skips the scaling when a mad is zero
#' @param data numeric measurement matrix
#' @param scales per-margin medians and median absolute deviations
#' @param preserve_mean restore the overall mean after scaling
#' @param dim matrix margin to scale
#' @param zero_mad_warning warning emitted when scaling cannot divide by the MAD
#' @noRd
.apply_robscales <- function(data, scales, preserve_mean, dim, zero_mad_warning) {
  data <- sweep(data, dim, scales$medians, "-")
  if (any(.is_zero(scales$mads))) {
    warning(zero_mad_warning)
  } else {
    data <- sweep(data, dim, scales$mads / mean(scales$mads), "/")
  }
  data + if (preserve_mean) mean(scales$medians) else 0
}


#' Apply function requiring a matrix to tidy table
#'
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param .func function taking and returning a matrix
#' @param .funcname name of function (used for creating new column)
#' @return list with `data` (data.frame) and `colname` (new column name)
#' @noRd
.apply_to_response_matrix <- function(lfqdata, .func, .funcname = NULL) {
  colname <- make.names(paste(lfqdata$response(), .funcname %||% deparse(match.call()$.func), sep = "_"))
  mat <- .tidy_to_wide_config(lfqdata, as.matrix = TRUE)$data
  mat <- .func(mat)
  data <- .response_matrix_as_tibble(
    mat,
    colname,
    lfqdata$get_config(),
    lfqdata$data_long()
  )
  list(data = data, colname = colname)
}

#' Scale data using a subset of the data
#'
#' this should reduce the overall variance.
#'
#' @param lfqdata LFQData object with full dataset
#' @param lfqsubset LFQData object with subset for computing scales
#' @param preserve_mean default FALSE - sets mean to zero
#' @param colname name of the scaled intensity column; default NULL uses \code{<response>_subset_scaled}
#' @return list with data, scales, and colname
#' @noRd
.scale_with_subset <- function(lfqdata, lfqsubset, preserve_mean = FALSE, colname = NULL) {
  scaled_name <- make.names(paste(lfqdata$response(), "subset_scaled", sep = "_"))
  scales <- .get_robscales(.tidy_to_wide_config(lfqsubset, as.matrix = TRUE)$data)
  mat <- .apply_robscales(
    .tidy_to_wide_config(lfqdata, as.matrix = TRUE)$data,
    scales,
    preserve_mean,
    2,
    "Skipping the scaling step in `scale_with_subset()`: one or more MAD values are zero."
  )
  data <- .response_matrix_as_tibble(mat, scaled_name, lfqdata$get_config(), lfqdata$data_long())
  if (!is.null(colname)) {
    data <- dplyr::rename(data, !!colname := dplyr::all_of(scaled_name))
  }
  list(data = data, scales = scales, colname = colname %||% scaled_name)
}


#' Function to normalize protein abundances by subtracting sample means of reference proteins
#' @param df long-format protein measurements to center
#' @param df_reference reference protein measurements
#' @param sample_name column identifying samples
#' @param abundance_column column containing intensities
#' @noRd
center_to_reference <- function(
  df,
  df_reference,
  sample_name,
  abundance_column = "normalized_abundance"
) {
  # Step 1: Calculate sample means for reference proteins
  sample_means <- df_reference |>
    dplyr::group_by(.data[[sample_name]]) |>
    dplyr::summarise(
      reference_mean = mean(.data[[abundance_column]], na.rm = TRUE),
      reference_median = median(.data[[abundance_column]], na.rm = TRUE),
      .groups = "drop"
    )
  # Step 2: Join back to original data and subtract sample means
  normalized_df <-
    dplyr::left_join(df, sample_means, by = sample_name, relationship = "many-to-one") |>
    dplyr::mutate(
      # Create normalized abundance column
      centered_abundance_by_mean = .data[[abundance_column]] - reference_mean,
      centered_abundance_by_median = .data[[abundance_column]] - reference_median
    ) |>
    dplyr::select(-reference_mean, -reference_median) # Remove the temporary column
  normalized_df
}

#' Center data to a reference subset
#'
#' takes the median of the lfqdareference per sample and subtracts it from a copy of lfqdata
#' @param lfqdata LFQData object containing the data to center
#' @param lfqdareference LFQData object containing the reference subset
#' @return The computed result.
#' @export
#' @examples
#' # example code
#'
#' bb <- sim_lfq_data_peptide_config(n_proteins = 100)
#' x <- LFQData$new(bb$data, bb$config)
#' xc <- x$get_copy()
#' xc$set_data(xc$data_long() |> dplyr::filter(protein_Id == "0EfVhX~3967"))
#' xxd <- center_to_reference_cfg(x, xc)
#' xxd$response()
#' x$response()
#'
center_to_reference_cfg <- function(lfqdata, lfqdareference) {
  resdata <- lfqdata$get_copy()
  cfg <- resdata$get_config()
  data <- center_to_reference(lfqdata$data_long(), lfqdareference$data_long(), cfg$sample_name, cfg$get_response())
  resdata$set_data(data)
  cfg$set_response("centered_abundance_by_median")
  invisible(resdata)
}
