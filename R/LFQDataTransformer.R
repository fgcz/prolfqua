# LFQDataTransformer ----

#' Decorate LFQData with Methods for transforming Intensities
#'
#' @return An R6 class generator.
#' @export
#'
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config(Nprot = 50)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#'
#' lfqcopy <- lfqdata$get_copy()
#' lfqTrans <- lfqcopy$get_Transformer()
#'
#' x <- lfqTrans$intensity_array(log2)
#'
#' x$lfq$is_transformed()
#' x <- x$intensity_matrix(robust_scale)
#' plotter <- x$lfq$get_Plotter()
#' plotter$intensity_distribution_density()
#'
#' # transform by asinh root and scale
#'
#' lfqcopy <- lfqdata$get_copy()
#' lfqTrans <- lfqcopy$get_Transformer()
#' x <- lfqTrans$intensity_array(asinh)
#' mads1 <- mean(x$get_scales()$mads)
#' x <- lfqTrans$intensity_matrix(robust_scale, force = TRUE)
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
#' if(require("vsn")){
#'  res <- lfqTrans$intensity_matrix( .func = vsn::justvsn)
#'  res$lfq$get_Plotter()$intensity_distribution_density()
#'  res$get_scales()
#' }
#' if(require("preprocessCore")){
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
      get_robscales(self$lfq)
    },
    #' @description
    #' robust scale data
    #' @param colname new name of transformed column
    #' @param preserve_mean should original mean value be preserved TRUE, if FALSE then center at zero
    #' @return LFQDataTransformer (self)
    robscale = function(preserve_mean = TRUE, colname = "transformedIntensity") {
      res <- self$robscale_subset(self$lfq, preserve_mean = preserve_mean, colname = colname)
      invisible(res)
    },
    #' @description
    #' log2 transform and robust scale data based on subset
    #' @param lfqsubset LFQData subset to use for normalization
    #' @param preserve_mean should original mean value be preserved TRUE, if FALSE then center at zero
    #' @param colname - how to name the transformed intensities, default transformedIntensity
    #' @return LFQDataTransformer (self)
    #'
    robscale_subset = function(lfqsubset, preserve_mean = TRUE, colname = "transformed_abundance") {
      message("data is : ", self$lfq$is_transformed())
      if (self$lfq$is_transformed() != lfqsubset$is_transformed()) {
        warning("the subset must have the same config as self")
      }
      cfg <- self$lfq$get_config()$clone(deep = TRUE)
      scaled <- prolfqua::scale_with_subset(self$lfq, lfqsubset, preserve_mean = preserve_mean, colname = colname)
      cfg$set_response(scaled$colname)
      self$lfq <- LFQData$new(scaled$data, cfg)
      invisible(self)
    },
    #' @description
    #' center data to a reference subset
    #' @param lfqsubset LFQData subset to use as reference
    #' @return LFQDataTransformer (self)
    #'
    center_to_reference = function(lfqsubset) {
      message("data is transformed: ", self$lfq$is_transformed())
      if (self$lfq$is_transformed() != lfqsubset$is_transformed()) {
        abort_bad_argument("lfqsubset", "have the same transformation state (is_transformed) as self")
      }
      if (!self$lfq$is_transformed()) {
        warning("data should be log2 transformed")
      }
      self$lfq <- center_to_reference_cfg(self$lfq, lfqsubset)
      invisible(self)
    },

    #' @description
    #' Transforms intensities
    #' @param .func transformation function working with arrays e.g. log2, log10, asinh etc.
    #' @param force transformation on already transformed data.
    #'
    #' @return LFQDataTransformer (self)
    #'
    intensity_array = function(.func = log2, force = FALSE) {
      .funcname <- deparse(match.call()$.func)
      private$.set_transformed(
        force,
        prolfqua::transform_work_intensity(self$lfq$data_long(), self$lfq$response(), .func, .funcname)
      )
    },
    #' @description
    #' pass a function which works with matrices, e.g., vsn::justvsn
    #' @param .func any function taking a matrix and returning a matrix
    #'   (columns sample, rows feature, e.g. `base::scale()`), default
    #'   `robust_scale`
    #' @param force transformation on data already transformed
    #' @return LFQDataTransformer (self)
    #'
    intensity_matrix = function(.func = robust_scale, force = FALSE) {
      .funcname <- deparse(match.call()$.func)
      private$.set_transformed(force, prolfqua::apply_to_response_matrix(self$lfq, .func, .funcname))
    }
  ),
  private = list(
    # `res` (list with data and colname) is a promise: the transformation only runs if it is applied
    .set_transformed = function(force, res) {
      if (self$lfq$is_transformed() && !force) {
        warning("data already transformed. If you still want to log2 tranform, set force = TRUE")
        return(invisible(self))
      }
      cfg <- self$lfq$get_config()$clone(deep = TRUE)
      cfg$set_response(res$colname)
      cfg$is_response_transformed <- TRUE
      self$lfq <- LFQData$new(res$data, cfg)
      invisible(self)
    }
  )
)

# Intensity transformation helpers ----

#' Transform intensity column by applying a function
#' @param pdata data.frame
#' @param response character, name of the response column to transform
#' @param .func function to transform intensities e.g. log2
#' @param .funcname name of function (used for creating new column name)
#' @return list with `data` (data.frame) and `colname` (new column name)
#' @export
#' @keywords internal
#' @examples
#'
#' dd <- prolfqua_data('data_spectronautDIA250_A')
#' config <- dd$config_f()
#' analysis <- dd$analysis(dd$data, config)
#' res <- transform_work_intensity(
#'   analysis, config$get_response(), .func = log2
#' )
#' stopifnot("log2_FG.Quantity" %in% colnames(res$data))
#' config <- dd$config_f()
#' res <- transform_work_intensity(
#'   analysis, config$get_response(), .func = asinh
#' )
#' stopifnot("asinh_FG.Quantity" %in% colnames(res$data))
#'
transform_work_intensity <- function(pdata, response, .func, .funcname = NULL) {
  newcol <- paste(.funcname %||% deparse(match.call()$.func), response, sep = "_")
  vals <- pdata[[response]]
  if (identical(.func, log2) || identical(.func, log) || identical(.func, log10)) {
    n_zero <- sum(vals == 0, na.rm = TRUE)
    n_neg <- sum(vals < 0, na.rm = TRUE)
    if (n_neg > 0) {
      warning(
        "log transform: ",
        n_neg,
        " negative values in '",
        response,
        "' will produce NaN. Consider filtering or using asinh."
      )
    }
    if (n_zero > 0) {
      warning(
        "log transform: ",
        n_zero,
        " zeros in '",
        response,
        "' will produce -Inf. Consider replacing zeros with NA first."
      )
    }
  }

  pdata <- pdata |>
    dplyr::mutate(
      !!sym(newcol) := .func(!!sym(response))
    )
  message("Column added : ", newcol)
  return(list(data = pdata, colname = newcol))
}

#' Takes matrix of responses and converts into tibble
#'
#' @param pdata matrix with rownames encoding hierarchy keys
#' @param value name of column to store values in
#' @param config AnalysisConfiguration (needed for column name mapping)
#' @param data optional data.frame to join back to
#' @param sep separator used to unite the hierarchy keys
#' @export
#'
#' @keywords internal
#' @examples
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' lfqdata <- prolfqua::LFQData$new(dd$data, dd$config)
#' res <- tidy_to_wide_config(lfqdata, as.matrix = TRUE)
#'
#' res_scaled <- scale(res$data)
#' xx <- response_matrix_as_tibble(
#'   res_scaled, "srm_intensityScaled", lfqdata$get_config()
#' )
#' xx <- response_matrix_as_tibble(
#'   res_scaled, "srm_intensityScaled",
#'   lfqdata$get_config(), lfqdata$data_long()
#' )
#'
response_matrix_as_tibble <- function(pdata, value, config, data = NULL, sep = "~lfq~") {
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
  pdata <- tidyr::separate(
    pdata,
    "row.names",
    unique(c(config$hierarchy_keys(), config$isotope_label)),
    sep = sep
  )
  if (!is.null(data)) {
    pdata <- dplyr::inner_join(data, pdata)
  }
  return(pdata)
}

#' compute median and mad on matrix
#' @keywords internal
#'
.get_robscales <- function(data, dim = 2) {
  medians <- apply(data, dim, median, na.rm = TRUE)
  data <- sweep(data, dim, medians, "-")
  mads <- apply(data, dim, mad, na.rm = TRUE)
  return(list(medians = medians, mads = mads))
}

#' compute median and standard deviation for each sample
#' @export
#' @keywords internal
#' @family preprocessing
#' @examples
#'
#'
#' bb <- prolfqua::sim_lfq_data_peptide_config()
#' lfqdata <- prolfqua::LFQData$new(bb$data, bb$config)
#' lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#' s1 <- get_robscales(lfqdata)
#'
#' res <- scale_with_subset(lfqdata, lfqdata)
#' lfqres <- prolfqua::LFQData$new(res$data, lfqdata$get_config()$clone(deep = TRUE))
#' s2 <- get_robscales(lfqres)
#' abs(mean(s1$mads) - mean(s2$mads)) < 0.1
#'
#'
get_robscales <- function(lfqdata) {
  data <- tidy_to_wide_config(lfqdata, as.matrix = TRUE)$data
  scales <- .get_robscales(data)
  return(scales)
}


#' robust scale wrapper
#' @keywords internal
#' @family preprocessing
#' @export
#' @examples
#' mat <- matrix(c(1, 2, 3, 4, 10, 12), nrow = 3)
#' scaled <- robust_scale(mat)
#' dim(scaled)
robust_scale <- function(data, dim = 2, preserve_mean = FALSE) {
  .apply_robscales(
    data,
    .get_robscales(data, dim = dim),
    preserve_mean,
    dim,
    "SKIPPING scaling step in robust_scale: one or more MAD values are zero."
  )
}

# center by the medians and scale by the relative mads (see .get_robscales); warns with `zero_mad_warning`
# and skips the scaling when a mad is zero
.apply_robscales <- function(data, scales, preserve_mean, dim, zero_mad_warning) {
  data <- sweep(data, dim, scales$medians, "-")
  if (!any(scales$mads == 0)) {
    data <- sweep(data, dim, scales$mads / mean(scales$mads), "/")
  } else {
    warning(zero_mad_warning)
  }
  data + if (preserve_mean) mean(scales$medians) else 0
}


#' Apply function requiring a matrix to tidy table
#'
#' @param lfqdata LFQData object
#' @param .func function taking and returning a matrix
#' @param .funcname name of function (used for creating new column)
#' @return list with `data` (data.frame) and `colname` (new column name)
#' @export
#' @keywords internal
#' @family preprocessing
#' @examples
#'
#' bb <- sim_lfq_data_peptide_config(Nprot = 100)
#' lfqdata <- LFQData$new(bb$data, bb$config)
#' res <- apply_to_response_matrix(lfqdata, .func = base::scale)
#' stopifnot("abundance_base..scale" %in% colnames(res$data))
#'
#' res <- apply_to_response_matrix(lfqdata, .func = robust_scale)
#'
#' if (require("vsn")) {
#'   res <- apply_to_response_matrix(lfqdata, .func = vsn::justvsn)
#' }
#'
apply_to_response_matrix <- function(lfqdata, .func, .funcname = NULL) {
  colname <- make.names(paste(lfqdata$response(), .funcname %||% deparse(match.call()$.func), sep = "_"))
  mat <- tidy_to_wide_config(lfqdata, as.matrix = TRUE)$data
  mat <- .func(mat)
  data <- response_matrix_as_tibble(
    mat,
    colname,
    lfqdata$get_config(),
    lfqdata$data_long()
  )
  return(list(data = data, colname = colname))
}

#' Scale data using a subset of the data
#'
#' this should reduce the overall variance.
#'
#' @export
#' @keywords internal
#' @param lfqdata LFQData object with full dataset
#' @param lfqsubset LFQData object with subset for computing scales
#' @param preserve_mean default FALSE - sets mean to zero
#' @param colname name of the scaled intensity column; default NULL uses \code{<response>_subset_scaled}
#' @family preprocessing
#' @return list with data, scales, and colname
#' @examples
#'
#' bb <- sim_lfq_data_peptide_config(Nprot = 100)
#' lfqdata <- LFQData$new(bb$data, bb$config)
#' lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#' s1 <- get_robscales(lfqdata)
#' res <- scale_with_subset(lfqdata, lfqdata)
#' cfg <- lfqdata$get_config()$clone(deep = TRUE)
#' cfg$set_response(res$colname)
#' lfqres <- LFQData$new(res$data, cfg)
#' s2 <- get_robscales(lfqres)
#' stopifnot(abs(mean(s1$mads) - mean(s2$mads)) < 1e-6)
scale_with_subset <- function(lfqdata, lfqsubset, preserve_mean = FALSE, colname = NULL) {
  scaled_name <- make.names(paste(lfqdata$response(), "subset_scaled", sep = "_"))
  scales <- .get_robscales(tidy_to_wide_config(lfqsubset, as.matrix = TRUE)$data)
  mat <- .apply_robscales(
    tidy_to_wide_config(lfqdata, as.matrix = TRUE)$data,
    scales,
    preserve_mean,
    2,
    "SKIPPING scaling step in scale_with_subset function."
  )
  data <- response_matrix_as_tibble(mat, scaled_name, lfqdata$get_config(), lfqdata$data_long())
  if (!is.null(colname)) {
    data <- dplyr::rename(data, !!colname := !!scaled_name)
  }
  return(list(data = data, scales = scales, colname = colname %||% scaled_name))
}


# Function to normalize protein abundances by subtracting sample means of reference proteins
center_to_reference <- function(
  df,
  df_reference,
  sample_name,
  abundance_column = "normalized_abundance"
) {
  # Step 1: Calculate sample means for reference proteins
  sample_means <- df_reference |>
    dplyr::group_by(!!rlang::sym(sample_name)) |>
    dplyr::summarise(
      reference_mean = mean(.data[[abundance_column]], na.rm = TRUE),
      reference_median = median(.data[[abundance_column]], na.rm = TRUE),
      .groups = "drop"
    )
  # Step 2: Join back to original data and subtract sample means
  normalized_df <-
    dplyr::left_join(df, sample_means, by = sample_name) |>
    dplyr::mutate(
      # Create normalized abundance column
      centered_abundance_by_mean = .data[[abundance_column]] - reference_mean,
      centered_abundance_by_median = .data[[abundance_column]] - reference_median
    ) |>
    dplyr::select(-reference_mean, -reference_median) # Remove the temporary column
  return(normalized_df)
}

#' center to reference
#'
#' takes the median of the lfqdareference per sample and subtracts it from a copy of lfqdata
#' @param lfqdata LFQData object containing the data to center
#' @param lfqdareference LFQData object containing the reference subset
#' @return The computed result.
#' @export
#' @examples
#' # example code
#'
#' bb <- sim_lfq_data_peptide_config(Nprot = 100)
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
