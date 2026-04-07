# LFQDataTransformer ----

#' Decorate LFQData with Methods for transforming Intensities
#'
#' @export
#'
#' @examples
#'
#' istar <- prolfqua_data('data_ionstar')$filtered()
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(data, istar$config)
#'
#' lfqcopy <- lfqdata$get_copy()
#' lfqTrans <- lfqcopy$get_Transformer()
#'
#' x <- lfqTrans$intensity_array(log2)
#'
#' x$lfq$config$is_response_transformed
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
      self$lfq = lfqdata$clone(deep = TRUE)
    },
    #' @description
    #' log2 transform data
    #' @param force if `FALSE`, data already log2 transformed will not be
    #'   transformed a second time. `TRUE` forces log transformation
    #' @return LFQDataTransformer
    log2 = function(force = FALSE) {
      if (self$lfq$is_transformed() == FALSE | force) {
        self$lfq$data <- prolfqua::transform_work_intensity(self$lfq$data, self$lfq$config, log2)
        self$lfq$is_transformed(TRUE)
      } else {
        warning("data already transformed. If you still want to log2 tranform, set force = TRUE")
      }
      invisible(self)
    },
    #' @description
    #' get mean and variance and standard deviation in each sample
    #' @return list with means and mads
    get_scales = function() {
      get_robscales(self$lfq$data, self$lfq$config)
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
        invisible(NULL)
      }
      scales <- prolfqua::scale_with_subset(
        self$lfq$data,
        lfqsubset$data,
        self$lfq$config,
        preserve_mean = preserve_mean
      )
      self$lfq$data <- scales$data
      if (!is.null(colname)) {
        self$lfq$data <- self$lfq$data |>
          dplyr::rename(!!colname := self$lfq$config$get_response())
        self$lfq$config$pop_response()
        self$lfq$config$set_response(colname)
      }
      invisible(self)
    },
    #' @description
    #' log2 transform and robust scale data based on subset
    #' @param lfqsubset LFQData subset to use for normalization
    #' @param preserve_mean should original mean value be preserved TRUE, if FALSE then center at zero
    #' @param colname - how to name the transformed intensities, default transformedIntensity
    #' @return LFQDataTransformer (self)
    #'
    center_to_reference = function(lfqsubset) {
      message("data is transoformed: ", self$lfq$is_transformed())
      if (self$lfq$is_transformed() != lfqsubset$is_transformed()) {
        stop("the subset must have the same config as self")
      }
      if (!self$lfq$is_transformed()) {
        warning("data should be log2 transformed")
      }
      center_to_reference_cfg(self$lfq, lfqsubset, copy = FALSE)
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
      if (!self$lfq$is_transformed() | force) {
        .call <- as.list(match.call())
        r <- prolfqua::transform_work_intensity(
          self$lfq$data,
          self$lfq$config,
          .func = .func,
          .funcname = deparse(.call$.func)
        )
        self$lfq$data <- r
        self$lfq$is_transformed(TRUE)
      } else {
        warning("data already transformed. If you still want to log2 tranform, set force = TRUE")
      }
      invisible(self)
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
      if (!self$lfq$is_transformed() | force) {
        .call <- as.list(match.call())
        r <- prolfqua::apply_to_response_matrix(
          self$lfq$data,
          self$lfq$config,
          .func = .func,
          .funcname = deparse(.call$.func)
        )
        self$lfq$data <- r
        self$lfq$is_transformed(TRUE)
      } else {
        warning("data already transformed. If you still want to log2 tranform, set force = TRUE")
      }
      invisible(self)
    }
  )
)

# Intensity transformation helpers ----

#' Transform intensity
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param .func function to transform intensities e.g. log2
#' @param .funcname generates new name from name of transformation and old working intensity column name.
#' @param intensity_new_name column name for new intensity, default NULL
#' @return data.frame
#' @export
#' @keywords internal
#' @examples
#'
#' dd <- prolfqua_data('data_spectronautDIA250_A')
#' config <- dd$config_f()
#' analysis <- dd$analysis(dd$data,config)
#' x <- transform_work_intensity(analysis, config, .func = log2)
#' stopifnot("log2_FG.Quantity" %in% colnames(x))
#' config <- dd$config_f()
#' x <- transform_work_intensity(analysis, config, .func = asinh)
#' stopifnot("asinh_FG.Quantity" %in% colnames(x))
#'
transform_work_intensity <- function(pdata, config, .func, .funcname = NULL, intensity_new_name = NULL, deep = FALSE) {
  if (deep) {
    config <- config$clone(deep = TRUE)
  }
  .call <- as.list(match.call())

  if (is.null(intensity_new_name)) {
    .funcname <- if (is.null(.funcname)) {
      deparse(.call$.func)
    } else {
      .funcname
    }
    newcol <- paste(.funcname, config$get_response(), sep = "_")
  } else {
    newcol <- intensity_new_name
  }

  #pdata <- pdata |> dplyr::mutate(across(all_of(config$get_response()),
  #                                    .fns = list(!!newcol := .func)))
  response_col <- config$get_response()
  vals <- pdata[[response_col]]
  if (identical(.func, log2) || identical(.func, log) || identical(.func, log10)) {
    n_zero <- sum(vals == 0, na.rm = TRUE)
    n_neg <- sum(vals < 0, na.rm = TRUE)
    if (n_neg > 0) {
      warning(
        "log transform: ",
        n_neg,
        " negative values in '",
        response_col,
        "' will produce NaN. Consider filtering or using asinh."
      )
    }
    if (n_zero > 0) {
      warning(
        "log transform: ",
        n_zero,
        " zeros in '",
        response_col,
        "' will produce -Inf. Consider replacing zeros with NA first."
      )
    }
  }

  pdata <- pdata |> dplyr::mutate(!!sym(newcol) := .func(!!sym(response_col)))

  config$set_response(newcol)
  message("Column added : ", newcol)
  config$is_response_transformed <- TRUE

  if (deep) {
    return(list(data = pdata, config = config))
  } else {
    return(pdata)
  }
}

#' Takes matrix of responses and converts into tibble
#'
#' @param pdata (matrix)
#' @param value name of column to store values in. (see `gather`)
#' @param config AnalysisConfiguration
#' @param data lfqdata
#' @param sep separater to unite the hierarchy keys.
#' @export
#'
#' @keywords internal
#' @examples
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' data <- dd$data
#' conf <- dd$config
#' res <- tidy_to_wide_config(data, conf, as.matrix = TRUE)
#'
#' res <- scale(res$data)
#' xx <- response_matrix_as_tibble(res,"srm_intensityScaled", conf)
#' xx <- response_matrix_as_tibble(res,"srm_intensityScaled", conf, data)
#' conf$get_response() == "srm_intensityScaled"
#'
response_matrix_as_tibble <- function(pdata, value, config, data = NULL, sep = "~lfq~") {
  pdata <- dplyr::bind_cols(
    tibble::tibble("row.names" := rownames(pdata)),
    tibble::as_tibble(pdata)
  )
  pdata <- tidyr::pivot_longer(pdata, cols = -1, names_to = config$sample_name, values_to = value)
  pdata <- tidyr::separate(
    pdata,
    "row.names",
    unique(c(config$hierarchy_keys(), config$isotope_label)),
    sep = sep
  )
  if (!is.null(data)) {
    pdata <- dplyr::inner_join(data, pdata)
    config$set_response(value)
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
#' conf <- bb$config
#' sample_analysis <- bb$data
#' pepIntensityNormalized <- transform_work_intensity(sample_analysis, conf, log2)
#' s1 <- get_robscales(pepIntensityNormalized, conf)
#'
#' res <- scale_with_subset(pepIntensityNormalized, pepIntensityNormalized, conf)
#' s2 <- get_robscales(res$data, conf)
#' abs(mean(s1$mads) - mean(s2$mads)) < 0.1
#'
#'
get_robscales <- function(data, config) {
  data <- tidy_to_wide_config(data, config, as.matrix = TRUE)$data
  scales <- .get_robscales(data)
  return(scales)
}


#' robust scale wrapper
#' @keywords internal
#' @family preprocessing
#' @export
robust_scale <- function(data, dim = 2, preserve_mean = FALSE) {
  scales <- .get_robscales(data, dim = dim)
  data <- sweep(data, dim, scales$medians, "-")
  if (!any(scales$mads == 0)) {
    mads <- scales$mads / mean(scales$mads)
    data = sweep(data, dim, mads, "/")
  } else {
    warning("SKIPPING scaling step in robust_scale: one or more MAD values are zero.")
  }
  meanmed <- mean(scales$medians)
  addmean <- if (preserve_mean) {
    meanmed
  } else {
    0
  }
  return(data + addmean)
}


#' Apply function requiring a matrix to tidy table
#'
#' @param data data.frame
#' @param config AnalysisConfiguration
#' @param .func function
#' @param .funcname name of function (used for creating new column)
#' @export
#' @keywords internal
#' @family preprocessing
#' @examples
#'
#'
#' bb <- sim_lfq_data_peptide_config(Nprot = 100)
#' data <- bb$data
#' conf <- bb$config
#' res <- apply_to_response_matrix(data, conf, .func = base::scale)
#'
#' stopifnot("abundance_base..scale" %in% colnames(res))
#' stopifnot("abundance_base..scale" == conf$get_response())
#' conf <- bb$config$clone(deep=TRUE)
#' conf$work_intensity <- "abundance"
#' res <- apply_to_response_matrix(data, conf$clone(deep=TRUE), .func = robust_scale)
#'
#' # Normalize data using the vsn method from bioconductor
#'
#' if( require("vsn")){
#'  res <- apply_to_response_matrix(data, conf$clone(deep=TRUE), .func = vsn::justvsn)
#' }
#'
apply_to_response_matrix <- function(data, config, .func, .funcname = NULL) {
  .call <- as.list(match.call())
  .funcname <- if (is.null(.funcname)) {
    deparse(.call$.func)
  } else {
    .funcname
  }
  colname <- make.names(paste(config$get_response(), .funcname, sep = "_"))
  mat <- tidy_to_wide_config(data, config, as.matrix = TRUE)$data
  mat <- .func(mat)
  data <- response_matrix_as_tibble(mat, colname, config, data)
  return(data)
}

#' Scale data using a subset of the data
#'
#' this should reduce the overall variance.
#'
#' @export
#' @keywords internal
#' @param data the whole dataset
#' @param subset a subset of the dataset
#' @param config configuration
#' @param preserve_mean default FALSE - sets mean to zero
#' @param get_scales return a list of transformed data and the scaling parameters
#' @family preprocessing
#' @examples
#'
#'
#'
#' bb <-sim_lfq_data_peptide_config(Nprot = 100)
#' conf <- bb$config$clone(deep=TRUE)
#' sample_analysis <- bb$data
#'
#' res <- transform_work_intensity(sample_analysis, conf, log2)
#' s1 <- get_robscales(res, conf)
#' res <- scale_with_subset(res, res, conf)
#' s2 <- get_robscales(res$data, conf)
#' stopifnot(abs(mean(s1$mads) - mean(s2$mads)) < 1e-6)
scale_with_subset <- function(data, subset, config, preserve_mean = FALSE, get_scales = TRUE) {
  colname <- make.names(paste(config$get_response(), "subset_scaled", sep = "_"))
  subset <- tidy_to_wide_config(subset, config, as.matrix = TRUE)$data

  scales <- .get_robscales(subset)
  mat <- tidy_to_wide_config(data, config, as.matrix = TRUE)$data
  mat <- sweep(mat, 2, scales$medians, "-")
  if (!any(scales$mads == 0)) {
    mads <- scales$mads / mean(scales$mads)
    mat <- sweep(mat, 2, mads, "/")
  } else {
    warning("SKIPPING scaling step in scale_with_subset function.")
  }

  meanmed <- mean(scales$medians)
  addmean <- if (preserve_mean) {
    meanmed
  } else {
    0
  }
  mat <- mat + addmean
  data <- response_matrix_as_tibble(mat, colname, config, data)
  if (get_scales) {
    return(list(data = data, scales = scales))
  } else {
    return(data)
  }
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
#' takes the mean or median of the lfqdareference per sample and subtracts from lfqdata
#' @param lfqdata LFQData object containing the data to center
#' @param lfqdareference LFQData object containing the reference subset
#' @param summary character, summary statistic to use ("median" or "mean")
#' @param copy logical, if TRUE return a copy, otherwise modify in place
#' @export
#' @examples
#' # example code
#'
#' bb <- sim_lfq_data_peptide_config(Nprot = 100)
#' x <- LFQData$new(bb$data, bb$config)
#' xc <- x$get_copy()
#' xc$data <- xc$data |> dplyr::filter(protein_Id == "0EfVhX~3967")
#' xxd <- center_to_reference_cfg(x, xc, summary="median")
#' xxd$response()
#' xxd$data
#' center_to_reference_cfg(x, xc, summary="median", copy=FALSE)
#' x$response()
#'
center_to_reference_cfg <- function(lfqdata, lfqdareference, summary = c("median", "mean"), copy = TRUE) {
  summary <- match.arg(summary)
  if (copy) {
    resdata <- lfqdata$get_copy()
  } else {
    resdata <- lfqdata
  }
  cfg <- resdata$config
  data <- center_to_reference(lfqdata$data, lfqdareference$data, cfg$sample_name, cfg$get_response())
  resdata$data <- data
  if (summary == "median") {
    cfg$set_response("centered_abundance_by_median")
  } else if (summary == "mean") {
    cfg$set_response("centered_abundance_by_mean")
  }
  invisible(resdata)
}
