


#' esitmate lod
#' @param data_matrix numeric matrix of abundance values
#' @param prop_na numeric, percentage threshold for NA proportion per row
#' @export
#' @examples
#' # example code
#' istar <- sim_lfq_data_peptide_config()
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' xx <- lfqdata$to_wide(as.matrix=TRUE)
#' stopifnot(length(estimate_lod_global(xx$data, prop_na = 90)) == 0)
#' stopifnot(length(estimate_lod_global(xx$data, prop_na = 10)) > 0)
estimate_lod_global <- function(data_matrix, prop_na = 90) {
  frac <- ceiling(ncol(data_matrix) / (100/prop_na))
  xx <- na.omit(as.numeric(data_matrix[apply(data_matrix, 1, function(x){sum(is.na(x))}) >= frac,]))
  return(xx)
}


#' get smallest values per sample
#' @param data_matrix numeric matrix or data.frame of abundance values
#' @param percent numeric, percentile cutoff for smallest values
#' @export
#' @examples
#' # example code
#'
#' istar <- sim_lfq_data_peptide_config()
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' xx <- lfqdata$to_wide(as.matrix=TRUE)
#' s <- function_lod_quantile(xx$data)
#' sapply(s, median)
#' sapply(s, mean)
function_lod_quantile <- function(data_matrix, percent = 10) {
  # Convert to data frame if it's a matrix
  if (is.matrix(data_matrix)) {
    data_matrix <- as.data.frame(data_matrix)
  }

  # Initialize list to store results
  result_list <- list()

  # Process each column
  for (col_name in names(data_matrix)) {
    col_data <- as.numeric(data_matrix[[col_name]])

    # Remove NA values
    col_data <- na.omit(col_data)

    if (length(col_data) == 0) {
      result_list[[col_name]] <- numeric(0)
      next
    }

    # Calculate number of values to select (x% of total)
    n_select <- ceiling(length(col_data) * percent / 100)

    # Sort and select the smallest values
    sorted_values <- sort(col_data)
    result_list[[col_name]] <- sorted_values[1:n_select]
  }

  return(result_list)
}
#' Impute missing values using zCompositions
#'
#' Replaces missing values in an LFQData object using methods from the
#' zCompositions package. The limit of detection (LOD) can be estimated
#' globally or per-sample using quantiles.
#'
#' @param lfqdata LFQData object containing the data to impute
#' @param method imputation method: "multRepl" (multiplicative replacement),
#'   "GBM", "SQ", "BL", or "CZM" (passed to zCompositions)
#' @param lod limit of detection strategy, either "global" or "quantile"
#' @return the modified LFQData object (lfqdata), with imputed values
#' @note This function assumes that missing values are Missing Completely
#'   At Random (MCAR) or Missing At Random (MAR). If missingness is
#'   abundance-dependent (MNAR, common in proteomics DDA), the imputed
#'   values and downstream statistics may be biased. For MNAR-aware
#'   analysis, consider packages such as proDA or DEqMS.
#' @export
#'
#' @examples
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' lfqdata <- LFQData$new(dd$data, dd$config)
#' if (requireNamespace("zCompositions", quietly = TRUE)) {
#'   wide_before <- lfqdata$to_wide(as.matrix = TRUE)
#'   has_na_before <- any(is.na(wide_before$data))
#'   lfqdata <- impute_with_zcomp(lfqdata, method = "multRepl", lod = "global")
#'   wide_after <- lfqdata$to_wide(as.matrix = TRUE)
#'   has_na_after <- any(is.na(wide_after$data))
#'   stopifnot(has_na_before || !has_na_after)
#'   stopifnot(!has_na_after)
#' }
#'
impute_with_zcomp <- function(lfqdata,
                              method = c("multRepl","GBM","SQ","BL","CZM") ,
                              lod = c("global","quantile")) {

  lod <- match.arg(lod)
  method <- match.arg(method)
  wide <- lfqdata$to_wide(as.matrix = TRUE)

  if (!any(is.na(wide$data))) {
    return(lfqdata)
  }

  data_matrix <- wide$data
  dl <- sapply(function_lod_quantile(data_matrix), median)

  if (lod == "global") {
    lod_value <- estimate_lod_global(data_matrix)
    if (length(lod_value) > 0) {
      dl <- rep(median(lod_value), ncol(data_matrix))
    }
  }

  if (method == "multRepl") {
    imputed_data <- zCompositions::multRepl(data_matrix, dl = dl, label = NA)
  } else if (method %in% c("GBM","SQ","BL","CZM")) {
    imputed_data <- zCompositions::cmultRepl(data_matrix, label = NA, method = method)
  }

  lfqdata$data <- response_matrix_as_tibble(
    imputed_data,
    paste0(lfqdata$response(),"_imputed"),
    lfqdata$config, lfqdata$data)

  return(lfqdata)
}

# LFQDataImp ----
#'
#' Decorates LFQData with methods to impute missing values
#'
#' Use `lfqdata$get_Imputer()` to create an instance. The imputer works on
#' a deep copy of the data, so the original LFQData is not modified.
#'
#' Typical workflow:
#' \preformatted{
#' imp <- lfqdata$get_Imputer()
#' imp$impute(method = "multRepl", lod = "global")
#' imp$lfq$get_Plotter()$pca()
#' }
#'
#' @export
#' @family LFQData
#'
#' @examples
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' lfqdata <- LFQData$new(dd$data, dd$config)
#' imp <- LFQDataImp$new(lfqdata)
#' if (requireNamespace("zCompositions", quietly = TRUE)) {
#'   imp$impute(method = "multRepl", lod = "global")
#'   wide <- imp$lfq$to_wide(as.matrix = TRUE)
#'   stopifnot(!any(is.na(wide$data)))
#' }
#'
LFQDataImp <- R6::R6Class(
  "LFQDataImp",
  public = list(
    #' @field lfq LFQData
    lfq = NULL,
    #' @description
    #' initialize
    #' @param lfqdata LFQData object
    initialize = function(lfqdata){
      self$lfq = lfqdata$clone(deep = TRUE)
    },
    #' @description
    #' Impute missing values using zCompositions
    #' @param method imputation method (default "multRepl")
    #' @param lod limit of detection strategy (default "global")
    #' @return invisible(self) for chaining
    impute = function(method = c("multRepl","GBM","SQ","BL","CZM"),
                      lod = c("global","quantile")) {
      self$lfq <- impute_with_zcomp(self$lfq, method = method, lod = lod)
      invisible(self)
    })
)
