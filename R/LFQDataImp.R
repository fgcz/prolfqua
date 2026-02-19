


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
#' impute using zCompositions
#' @param lfqdata LFQData object containing the data to impute
#' @param method imputation method passed to zCompositions (default "multRepl")
#' @param lod limit of detection strategy, either "global" or "quantile"
#' @export
#'
#' @examples
#' # example code
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' lfqdata <- LFQData$new(dd$data, dd$config)
#'

impute_with_zcomp <- function(lfqdata,
                              method = c("multRepl","GBM","SQ","BL","CZM" ) ,
                              lod = c("global","quantile")) {

  lod <- match.arg(lod)
  method <- match.arg(method)
  method <- "SQ"
  lod <- "global"
  wide <- lfqdata$to_wide(as.matrix = TRUE)

  if (!any(is.na(wide$data))) {
    return(t(wide$data))
  }
  # TODO: check
  data_matrix <- t(wide$data)
  data_matrix <- wide$data

  #if (method %in% c("lrEM", "lrDA")) {
  dl <- sapply(function_lod_quantile(data_matrix), median)

  if (lod == "global") {
    lod_value <- estimate_lod_global(data_matrix)
    if (length(lod_value) > 0) {
      dl <- rep(median(lod_value), ncol(data_matrix))
    }
  }

  data_matrix[is.na(data_matrix)] <- NA

  # Pre-process for methods that need complete data
  if (method == "multRepl") {
    imputed_data <- zCompositions::multRepl(data_matrix, dl = dl, label = NA)
  } else if (method %in% c("GBM","SQ","BL","CZM")) {
    imputed_data <- zCompositions::cmultRepl(data_matrix, label = NA, method = method)
  }


  lfqdata$data <- response_matrix_as_tibble(
    imputed_data,
    paste0(lfqdata$response(),"_imputed"),
    lfqdata$conf, lfqdata$data)

  return(imputed_data)


}


# Example usage:
#
# # Simple usage - no detection limits needed for GS method
# imputed_matrix <- impute_with_zcomp(data, config, method = "GS")
#
# # lrEM method with auto-estimated detection limits
# imputed_matrix <- impute_with_zcomp(data, config, method = "lrEM")
#
# # Custom detection limit estimation
# imputed_matrix <- impute_with_zcomp(data, config, method = "lrEM",
#                                     dl_method = "quantile")
#
# # Provide your own detection limits
# custom_dl <- c(0.01, 0.02, 0.005, ...)  # One per variable
# imputed_matrix <- impute_with_zcomp(data, config, method = "lrEM", dl = custom_dl)
#
# # For quick testing without worrying about detection limits:
# imputed_matrix <- impute_with_zcomp(data, config, method = "GS", verbose = FALSE)

# LFQDataImp ----
#'
#' Decorates LFQData with methods to impute missing values
#'
#' @export
#' @family LFQData
#'
LFQDataImp <- R6::R6Class(
  "LFQDataImp",
  public = list(
    #' @field lfq LFQData
    lfq = NULL,
    #' @description
    #' initialize
    #' @param lfq lfqdata
    #' @param prefix datasetname
    initialize = function(lfq, prefix = "protein"){
      self$lfq=lfq$get_copy()
      self$prefix = prefix
    })
)
