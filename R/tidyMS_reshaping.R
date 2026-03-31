# Intensities to wide ----

#' Transform tidy table into a table with a column of responses for each sample
#'
#' @export
#' @keywords internal
tidy_to_wide <- function(data, rowIDs, columnLabels, value) {
  wide <- data |>
    dplyr::select(all_of(c(rowIDs, columnLabels, value)))

  wide_spread <- wide |>
    tidyr::pivot_wider(names_from = all_of(columnLabels), values_from = all_of(value))

  return(wide_spread)
}

#' transform long to wide
#' @export
#' @keywords internal
#' @return list with data, rowdata, and annotation (colData)
#' @examples
#'
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' config <- dd$config
#' data <- dd$data
#' res <- tidy_to_wide_config(data, config)
#' testthat::expect_equal(nrow(res$rowdata), nrow(res$data))
#' testthat::expect_equal(ncol(res$data) - ncol(res$rowdata) , nrow(res$annotation))
#' res <- tidy_to_wide_config(data, config, as.matrix = TRUE)
#' stopifnot(all(dim(res$data) == c(28,  12)))
#' stopifnot(all(dim(res$annotation) == c(12,  4)))
#' stopifnot(all(dim(res$rowdata) == c(28, 3)))
#'
#' res <- scale(res$data)
#' tidy_to_wide_config(data, config,  value = config$nr_children)
#'
#'
#' xt <- prolfqua::LFQData$new(dd$data, dd$config)
#' xt$data$nr_children
#' #xt$config$is_response_transformed <- TRUE
#' res <- xt$get_Aggregator("medpolish")
#' x <- res$aggregate()
#' towide <- tidy_to_wide_config(x$data, x$config,  value = x$config$nr_children)
#'
#' dd <- prolfqua::sim_lfq_data_protein_config()
#' dd$config$nr_children
#' dd$data
#' xt <- tidy_to_wide_config(dd$data, dd$config,  value = dd$config$nr_children)
#' xt$data
#'
tidy_to_wide_config <- function(
  data,
  config,
  as.matrix = FALSE,
  fileName = FALSE,
  sep = "~lfq~",
  value = config$get_response()
) {
  if (fileName) {
    newcolname <- config$file_name
  } else {
    newcolname <- config$sample_name
  }

  ids <- dplyr::select(
    data,
    all_of(c(config$sample_name, config$file_name, config$factor_keys(), config$isotope_label))
  ) |>
    dplyr::distinct() |>
    dplyr::arrange_at(newcolname)

  res <- tidy_to_wide(data, c(config$hierarchy_keys(), config$isotope_label), newcolname, value = value)
  rowdata <- res |> dplyr::select(all_of(c(config$hierarchy_keys(), config$isotope_label)))
  if (as.matrix) {
    resMat <- as.matrix(dplyr::select(res, -dplyr::all_of(c(config$hierarchy_keys(), config$isotope_label))))
    names <- rowdata |>
      tidyr::unite("newID", !!!dplyr::syms(c(config$hierarchy_keys(), config$isotope_label)), sep = sep) |>
      dplyr::pull("newID")
    rownames(resMat) <- names
    res <- resMat
  }
  return(list(data = res, annotation = ids, rowdata = rowdata))
}
