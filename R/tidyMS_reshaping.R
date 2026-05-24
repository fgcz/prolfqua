# Intensities to wide ----

#' Transform tidy table into a table with a column of responses for each sample
#'
#' @export
#' @keywords internal
#' @examples
#' pdata <- data.frame(
#'   protein_Id = c("P1", "P1", "P2", "P2"),
#'   sampleName = c("S1", "S2", "S1", "S2"),
#'   abundance = c(10, 12, 20, 25)
#' )
#' tidy_to_wide(pdata, row_ids = "protein_Id", column_labels = "sampleName", value = "abundance")
tidy_to_wide <- function(data, row_ids, column_labels, value) {
  wide <- data |>
    dplyr::select(all_of(c(row_ids, column_labels, value)))

  wide_spread <- wide |>
    tidyr::pivot_wider(names_from = all_of(column_labels), values_from = all_of(value))

  return(wide_spread)
}

#' transform long to wide
#' @param lfqdata LFQData object
#' @param as.matrix if TRUE return matrix, otherwise data.frame
#' @param file_name if TRUE use file_name as column labels, otherwise sample_name
#' @param sep separator for row IDs when as.matrix = TRUE
#' @param value column name to pivot (default: lfqdata$response())
#' @export
#' @keywords internal
#' @return list with data, rowdata, and annotation (colData)
#' @examples
#'
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' lfqdata <- prolfqua::LFQData$new(dd$data, dd$config)
#' res <- tidy_to_wide_config(lfqdata)
#' testthat::expect_equal(nrow(res$rowdata), nrow(res$data))
#' testthat::expect_equal(ncol(res$data) - ncol(res$rowdata) , nrow(res$annotation))
#' res <- tidy_to_wide_config(lfqdata, as.matrix = TRUE)
#' stopifnot(all(dim(res$data) == c(28,  12)))
#' stopifnot(all(dim(res$annotation) == c(12,  4)))
#' stopifnot(all(dim(res$rowdata) == c(28, 3)))
#'
#' res <- scale(res$data)
#' tidy_to_wide_config(lfqdata, value = lfqdata$nr_children_col())
#'
#' res <- lfqdata$get_Aggregator("medpolish")
#' x <- res$aggregate()
#' towide <- tidy_to_wide_config(x, value = x$nr_children_col())
#'
#' dd <- prolfqua::sim_lfq_data_protein_config()
#' lfqprot <- prolfqua::LFQData$new(dd$data, dd$config)
#' xt <- tidy_to_wide_config(lfqprot, value = lfqprot$nr_children_col())
#' xt$data
#'
tidy_to_wide_config <- function(
  lfqdata,
  as.matrix = FALSE,
  file_name = FALSE,
  sep = "~lfq~",
  value = lfqdata$response()
) {
  data <- lfqdata$data_long()
  if (file_name) {
    newcolname <- lfqdata$file_name()
  } else {
    newcolname <- lfqdata$sample_name()
  }

  ids <- dplyr::select(
    data,
    all_of(c(lfqdata$sample_name(), lfqdata$file_name(), lfqdata$factor_keys(), lfqdata$isotope_label()))
  ) |>
    dplyr::distinct() |>
    dplyr::arrange_at(newcolname)

  hierarchy_isotope <- c(lfqdata$hierarchy_keys(), lfqdata$isotope_label())
  res <- tidy_to_wide(data, hierarchy_isotope, newcolname, value = value)
  sample_cols <- ids[[newcolname]]
  res <- res |>
    dplyr::select(dplyr::all_of(c(hierarchy_isotope, sample_cols)))
  rowdata <- res |> dplyr::select(all_of(hierarchy_isotope))
  if (as.matrix) {
    res_mat <- as.matrix(dplyr::select(res, -dplyr::all_of(hierarchy_isotope)))
    names <- rowdata |>
      tidyr::unite("newID", !!!dplyr::syms(hierarchy_isotope), sep = sep) |>
      dplyr::pull("newID")
    rownames(res_mat) <- names
    res <- res_mat
  }
  return(list(data = res, annotation = ids, rowdata = rowdata))
}
