# Intensities to wide ----

#' transform long to wide
#' @param lfqdata LFQData object
#' @param as.matrix if TRUE return matrix, otherwise data.frame
#' @param file_name if TRUE use file_name as column labels, otherwise sample_name
#' @param sep separator for row IDs when as.matrix = TRUE
#' @param value column name to pivot (default: lfqdata$response())
#' @param row_keys columns identifying a row of the wide table (default: the
#'   hierarchy keys and the isotope label). They must identify one value per
#'   sample; e.g. drop the isotope label only for single-label data.
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
#' # row names without the isotope label
#' xh <- tidy_to_wide_config(lfqprot, as.matrix = TRUE, row_keys = lfqprot$hierarchy_keys())
#' stopifnot(!any(grepl("light", rownames(xh$data))))
#'
tidy_to_wide_config <- function(
  lfqdata,
  as.matrix = FALSE,
  file_name = FALSE,
  sep = "~lfq~",
  value = lfqdata$response(),
  row_keys = c(lfqdata$hierarchy_keys(), lfqdata$isotope_label())
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

  long <- dplyr::select(data, all_of(c(row_keys, newcolname, value)))
  # the default keys are unique by setup_analysis(); a reduced key set may not be
  if (
    !setequal(row_keys, c(lfqdata$hierarchy_keys(), lfqdata$isotope_label())) &&
      anyDuplicated(long[c(row_keys, newcolname)]) > 0
  ) {
    abort_bad_argument("row_keys", "identify one value per sample", not = paste(row_keys, collapse = ", "))
  }
  res <- long |>
    tidyr::pivot_wider(names_from = all_of(newcolname), values_from = all_of(value))
  sample_cols <- ids[[newcolname]]
  res <- res |>
    dplyr::select(dplyr::all_of(c(row_keys, sample_cols)))
  rowdata <- res |> dplyr::select(all_of(row_keys))
  if (as.matrix) {
    res_mat <- as.matrix(dplyr::select(res, -dplyr::all_of(row_keys)))
    names <- rowdata |>
      tidyr::unite("newID", !!!dplyr::syms(row_keys), sep = sep) |>
      dplyr::pull("newID")
    rownames(res_mat) <- names
    res <- res_mat
  }
  return(list(data = res, annotation = ids, rowdata = rowdata))
}
