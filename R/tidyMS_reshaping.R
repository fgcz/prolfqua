# Intensities to wide ----

#' Transform long to wide format
#' @param lfqdata [LFQData] object with the data and its configuration
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
    dplyr::arrange(dplyr::pick(dplyr::all_of(newcolname)))

  long <- dplyr::select(data, all_of(c(row_keys, newcolname, value)))
  # the default keys are unique by setup_analysis(); a reduced key set may not be
  if (
    !setequal(row_keys, c(lfqdata$hierarchy_keys(), lfqdata$isotope_label())) &&
      anyDuplicated(long[c(row_keys, newcolname)]) > 0
  ) {
    abort_bad_argument("row_keys", "identify one value per sample", not = toString(row_keys))
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

#' Input checks of lfqdata_from_wide(): `values` is a named list of matrices of one shape.
#' @param values named list of feature-by-sample matrices with matching dimensions and column names
#' @noRd
.check_wide_values <- function(values) {
  is_named_matrix_list <- is.list(values) &&
    !is.null(names(values)) &&
    all(nzchar(names(values))) &&
    all(vapply(values, is.matrix, logical(1)))
  if (!is_named_matrix_list) {
    abort_bad_argument("values", "be a matrix or a named list of matrices")
  }
  first <- values[[1]]
  same_shape <- vapply(
    values,
    function(m) identical(dim(m), dim(first)) && identical(colnames(m), colnames(first)),
    logical(1)
  )
  if (!all(same_shape)) {
    abort_bad_argument("values", "hold matrices with identical dimensions and column names")
  }
  invisible(TRUE)
}

#' Input checks of lfqdata_from_wide(): rowdata and annotation match the matrices.
#' @param values named list of feature-by-sample matrices with matching dimensions and column names
#' @param rowdata one row of identifiers per feature
#' @param annotation one row of annotations per sample
#' @param config AnalysisConfiguration defining the column roles
#' @noRd
.check_wide_input <- function(values, rowdata, annotation, config) {
  .check_wide_values(values)
  first <- values[[1]]
  missing_rowdata <- setdiff(config$hierarchy_keys(), colnames(rowdata))
  if (length(missing_rowdata) > 0) {
    abort_missing_columns(missing_rowdata, data_nm = "rowdata")
  }
  if (!config$sample_name %in% colnames(annotation)) {
    abort_missing_columns(config$sample_name, data_nm = "annotation")
  }
  if (nrow(rowdata) != nrow(first)) {
    abort_bad_argument("rowdata", "have one row per matrix row", not = paste(nrow(rowdata), "rows"))
  }
  samples <- as.character(annotation[[config$sample_name]])
  if (anyDuplicated(samples) > 0 || !setequal(colnames(first), samples)) {
    abort_bad_argument("annotation", "have one row per matrix column, named by `config$sample_name`")
  }
  invisible(TRUE)
}

#' One feature x sample matrix to long format, keyed by row position and sample.
#' @param mat feature-by-sample matrix
#' @param value name of the output measurement column
#' @param sample_name column identifying samples
#' @noRd
.wide_matrix_to_long <- function(mat, value, sample_name) {
  long <- tibble::as_tibble(mat, .name_repair = "minimal")
  long$.row <- seq_len(nrow(long))
  tidyr::pivot_longer(long, -".row", names_to = sample_name, values_to = value)
}

#' Build LFQData from wide matrices
#'
#' The inverse of [tidy_to_wide_config()]: feature x sample matrices, their row
#' annotation and their sample annotation back to long-format LFQData. Row `i`
#' of every matrix belongs to row `i` of `rowdata`; matrix column names are the
#' sample names.
#'
#' @param values a matrix (features in rows, samples in columns), stored under
#'   `config$get_response()`, or a named list of such matrices with identical
#'   dimensions and column names; list names are the value columns
#' @param rowdata data.frame with one row per matrix row and the hierarchy key
#'   columns; the isotope label column is optional and set to `"light"` if absent
#' @param annotation data.frame with one row per sample: `config$sample_name`
#'   plus the file name and factor columns
#' @param config AnalysisConfiguration describing the result
#' @param prefix prefix passed to [LFQData]
#' @return LFQData
#' @export
#' @family LFQData
#' @examples
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' lfqdata <- prolfqua::LFQData$new(dd$data, dd$config)
#' wide <- tidy_to_wide_config(lfqdata, as.matrix = TRUE)
#' back <- lfqdata_from_wide(wide$data, wide$rowdata, wide$annotation, lfqdata$get_config())
#' stopifnot(nrow(back$data_long()) == nrow(lfqdata$data_long()))
#' @srrstatsTODO {G2.7} Long data frames and wide matrices are supported; limit: no SummarizedExperiment input adapter.
#' @srrstatsTODO {G2.10} Most extraction uses names and [[; limit: cross-class data.frame/tibble/data.table tests are
#'   incomplete.
lfqdata_from_wide <- function(values, rowdata, annotation, config, prefix = "ms") {
  if (is.matrix(values)) {
    values <- stats::setNames(list(values), config$get_response())
  }
  .check_wide_input(values, rowdata, annotation, config)
  sample_name <- config$sample_name
  samples <- as.character(annotation[[sample_name]])

  long <- Reduce(
    function(acc, one) dplyr::inner_join(acc, one, by = c(".row", sample_name), relationship = "one-to-one"),
    Map(.wide_matrix_to_long, values, names(values), sample_name)
  )
  rowdata <- tibble::as_tibble(rowdata)
  rowdata$.row <- seq_len(nrow(rowdata))
  long <- dplyr::inner_join(rowdata, long, by = ".row", relationship = "one-to-many")
  long$.row <- NULL

  annotation <- tibble::as_tibble(annotation)
  annotation[[sample_name]] <- samples
  # columns already taken from rowdata (e.g. the isotope label) are not joined twice
  annotation <- annotation[, c(sample_name, setdiff(colnames(annotation), colnames(long))), drop = FALSE]
  long <- dplyr::inner_join(long, annotation, by = sample_name, relationship = "many-to-one")
  if (!config$isotope_label %in% colnames(long)) {
    long[[config$isotope_label]] <- "light"
  }
  LFQData$new(long, config, prefix = prefix)
}
