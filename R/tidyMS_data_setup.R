.validate_setup_file_name <- function(data, configuration) {
  if (is.null(configuration$file_name)) {
    abort_invalid_config("`file_name` column is not specified in the AnalysisConfiguration.")
  }
  if (!configuration$file_name %in% colnames(data)) {
    abort_missing_columns(configuration$file_name)
  }
}

.setup_hierarchy_columns <- function(data, configuration) {
  for (i in seq_along(configuration$hierarchy)) {
    data <- tidyr::unite(
      data,
      !!sym(configuration$hierarchy_keys()[i]),
      configuration$hierarchy[[i]],
      remove = FALSE,
      sep = configuration$sep
    )
  }
  dplyr::select(
    data,
    -dplyr::all_of(dplyr::setdiff(
      unlist(configuration$hierarchy),
      configuration$hierarchy_keys()
    ))
  )
}

.setup_factor_columns <- function(data, configuration) {
  if (length(configuration$factors) == 0) {
    abort_invalid_config(
      paste0(
        "No factors (explanatory variables) specified in the AnalysisConfiguration.\n",
        "Please use config$factors[\"Condition\"] = \"columnName\".\n",
        "where Condition is the new name of the variable and\n",
        "columnName is the name of the column containing the variable."
      )
    )
  }

  for (i in seq_along(configuration$factors)) {
    if (length(configuration$factors[[i]]) > 1) {
      data <- tidyr::unite(
        data,
        !!sym(configuration$factor_keys()[i]),
        configuration$factors[[i]],
        remove = FALSE,
        sep = configuration$sep
      )
    } else {
      newname <- configuration$factor_keys()[i]
      data <- dplyr::mutate(
        data,
        !!newname := as.character(!!sym(configuration$factors[[i]]))
      )
    }
  }
  data
}

.setup_sample_name <- function(data, configuration) {
  sample_name <- configuration$sample_name
  if (!sample_name %in% names(data)) {
    message("creating sampleName from file_name column")
    data[[sample_name]] <- tools::file_path_sans_ext(basename(data[[configuration$file_name]]))
    return(data)
  }

  message("column sampleName already exists, using :", sample_name)
  data
}

.add_setup_default_column <- function(data, column, value, warning_message) {
  if (column %in% colnames(data)) {
    return(data)
  }
  warning(warning_message)
  data[[column]] <- value
  data
}

.add_setup_default_columns <- function(data, configuration) {
  data <- .add_setup_default_column(
    data,
    configuration$isotope_label,
    "light",
    "no isotopeLabel column specified in the data, adding column isotopeLabel automatically and setting to 'light'."
  )
  data <- .add_setup_default_column(
    data,
    configuration$ident_q_value,
    0,
    "no qValue column specified in the data. Creating column qValue and setting qValues to 0."
  )
  .add_setup_default_column(
    data,
    configuration$nr_children,
    1,
    "no nr_children column specified in the data, adding column nr_children and setting to 1."
  )
}

.setup_duplicate_counts <- function(data, configuration) {
  data |>
    group_by(
      !!!syms(c(
        configuration$file_name,
        configuration$hierarchy_keys(),
        configuration$isotope_label
      ))
    ) |>
    summarize(n = n())
}

.setup_duplicate_message <- function(offending, configuration) {
  paste0(
    "setup_analysis: there is more than ONE observation for ",
    nrow(offending),
    " ",
    paste(configuration$hierarchy_keys(), collapse = ", "),
    " / sample (column '",
    configuration$sample_name,
    "', file '",
    configuration$file_name,
    "') combination(s).\n",
    "Each hierarchy key must be observed at most once per sample.\n",
    "First offending keys:\n",
    paste(utils::capture.output(utils::head(offending, 5)), collapse = "\n"),
    "\nRe-run setup_analysis(..., debug = TRUE) to return the full count table for inspection ",
    "(rows where n > 1)."
  )
}

.handle_setup_duplicates <- function(counts, configuration, debug) {
  offending <- dplyr::filter(dplyr::ungroup(counts), .data[["n"]] > 1)
  if (nrow(offending) == 0) {
    return(NULL)
  }

  message <- .setup_duplicate_message(offending, configuration)
  if (debug) {
    warning(message)
    return(counts)
  }
  rlang::abort(
    message,
    class = c("prolfqua_error_duplicate_keys", "prolfqua_error")
  )
}

#' Setup a tidy table compatible with a \code{\link{AnalysisConfiguration}}
#'
#' Extracts columns relevant for a configuration from a data frame
#' and create new columns e.g. sampleName column etc.
#' @param data data.frame
#' @param configuration AnalysisConfiguration
#' @param debug if FALSE (default) and any hierarchy-key/sample combination has
#'   more than one observation, stop with an informative error. If TRUE, warn
#'   and return the diagnostic count table (inspect rows where n > 1) instead.
#' @return A tibble with the columns relevant for \code{configuration}
#'   (sample, hierarchy, factor and response columns). Stops if any
#'   hierarchy-key/sample combination is observed more than once (unless
#'   \code{debug = TRUE}).
#' @export
#' @family configuration
#' @examples
#'
#' # Example with a normValue column (e.g., creatinine)
#' set.seed(1234)
#' data <- sim_lfq_data(Nprot = 10, PEPTIDE = TRUE, N = 4)
#' data$nr_children <- 1
#' data$isotopeLabel <- "light"
#' data$qValue <- 0
#'
#' # Add creatinine values per sample (not per protein/peptide)
#' sample_creatinine <- data |> dplyr::select(sample) |> dplyr::distinct() |>
#'   dplyr::mutate(creatinine = rnorm(dplyr::n(), mean = 100, sd = 10))
#' data <- dplyr::inner_join(data, sample_creatinine, by = "sample")
#'
#' config <- AnalysisConfiguration$new()
#' config$file_name = "sample"
#' config$factors["group_"] = "group"
#' config$hierarchy[["protein_Id"]] = c("proteinID", "idtype2")
#' config$hierarchy[["peptide_Id"]] = "peptideID"
#' config$set_response("abundance")
#' config$norm_value = "creatinine"
#'
#' adata <- setup_analysis(data, config)
#'
setup_analysis <- function(data, configuration, debug = FALSE) {
  configuration <- configuration$clone(deep = TRUE)
  .validate_setup_file_name(data, configuration)
  data <- .setup_hierarchy_columns(data, configuration)
  data <- .setup_factor_columns(data, configuration)
  data <- .setup_sample_name(data, configuration)
  data <- data |>
    dplyr::select(-dplyr::all_of(dplyr::setdiff(unlist(configuration$factors), configuration$factor_keys())))
  data <- .add_setup_default_columns(data, configuration)
  data <- data |> dplyr::select(c(configuration$id_vars(), configuration$value_vars()))

  duplicate_counts <- .setup_duplicate_counts(data, configuration)
  duplicate_result <- .handle_setup_duplicates(
    duplicate_counts,
    configuration,
    debug
  )
  if (!is.null(duplicate_result)) {
    return(duplicate_result)
  }
  data <- .complete_cases_impl(
    data,
    configuration$file_name,
    configuration$sample_name,
    configuration$factor_keys(),
    configuration$isotope_label,
    configuration$hierarchy_keys()
  )
  message("completing cases done")
  message("setup done")
  return(data)
}

#' Separates hierarchy columns into starting columns
#'
#'
#' @export
#' @param data data.frame
#' @param config AnalysisConfiguration
#' @family configuration
#' @keywords internal
#' @examples
#'
#'
#' bb <- sim_lfq_data_protein_config()
#' dt <- separate_hierarchy(bb$data, bb$config)
#' base::setdiff(colnames(dt) ,colnames(bb$data))
#' stopifnot(ncol(dt) >= ncol(bb$data))
#'
separate_hierarchy <- function(data, config) {
  for (hkey in config$hierarchy_keys_depth()) {
    if (length(config$hierarchy[[hkey]]) == 1 && hkey == config$hierarchy[[hkey]]) {} else {
      data <- data |> tidyr::separate(hkey, config$hierarchy[[hkey]], sep = config$sep, remove = FALSE)
    }
  }
  return(data)
}

# Internal implementation for complete_cases (used by setup_analysis before LFQData exists)
.complete_cases_impl <- function(pdata, file_name, sample_name, factor_keys, isotope_label, hierarchy_keys) {
  message("completing cases")
  fkeys <- unique(c(file_name, sample_name, factor_keys))
  hkeys <- c(isotope_label, hierarchy_keys)
  res <- tidyr::complete(
    pdata,
    tidyr::nesting(!!!syms(fkeys)),
    tidyr::nesting(!!!syms(hkeys))
  )
  return(res)
}

#' Complete cases
#'
#' The tidy table does not need to contain missing data.
#' This function re-establishes the missing observations in a sample.
#'
#' @export
#' @param lfqdata LFQData object
#' @keywords internal
#' @family configuration
#' @examples
#'
#' bb <- sim_lfq_data_protein_config()
#' lfq <- LFQData$new(bb$data, bb$config)
#' xx <- complete_cases(lfq)
#' stopifnot(nrow(bb$data) <= nrow(xx))
#'
complete_cases <- function(lfqdata) {
  .complete_cases_impl(
    lfqdata$data_long(),
    lfqdata$file_name(),
    lfqdata$sample_name(),
    lfqdata$factor_keys(),
    lfqdata$isotope_label(),
    lfqdata$hierarchy_keys()
  )
}


#' Sample subset of proteins/peptides/precursors
#' @param size size of sample
#' @param pdata tidy table
#' @param hierarchy_keys_depth character vector — hierarchy columns at current depth
#' @export
#' @keywords internal
#' @family configuration
#'
#' @examples
#' bb <- sim_lfq_data_peptide_config(Nprot = 5)
#' subset <- sample_subset(2, bb$data, bb$config$hierarchy_keys_depth())
#' length(unique(subset[[bb$config$hierarchy_keys_depth()[1]]]))
sample_subset <- function(size, pdata, hierarchy_keys_depth) {
  keys <- paste(hierarchy_keys_depth, collapse = ",")
  message("Sampling ", size, " ", keys)
  hkdf <- pdata |> select(all_of(hierarchy_keys_depth)) |> distinct() |> sample_n(size = size)
  sdata <- inner_join(hkdf, pdata)
  return(sdata)
}

#' Count the children of each parent in the hierarchy
#'
#' The number of distinct `level_b` entries, e.g. peptides, within each
#' `level_a` entry, e.g. protein.
#' @param pdata data.frame
#' @param level_a parent key columns
#' @param level_b child key columns
#' @param name name of the count column
#' @return data.frame of the `level_a` keys and the count column
#' @keywords internal
#' @family preprocessing
#' @examples
#' istar <- prolfqua::sim_lfq_data_peptide_config()
#' head(prolfqua:::nr_B_in_A(istar$data, "protein_Id", "peptide_Id"))
nr_B_in_A <- function(
  pdata,
  level_a,
  level_b,
  name = paste0("nr_", paste(level_b, collapse = "_"), "_IN_", paste(level_a, collapse = "_"))
) {
  dplyr::distinct(dplyr::select(pdata, dplyr::all_of(c(level_a, level_b)))) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(level_a))) |>
    dplyr::summarize(!!name := dplyr::n(), .groups = "drop")
}

#' Keep only the proteins with enough peptides
#'
#' Keeps the rows of the proteins with at least `config$min_peptides_protein`
#' children at the next hierarchy level, counted by [nr_B_in_A()], and adds
#' that count. The package
#' dataset `data_ionstar` calls it from its `filtered()` method.
#' @param pdata data.frame in the long format of `config`
#' @param config AnalysisConfiguration
#' @return list with the filtered `data` and the `name` of the count column,
#'   `NULL` when the hierarchy has no child level
#' @keywords internal
#' @family preprocessing
#' @examples
#' istar <- prolfqua::sim_lfq_data_peptide_config()
#' filtered <- prolfqua:::filter_proteins_by_peptide_count(istar$data, istar$config)
#' stopifnot(all(filtered$data[[filtered$name]] >= istar$config$min_peptides_protein))
filter_proteins_by_peptide_count <- function(pdata, config) {
  level_a <- config$hierarchy_keys_depth()
  level_b <- config$hierarchy_keys()[length(level_a) + 1]
  if (is.na(level_b)) {
    warning("here is no B in A")
    return(list(data = pdata, name = NULL))
  }
  counts <- nr_B_in_A(pdata, level_a, level_b)
  name <- setdiff(names(counts), level_a)
  data <- dplyr::inner_join(pdata, counts, by = level_a)
  list(data = dplyr::filter(data, .data[[name]] >= config$min_peptides_protein), name = name)
}
