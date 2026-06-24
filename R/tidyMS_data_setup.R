#' Setup a tidy table compatible with a \code{\link{AnalysisConfiguration}}
#'
#' Extracts columns relevant for a configuration from a data frame
#' and create new columns e.g. sampleName column etc.
#' @param data data.frame
#' @param configuration AnalysisConfiguration
#' @param cc complete cases default TRUE
#' @param from_factors if TRUE, create sampleName from factor columns
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
#' skylineconfig <- AnalysisConfiguration$new()
#' skylineconfig$file_name <- "Replicate.Name"
#' skylineconfig$hierarchy[["protein_Id"]] <- "Protein.Name"
#' skylineconfig$hierarchy[["peptide_Id"]] <- "Peptide.Sequence"
#' skylineconfig$hierarchy[["precursor_Id"]] <- c("Peptide.Sequence", "Precursor.Charge")
#' skylineconfig$hierarchy[["fragment_Id"]] <- c(
#'   "Peptide.Sequence", "Precursor.Charge",
#'   "Fragment.Ion", "Product.Charge"
#' )
#' skylineconfig$ident_q_value <- "Detection.Q.Value"
#' skylineconfig$set_response("Area")
#' skylineconfig$isotope_label <- "Isotope.Label.Type"
#' skylineconfig$factors[["Time"]] = "Sampling.Time.Point"
#' sample_analysis <- setup_analysis(prolfqua_data('data_skylinePRMSample_A')$data, skylineconfig)
#'
#' # Example with normValue column (e.g., creatinine)
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
setup_analysis <- function(data, configuration, cc = TRUE, from_factors = FALSE, debug = FALSE) {
  configuration <- configuration$clone(deep = TRUE)
  if (is.null(configuration$file_name)) {
    stop("file_name column is not specified in configuration.")
  }
  if (!configuration$file_name %in% colnames(data)) {
    stop("File name column :", configuration$file_name, ", is missing in data.")
  }

  # extract hierarchy columns
  for (i in seq_along(configuration$hierarchy)) {
    data <- tidyr::unite(
      data,
      !!sym(configuration$hierarchy_keys()[i]),
      configuration$hierarchy[[i]],
      remove = FALSE,
      sep = configuration$sep
    )
  }
  data <- dplyr::select(
    data,
    -dplyr::all_of(dplyr::setdiff(unlist(configuration$hierarchy), configuration$hierarchy_keys()))
  )

  # extract factors
  if (length(configuration$factors) == 0) {
    stop(
      "No factors (explanatory variables) specified in the AnalysisConfiguration.\n",
      "Please use config$factors[\"Condition\"] = \"columnName\".\n",
      "where Condition is the new name of the variable and\n",
      "columnName is the name of the column containing the variable."
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
      data <- dplyr::mutate(data, !!newname := as.character(!!sym(configuration$factors[[i]])))
    }
  }

  sample_name <- configuration$sample_name

  if (from_factors && !sample_name %in% names(data)) {
    message("creating sampleName from factor columns")
    data <- data |>
      tidyr::unite(
        !!sym(sample_name),
        unique(unlist(configuration$factors)),
        remove = TRUE,
        sep = configuration$sep
      ) |>
      dplyr::select(sample_name, configuration$file_name) |>
      dplyr::distinct() |>
      dplyr::mutate(across(all_of(sample_name), function(x) {
        make.unique(x, sep = configuration$sep)
      })) |>
      dplyr::inner_join(data, by = configuration$file_name)
  } else if (!sample_name %in% names(data)) {
    message("creating sampleName from file_name column")
    data[[configuration$sample_name]] <- tools::file_path_sans_ext(basename(data[[configuration$file_name]]))
  } else {
    message("column sampleName already exists, using :", sample_name)
  }

  data <- data |>
    dplyr::select(-dplyr::all_of(dplyr::setdiff(unlist(configuration$factors), configuration$factor_keys())))

  # Make implicit NA's explicit
  if (!(configuration$isotope_label %in% colnames(data))) {
    warning(
      "no isotopeLabel column specified in the data, adding column isotopeLabel automatically and setting to 'light'."
    )
    data[[configuration$isotope_label]] <- "light"
  }
  if (!(configuration$ident_q_value %in% colnames(data))) {
    warning("no qValue column specified in the data. Creating column qValue and setting qValues to 0.")
    data[[configuration$ident_q_value]] <- 0
  }

  # Make implicit NA's explicit
  if (!(configuration$nr_children %in% colnames(data))) {
    warning("no nr_children column specified in the data, adding column nr_children and setting to 1.")
    data[[configuration$nr_children]] <- 1
  }

  # TODO add better warning....
  data <- data |> dplyr::select(c(configuration$id_vars(), configuration$value_vars()))

  txd <- data |>
    group_by(!!!syms(c(configuration$file_name, configuration$hierarchy_keys(), configuration$isotope_label))) |>
    summarize(n = n())
  if (any(txd$n > 1)) {
    offending <- dplyr::filter(dplyr::ungroup(txd), .data[["n"]] > 1)
    str <- paste0(
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
      paste(utils::capture.output(print(utils::head(offending, 5))), collapse = "\n"),
      "\nRe-run setup_analysis(..., debug = TRUE) to return the full count table for inspection ",
      "(rows where n > 1)."
    )
    if (debug) {
      warning(str)
      return(txd)
    }
    stop(str, call. = FALSE)
  }
  if (cc) {
    data <- .complete_cases_impl(
      data,
      configuration$file_name,
      configuration$sample_name,
      configuration$factor_keys(),
      configuration$isotope_label,
      configuration$hierarchy_keys()
    )
  }
  message("completing cases done")
  message("setup done")
  return(data)
}

#' Separates hierarchy columns into starting columns
#'
#'
#' @export
#' @param data data.frame
#' @param config AnlalysisConfiguration
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
  fkeys <- c(file_name, sample_name, factor_keys)
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
  message("Sampling ", size, paste(hierarchy_keys_depth, collapse = ","))
  hkdf <- pdata |> select(all_of(hierarchy_keys_depth)) |> distinct() |> sample_n(size = size)
  sdata <- inner_join(hkdf, pdata)
  return(sdata)
}
