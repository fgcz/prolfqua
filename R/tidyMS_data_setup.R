#' Setup a tidy table compatible with a \code{\link{AnalysisConfiguration}}
#'
#' Extracts columns relevant for a configuration from a data frame
#' and create new columns e.g. sampleName column etc.
#' @param data data.frame
#' @param configuration AnalysisConfiguration
#' @param cc complete cases default TRUE
#' @param from_factors if TRUE, create sampleName from factor columns
#' @export
#' @family configuration
#' @examples
#'
#' skylineconfig <- AnalysisConfiguration$new()
#' skylineconfig$fileName <- "Replicate.Name"
#' skylineconfig$hierarchy[["protein_Id"]] <- "Protein.Name"
#' skylineconfig$hierarchy[["peptide_Id"]] <- "Peptide.Sequence"
#' skylineconfig$hierarchy[["precursor_Id"]] <- c("Peptide.Sequence", "Precursor.Charge")
#' skylineconfig$hierarchy[["fragment_Id"]] <- c(
#'   "Peptide.Sequence", "Precursor.Charge",
#'   "Fragment.Ion", "Product.Charge"
#' )
#' skylineconfig$ident_qValue <- "Detection.Q.Value"
#' skylineconfig$set_response("Area")
#' skylineconfig$isotopeLabel <- "Isotope.Label.Type"
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
#' config$fileName = "sample"
#' config$factors["group_"] = "group"
#' config$hierarchy[["protein_Id"]] = c("proteinID", "idtype2")
#' config$hierarchy[["peptide_Id"]] = "peptideID"
#' config$set_response("abundance")
#' config$normValue = "creatinine"
#'
#' adata <- setup_analysis(data, config)
#'
setup_analysis <- function(data, configuration, cc = TRUE, from_factors = FALSE) {
  configuration <- configuration$clone(deep = TRUE)
  if (is.null(configuration$fileName)) {
    stop("fileName column is not specified in configuration.")
  }
  if (!configuration$fileName %in% colnames(data)) {
    stop("File name column :", configuration$fileName, ", is missing in data.")
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

  sampleName <- configuration$sampleName

  if (from_factors && !sampleName %in% names(data)) {
    message("creating sampleName from factor columns")
    data <- data |>
      tidyr::unite(
        !!sym(sampleName),
        unique(unlist(configuration$factors)),
        remove = TRUE,
        sep = configuration$sep
      ) |>
      dplyr::select(sampleName, configuration$fileName) |>
      dplyr::distinct() |>
      dplyr::mutate(across(all_of(sampleName), function(x) {
        make.unique(x, sep = configuration$sep)
      })) |>
      dplyr::inner_join(data, by = configuration$fileName)
  } else if (!sampleName %in% names(data)) {
    message("creating sampleName from fileName column")
    data[[configuration$sampleName]] <- tools::file_path_sans_ext(basename(data[[configuration$fileName]]))
  } else {
    message("column sampleName already exists, using :", sampleName)
  }

  data <- data |>
    dplyr::select(-dplyr::all_of(dplyr::setdiff(unlist(configuration$factors), configuration$factor_keys())))

  # Make implicit NA's explicit
  if (!(configuration$isotopeLabel %in% colnames(data))) {
    warning(
      "no isotopeLabel column specified in the data, adding column isotopeLabel automatically and setting to 'light'."
    )
    data[[configuration$isotopeLabel]] <- "light"
  }
  if (!(configuration$ident_qValue %in% colnames(data))) {
    warning("no qValue column specified in the data. Creating column qValue and setting qValues to 0.")
    data[[configuration$ident_qValue]] <- 0
  }

  # Make implicit NA's explicit
  if (!(configuration$nr_children %in% colnames(data))) {
    warning("no nr_children column specified in the data, adding column nr_children and setting to 1.")
    data[[configuration$nr_children]] <- 1
  }

  # TODO add better warning....
  data <- data |> dplyr::select(c(configuration$id_vars(), configuration$value_vars()))

  txd <- data |>
    group_by(!!!syms(c(configuration$fileName, configuration$hierarchy_keys(), configuration$isotopeLabel))) |>
    summarize(n = n())
  if (any(txd$n > 1)) {
    str <- paste(
      "There is more than ONE observations for each : ",
      paste(configuration$hierarchy_keys(), collapse = ", "),
      ",\n",
      "and sample : ",
      configuration$sampleName,
      "; (filename) : ",
      configuration$fileName,
      "\n"
    )
    warning(str)
    warning("Please inspect the returned dataframe. Check for rows where n > 1\n e.g. res |> dplyr::filter(n > 1)")
    return(txd)
  }
  if (cc) {
    data <- complete_cases(data, configuration)
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

#' Complete cases
#'
#' The tidy table does not need to contain missing data.
#' This function re-establishes the missing observations in a sample.
#'
#' @export
#' @param pdata data.frame
#' @param config AnlalysisConfiguration
#' @keywords internal
#' @family configuration
#' @examples
#'
#' bb <- sim_lfq_data_protein_config()
#'
#' xx <- complete_cases(bb$data, bb$config)
#' stopifnot(nrow(bb$data) <= nrow(xx))
#'
complete_cases <- function(pdata, config) {
  message("completing cases")
  fkeys <- c(config$fileName, config$sampleName, config$factor_keys())
  hkeys <- c(config$isotopeLabel, config$hierarchy_keys())
  res <- tidyr::complete(
    pdata,
    tidyr::nesting(!!!syms(fkeys)),
    tidyr::nesting(!!!syms(hkeys))
  )
  return(res)
}


#' Sample subset of proteins/peptides/precursors
#' @param size size of sample
#' @param pdata tidy table
#' @param config \code{\link{AnalysisConfiguration}}
#' @export
#' @keywords internal
#' @family configuration
#'
sample_subset <- function(size, pdata, config) {
  hk <- config$hierarchy_keys_depth()
  message("Sampling ", size, paste(hk, collapse = ","))
  hkdf <- pdata |> select(all_of(hk)) |> distinct() |> sample_n(size = size)
  sdata <- inner_join(hkdf, pdata)
  return(sdata)
}
