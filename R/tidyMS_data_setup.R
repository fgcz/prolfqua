#' Require a configured and present sample-file column
#' @param data long-format measurement data frame
#' @param configuration AnalysisConfiguration defining the column roles
#' @noRd
.validate_setup_file_name <- function(data, configuration) {
  if (is.null(configuration$file_name)) {
    abort_invalid_config("`file_name` column is not specified in the AnalysisConfiguration.")
  }
  if (!configuration$file_name %in% colnames(data)) {
    abort_missing_columns(configuration$file_name)
  }
}

#' Combine configured source columns into hierarchy identifiers
#' @param data long-format measurement data frame
#' @param configuration AnalysisConfiguration defining the column roles
#' @noRd
.setup_hierarchy_columns <- function(data, configuration) {
  for (i in seq_along(configuration$hierarchy)) {
    .inform_united(configuration$hierarchy_keys()[i], configuration$hierarchy[[i]], configuration$sep)
    data <- tidyr::unite(
      data,
      !!sym(configuration$hierarchy_keys()[i]),
      dplyr::all_of(configuration$hierarchy[[i]]),
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

#' Construct categorical factor columns from their configured sources
#' @param data long-format measurement data frame
#' @param configuration AnalysisConfiguration defining the column roles
#' @noRd
#' @srrstats {G2.4c} Configured factor columns are explicitly converted with as.character.
#' @srrstats {G2.4e} Factor labels are explicitly converted with as.character instead of using integer codes.
#' @srrstats {G2.5} Setup documents and enforces character labels for both ordered and unordered predictor factors.
.setup_factor_columns <- function(data, configuration) {
  if (length(configuration$factors) == 0) {
    abort_invalid_config(c(
      "The AnalysisConfiguration specifies no factors (explanatory variables).",
      i = paste0(
        "Set them with `config$factors[\"Condition\"] <- \"columnName\"`, where `Condition` names the variable ",
        "and `columnName` is the data column holding it."
      )
    ))
  }

  for (i in seq_along(configuration$factors)) {
    if (length(configuration$factors[[i]]) > 1) {
      .inform_united(configuration$factor_keys()[i], configuration$factors[[i]], configuration$sep)
      data <- tidyr::unite(
        data,
        !!sym(configuration$factor_keys()[i]),
        dplyr::all_of(configuration$factors[[i]]),
        remove = FALSE,
        sep = configuration$sep
      )
    } else {
      newname <- configuration$factor_keys()[i]
      source <- data[[configuration$factors[[i]]]]
      if (!is.character(source)) {
        inform_prolfqua(
          "Factor `",
          newname,
          "`: ",
          class(source)[1],
          " column `",
          configuration$factors[[i]],
          "` converted to character labels."
        )
      }
      data <- dplyr::mutate(
        data,
        !!newname := as.character(!!sym(configuration$factors[[i]]))
      )
    }
  }
  data
}

#' Announce that several source columns are united into one key column
#' @param key resulting key column
#' @param sources configured source columns
#' @param sep separator
#' @noRd
.inform_united <- function(key, sources, sep) {
  if (length(sources) > 1L) {
    inform_prolfqua("`", key, "`: columns ", toString(paste0("`", sources, "`")), " united with separator '", sep, "'.")
  }
}

#' Create sample names from file names when absent
#' @param data long-format measurement data frame
#' @param configuration AnalysisConfiguration defining the column roles
#' @noRd
.setup_sample_name <- function(data, configuration) {
  sample_name <- configuration$sample_name
  if (!sample_name %in% names(data)) {
    inform_prolfqua("creating sampleName from file_name column")
    data[[sample_name]] <- tools::file_path_sans_ext(basename(data[[configuration$file_name]]))
    return(data)
  }

  data
}

#' Insert a missing configured column and report the default
#' @param data long-format measurement data frame
#' @param column name of the column to insert
#' @param value default value assigned to the new column
#' @param text diagnostic message passed to notify
#' @param notify function reporting an inserted column
#' @noRd
.add_setup_default_column <- function(data, column, value, text, notify) {
  if (column %in% colnames(data)) {
    return(data)
  }
  notify(text)
  data[[column]] <- value
  data
}

#' isotopeLabel and qValue are optional and only reported. A missing nr_children
#'
#' warns: the default of 1 changes the fitting weights of aggregated data.
#' @param data long-format measurement data frame
#' @param configuration AnalysisConfiguration defining the column roles
#' @noRd
.add_setup_default_columns <- function(data, configuration) {
  data <- .add_setup_default_column(
    data,
    configuration$isotope_label,
    "light",
    "no isotopeLabel column specified in the data, adding column isotopeLabel automatically and setting to 'light'.",
    inform_prolfqua
  )
  data <- .add_setup_default_column(
    data,
    configuration$ident_q_value,
    0,
    "no qValue column specified in the data. Creating column qValue and setting qValues to 0.",
    inform_prolfqua
  )
  .add_setup_default_column(
    data,
    configuration$nr_children,
    1,
    "no nr_children column specified in the data, adding column nr_children and setting to 1.",
    warning
  )
}

#' Count rows per hierarchy identifier and sample
#' @param data long-format measurement data frame
#' @param configuration AnalysisConfiguration defining the column roles
#' @noRd
.setup_duplicate_counts <- function(data, configuration) {
  data |>
    group_by(
      !!!syms(c(
        configuration$file_name,
        configuration$hierarchy_keys(),
        configuration$isotope_label
      ))
    ) |>
    summarize(n = n(), .groups = "drop")
}

#' Describe duplicate measurements and their offending keys
#' @param offending duplicate measurement keys
#' @param configuration AnalysisConfiguration defining the column roles
#' @noRd
.setup_duplicate_message <- function(offending, configuration) {
  paste0(
    "setup_analysis: there is more than one observation for ",
    nrow(offending),
    " ",
    toString(configuration$hierarchy_keys()),
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

#' Reject duplicate measurements or return their diagnostic counts
#' @param counts counts of measurements per subject and sample
#' @param configuration AnalysisConfiguration defining the column roles
#' @param debug return diagnostic duplicate counts instead of aborting
#' @noRd
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
#' @details `setup_analysis()` is optional: `LFQData$new(setup = FALSE)`, the
#'   default, applies none of the steps below. Each step that changes the data
#'   says so with a message (class `prolfqua_message`) or warning:
#'   \itemize{
#'   \item several configured hierarchy or factor source columns are united into
#'     one key column (separator `configuration$sep`); avoid by configuring one
#'     column per key;
#'   \item factor source columns that are not character (numeric, factor,
#'     ordered) are converted to character labels; ordered levels and their
#'     ordering are not retained; avoid by supplying character columns;
#'     names of \code{configuration$factors} become the predictor column names;
#'   \item a missing sample name is derived from the file name; avoid by
#'     supplying the sample name column;
#'   \item missing isotope labels, identification q-values and child counts
#'     receive default values (child counts of one with a warning); avoid by
#'     supplying these columns;
#'   \item unconfigured columns and their attributes are dropped; the returned
#'     table keeps only configured analysis columns;
#'   \item feature and sample combinations absent from the input are added as
#'     rows with NA response, so every feature has a row in every sample.
#'   }
#' @param data a data.frame with the configured source hierarchy, factor,
#'   file-name and numeric response columns
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
#' data <- sim_lfq_data(n_proteins = 10, peptide = TRUE, n_per_group = 4)
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
#' @srrstatsTODO {G2.4} Setup explicitly converts categorical predictors; specialized input conversion is not fully
#'   audited.
#' @srrstatsTODO {G2.4d} Design matrices factor categorical labels as documented; limit: setup discards supplied level
#'   order.
#' @srrstats {G2.9} Setup announces conversion of factor columns to character, united key columns, inserted default
#'   columns and dropped unconfigured columns; `LFQData$new()` announces dropped response attributes.
#' @srrstats {RE1.2} Setup help documents categorical predictors, numeric coercion and loss of ordered-factor
#'   semantics.
#' @srrstats {RE1.3a} Setup help states retained columns, discarded attributes, categorical conversion and inserted
#'   defaults.
#' @srrstats {RE2.0} Every setup transformation is documented in Details, announced by a message or warning when
#'   it changes the data, and avoidable: setup is opt-in (`LFQData$new(setup = FALSE)` is the default), and each step
#'   is skipped when the input already has the configured form (test-setup-messages.R).
setup_analysis <- function(data, configuration, debug = FALSE) {
  .assert_response_data(data, configuration)
  .validate_setup_file_name(data, configuration)
  required <- unique(c(unlist(configuration$hierarchy), unlist(configuration$factors)))
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    abort_missing_columns(missing)
  }
  configuration <- configuration$clone(deep = TRUE)
  data <- .setup_hierarchy_columns(data, configuration)
  data <- .setup_factor_columns(data, configuration)
  data <- .setup_sample_name(data, configuration)
  data <- data |>
    dplyr::select(-dplyr::all_of(dplyr::setdiff(unlist(configuration$factors), configuration$factor_keys())))
  data <- .add_setup_default_columns(data, configuration)
  dropped <- setdiff(names(data), c(configuration$id_vars(), configuration$value_vars()))
  if (length(dropped) > 0L) {
    inform_prolfqua("setup_analysis: dropped ", length(dropped), " unconfigured columns: ", toString(dropped), ".")
  }
  data <- data |> dplyr::select(dplyr::all_of(c(configuration$id_vars(), configuration$value_vars())))

  duplicate_counts <- .setup_duplicate_counts(data, configuration)
  duplicate_result <- .handle_setup_duplicates(
    duplicate_counts,
    configuration,
    debug
  )
  if (!is.null(duplicate_result)) {
    return(duplicate_result)
  }
  n_input <- nrow(data)
  data <- .complete_cases_impl(
    data,
    configuration$file_name,
    configuration$sample_name,
    configuration$factor_keys(),
    configuration$isotope_label,
    configuration$hierarchy_keys()
  )
  if (nrow(data) > n_input) {
    inform_prolfqua(
      "setup_analysis: added ",
      nrow(data) - n_input,
      " rows with NA response for feature and sample combinations absent from the input."
    )
  }
  data
}

#' Separates hierarchy columns into starting columns
#'
#'
#' @export
#' @param data data.frame
#' @param config AnalysisConfiguration
#' @return `data` with each hierarchy key at the current depth split at `config$sep` into its component columns,
#'   which are added after the key column
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
      data <- data |>
        tidyr::separate_wider_delim(
          dplyr::all_of(hkey),
          delim = config$sep,
          names = config$hierarchy[[hkey]],
          too_few = "align_start",
          too_many = "drop",
          cols_remove = FALSE
        ) |>
        dplyr::relocate(dplyr::all_of(hkey), .before = dplyr::all_of(config$hierarchy[[hkey]][1]))
    }
  }
  data
}

#' Internal implementation for complete_cases (used by setup_analysis before LFQData exists)
#' @param pdata long-format measurement table
#' @param file_name column identifying the input measurement files
#' @param sample_name column identifying samples
#' @param factor_keys configured explanatory-variable columns
#' @param isotope_label column identifying isotope labels
#' @param hierarchy_keys hierarchy identifier columns
#' @noRd
.complete_cases_impl <- function(pdata, file_name, sample_name, factor_keys, isotope_label, hierarchy_keys) {
  fkeys <- unique(c(file_name, sample_name, factor_keys))
  hkeys <- c(isotope_label, hierarchy_keys)
  res <- tidyr::complete(
    pdata,
    tidyr::nesting(!!!syms(fkeys)),
    tidyr::nesting(!!!syms(hkeys))
  )
  res
}

#' Complete cases
#'
#' The tidy table does not need to contain missing data.
#' This function re-establishes the missing observations in a sample.
#'
#' @export
#' @param lfqdata [LFQData] object with the data and its configuration
#' @return a tibble of the long-format data with one row for every combination of sample and hierarchy element;
#'   added rows have `NA` intensities
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
#' @param hierarchy_keys_depth character vector - hierarchy columns at current depth
#' @return a tibble with all rows of `pdata` for `size` randomly drawn hierarchy elements
#' @export
#' @keywords internal
#' @family configuration
#'
#' @examples
#' bb <- sim_lfq_data_peptide_config(n_proteins = 5)
#' subset <- sample_subset(2, bb$data, bb$config$hierarchy_keys_depth())
#' length(unique(subset[[bb$config$hierarchy_keys_depth()[1]]]))
sample_subset <- function(size, pdata, hierarchy_keys_depth) {
  hkdf <- pdata |> select(all_of(hierarchy_keys_depth)) |> distinct() |> dplyr::slice_sample(n = size)
  sdata <- inner_join(hkdf, pdata, by = hierarchy_keys_depth, relationship = "one-to-many")
  sdata
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

#' Keep only the parents with enough children
#'
#' Keeps the rows of the `parent` entries (e.g. proteins) with at least
#' `threshold` distinct `child` entries (e.g. peptides), counted by
#' [nr_B_in_A()], and adds that count. By default it counts the second
#' hierarchy level within the first, independent of `hierarchy_depth`. The
#' package dataset `data_ionstar` calls it as `(pdata, config)` from its
#' `filtered()` method.
#' @param pdata data.frame in the long format of `config`
#' @param config AnalysisConfiguration
#' @param parent parent key column, default the first hierarchy key
#' @param child child key column, default the second hierarchy key
#' @param threshold minimum number of children, default `config$min_peptides_protein`
#' @return list with the filtered `data` and the `name` of the count column
#' @keywords internal
#' @family preprocessing
#' @examples
#' istar <- prolfqua::sim_lfq_data_peptide_config()
#' filtered <- prolfqua:::filter_proteins_by_peptide_count(
#'   istar$data, istar$config, parent = "protein_Id", child = "peptide_Id", threshold = 2
#' )
#' stopifnot(all(filtered$data[[filtered$name]] >= 2))
filter_proteins_by_peptide_count <- function(
  pdata,
  config,
  parent = config$hierarchy_keys()[1],
  child = config$hierarchy_keys()[2],
  threshold = config$min_peptides_protein
) {
  if (is.na(child)) {
    abort_bad_argument("child", "name a hierarchy level below `parent`", not = "the hierarchy has one level")
  }
  missing <- setdiff(c(parent, child), colnames(pdata))
  if (length(missing) > 0) {
    abort_missing_columns(missing, data_nm = "pdata")
  }
  counts <- nr_B_in_A(pdata, parent, child)
  name <- setdiff(names(counts), parent)
  data <- dplyr::inner_join(pdata, counts, by = parent, relationship = "many-to-one")
  list(data = dplyr::filter(data, .data[[name]] >= .env$threshold), name = name)
}
