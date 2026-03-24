# AnalysisConfiguration ----
#'
#' Analysis Configuration
#' @description
#' Analysis Configuration — holds all table annotations, hierarchy
#' definitions, factor definitions, and analysis parameters in a single
#' flat object.
#'
#' @family configuration
#' @export
#' @examples
#' istar <- sim_lfq_data_peptide_config()
#' config <- istar$config
#' stopifnot("AnalysisConfiguration" %in% class(config))
#' stopifnot(length(config$hierarchy_keys()) > 0)
#' stopifnot(length(config$factor_keys()) > 0)
#' stopifnot(length(config$get_response()) == 1)
#'
AnalysisConfiguration <- R6::R6Class(
  "AnalysisConfiguration",
  public = list(
    #' @field sep separator to use when uniting columns is necessary
    sep = "~",

    # --- fields from AnalysisTableAnnotation ---
    #' @field fileName column name of column containing raw file names
    fileName = NULL,
    #' @field sampleName (will be generated from factors or fileName)
    sampleName = "sampleName",
    #' @field normValue optional column with normalization values (e.g., Creatinine)
    normValue = NULL,

    #' @field isotopeLabel which column contains the isotope label (e.g. heavy or light), or light only if LFQ.
    isotopeLabel = "isotopeLabel",
    #' @field ident_qValue column name with identification QValues (smaller better)
    ident_qValue = "qValue",
    #' @field ident_Score column with identification score (larger better)
    ident_Score = character(),

    #' @field opt_rt optional column with rt information
    opt_rt = character(),
    #' @field opt_mz optional column with mz information
    opt_mz = character(),
    #' @field nr_children optional column containing for instance the number of peptides
    nr_children = "nr_children",

    #' @field workIntensity column which contains the intensities
    workIntensity = NULL,
    #' @field is_response_transformed are the intensities transformed for constant variance
    is_response_transformed = FALSE,
    #' @field bin_resp column with encoded missing information
    bin_resp = character(),

    #' @field factors Names of columns containing factors (annotations)
    factors = list(),
    #' @field factorDepth number of relevant factors (used by plotting functions etc)
    factorDepth = 1,

    #' @field hierarchy list with columns describing the measurement hierarchy (i.e. protein peptide precursor fragment)
    hierarchy = list(),
    #' @field hierarchyDepth At which depth do you want to model i.e. protein than 1
    hierarchyDepth = 1,

    # --- from AnalysisParameters ---
    #' @field min_peptides_protein minimum number of peptides per protein
    min_peptides_protein = 2,

    #' @description
    #' create AnalysisConfiguration
    #' @param analysisTableAnnotation optional, for backwards compatibility with old constructor
    #' @param analysisParameter optional, for backwards compatibility with old constructor
    initialize = function(analysisTableAnnotation = NULL, analysisParameter = NULL) {
      if (!is.null(analysisTableAnnotation)) {
        # Copy all fields from an AnalysisTableAnnotation or AnalysisConfiguration
        fields <- c(
          "fileName",
          "sampleName",
          "normValue",
          "isotopeLabel",
          "ident_qValue",
          "ident_Score",
          "opt_rt",
          "opt_mz",
          "nr_children",
          "workIntensity",
          "is_response_transformed",
          "bin_resp",
          "factors",
          "factorDepth",
          "hierarchy",
          "hierarchyDepth"
        )
        for (f in fields) {
          val <- analysisTableAnnotation[[f]]
          if (!is.null(val)) self[[f]] <- val
        }
        # Also copy min_peptides_protein if present on the source
        if (!is.null(analysisTableAnnotation[["min_peptides_protein"]])) {
          self$min_peptides_protein <- analysisTableAnnotation$min_peptides_protein
        }
      }
      if (!is.null(analysisParameter)) {
        if (!is.null(analysisParameter[["min_peptides_protein"]])) {
          self$min_peptides_protein <- analysisParameter$min_peptides_protein
        }
      }
    },

    # --- methods from AnalysisTableAnnotation ---
    #' @description
    #' Add name of intensity column
    #' @param colName name of intensity column
    set_response = function(colName) {
      self$workIntensity <- c(self$workIntensity, colName)
    },
    #' @description
    #' Get name of working intensity column
    get_response = function() {
      return(tail(self$workIntensity, n = 1))
    },
    #' @description
    #' Remove last name in array of working intensity column names
    pop_response = function() {
      res <- self$workIntensity[length(self$workIntensity)]
      self$workIntensity <- self$workIntensity[-length(self$workIntensity)]
      return(res)
    },

    #' @description
    #' Get factor keys
    #' @return array with keys
    factor_keys = function() {
      return(names(self$factors))
    },
    #' @description
    #' Get factor keys till factorDepth
    factor_keys_depth = function() {
      res <- head(self$factors, n = self$factorDepth)
      return(names(res))
    },

    #' @description
    #' get hierarchy keys
    #' @param rev return in reverse order
    #' @return array of column names
    hierarchy_keys = function(rev = FALSE) {
      if (rev) {
        return(rev(names(self$hierarchy)))
      } else {
        return(names(self$hierarchy))
      }
    },
    #' @description
    #' get hierarchy keys (deprecated alias for hierarchy_keys)
    #' @param rev return in reverse order
    #' @return array of column names
    hierarchyKeys = function(rev = FALSE) {
      self$hierarchy_keys(rev = rev)
    },

    #' @description
    #' get hierarchy keys up to depth
    #' @param names if TRUE names only if FALSE key value pairs
    #' @return array of column names
    hierarchy_keys_depth = function(names = TRUE) {
      res <- head(self$hierarchy, n = self$hierarchyDepth)
      res <- if (names) {
        names(res)
      } else {
        res
      }
      return(res)
    },
    #' @description
    #' get hierarchy keys up to depth (deprecated alias for hierarchy_keys_depth)
    #' @param names if TRUE names only if FALSE key value pairs
    #' @return array of column names
    hkeysDepth = function(names = TRUE) {
      self$hierarchy_keys_depth(names = names)
    },

    #' @description
    #' Id Columns which must be in the input data frame
    #' @return character array
    id_required = function() {
      id_vars <- c(
        self$fileName,
        unlist(self$factors),
        unlist(self$hierarchy),
        self$isotopeLabel
      )
      return(id_vars)
    },
    #' @description
    #' get names of columns annotating values (e.g. intensities)
    #' @return character array
    id_vars = function() {
      "Id Columns which must be in the output data frame"
      id_vars <- c(
        self$fileName,
        names(self$factors),
        names(self$hierarchy),
        self$isotopeLabel,
        self$sampleName,
        self$normValue
      )
      return(id_vars)
    },
    #' @description
    #' get names of columns containing observations e.g. (intensity, qValue, mz or rt)
    value_vars = function() {
      "Columns containing values"
      valueVars <- c(
        self$get_response(),
        self$ident_qValue,
        self$ident_Score,
        self$opt_mz,
        self$opt_rt,
        self$nr_children
      )
      return(valueVars)
    },
    #' @description
    #' get names of columns with sample annotations
    #'
    annotation_vars = function() {
      annotationVars <- c(self$fileName, self$sampleName, self$factor_keys(), self$normValue)
      return(annotationVars)
    }
  ),
  active = list(
    #' @field table deprecated, use config directly. Returns self for backwards compatibility.
    table = function(value) {
      if (!missing(value)) {
        return(invisible(self))
      }
      message("config$table is deprecated, use config directly")
      self
    },
    #' @field parameter deprecated, use config directly. Returns self for backwards compatibility.
    parameter = function(value) {
      if (!missing(value)) {
        return(invisible(self))
      }
      message("config$parameter is deprecated, use config directly")
      self
    }
  )
)


#' read minimal yaml to reconstruct configuration
#' @param dd list with table and parameter elements as produced by R6_extract_values
#' @export
#' @examples
#'
#' DEAconfig <- AnalysisConfiguration$new()
#' DEAconfig$fileName <- "Replicate.Name"
#' DEAconfig$hierarchy[["protein_Id"]] <- "Protein.Name"
#' DEAconfig$hierarchy[["peptide_Id"]] <- "Peptide.Sequence"
#' DEAconfig$hierarchy[["precursor_Id"]] <- c("Peptide.Sequence", "Precursor.Charge")
#' DEAconfig$hierarchy[["fragment_Id"]] <- c(
#'   "Peptide.Sequence", "Precursor.Charge",
#'   "Fragment.Ion", "Product.Charge"
#' )
#' DEAconfig$ident_qValue <- "annotation_QValue"
#' DEAconfig$set_response("Area")
#' DEAconfig$isotopeLabel <- "Isotope.Label"
#' configList <- prolfqua::R6_extract_values(DEAconfig)
#' stopifnot(class(configList) == "list")
#' config <- list_to_AnalysisConfiguration(configList)
#' all.equal(prolfqua::R6_extract_values(config), configList)
list_to_AnalysisConfiguration <- function(dd) {
  config <- AnalysisConfiguration$new()
  if (!is.null(dd$table)) {
    # Old nested format (backwards compat)
    for (i in names(dd$table)) {
      config[[i]] <- dd$table[[i]]
    }
    if (!is.null(dd$parameter)) {
      for (i in names(dd$parameter)) {
        if (i %in% names(config)) {
          config[[i]] <- dd$parameter[[i]]
        }
      }
    }
  } else {
    # New flat format
    for (i in names(dd)) {
      if (i != "sep") {
        config[[i]] <- dd[[i]]
      }
    }
  }
  if (!is.null(dd$sep)) {
    config$sep <- dd$sep
  }
  return(config)
}


#' Make reduced hierarchy configuration
#' @export
#' @keywords internal
#' @param config AnalysisConfiguration
#' @param workIntensity work intensity column
#' @param hierarchy new reduced hierarchy
#' @family configuration
#' @return AnalysisConfiguration with reduced hieararchy
#' @examples
#' bb <- prolfqua::sim_lfq_data_peptide_config()
#'
#' red <- make_reduced_hierarchy_config(bb$config,
#'  "testintensity",
#'  bb$config$hierarchy[1])
#' stopifnot(red$get_response() == "testintensity")
#' stopifnot(length(red$hierarchy) == 1)
make_reduced_hierarchy_config <- function(config, workIntensity, hierarchy) {
  newConfig <- config$clone(deep = TRUE)
  newConfig$hierarchy <- hierarchy
  newConfig$workIntensity <- workIntensity
  return(newConfig)
}

#' create interaction column from factors
#' @export
#' @keywords internal
#' @family configuration
#' @examples
#' xx <- data.frame(A = c("a","a","a"), B = c("d","d","e"))
#' x <- make_interaction_column(xx, c("B","A"))
#' x <- make_interaction_column(xx, c("A"))
#' bb <- prolfqua::sim_lfq_data_protein_config()
#' config <- bb$config
#' analysis <- bb$data
#'
#' config$factorDepth <- 1
#' make_interaction_column(analysis,
#'    config$factor_keys_depth())
#'
make_interaction_column <- function(data, columns, sep = ".") {
  intr <- dplyr::select(data, dplyr::all_of(columns))
  intr <- purrr::map_dfc(intr, factor)
  names(columns) <- columns
  newlev <- purrr::map2(columns, intr, function(x, y) {
    paste0(x, levels(y))
  })
  intr <- purrr::map2_dfc(columns, intr, paste0)
  intr <- purrr::map2_dfc(intr, newlev, forcats::fct_relevel)

  colnames(intr) <- paste0("interaction_", columns)
  colname <- "interaction"
  data <- data |> dplyr::mutate(!!colname := interaction(intr, sep = sep))
  return(data)
}


# Functions - Configuration ----
#' Extract all value slots in an R6 object
#' @param r6class r6 class
#' @family configuration
#' @export
R6_extract_values <- function(r6class) {
  tmp <- sapply(r6class, class)
  slots <- tmp[!tmp %in% c("environment", "function")]
  res <- list()
  for (i in names(slots)) {
    val <- r6class[[i]]
    if ("R6" %in% class(val)) {
      if (identical(val, r6class)) {
        next
      }
      res[[i]] <- R6_extract_values(val)
    } else {
      res[[i]] <- val
    }
  }
  return(res)
}


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

# Functions - summary ----

# Functions - summarize factors ----

#' Table of distinct factors (sample annotation)
#' @param pdata data.frame
#' @param configuration AnalysisConfiguration
#'
#' @export
#' @keywords internal
#' @family configuration
#' @examples
#'
#'
#' istar <- sim_lfq_data_peptide_config()
#'
#'
#' xx <- table_factors(istar$data,istar$config )
#' xx
#' xt <- xx |> dplyr::group_by(!!!rlang::syms(istar$config$factor_keys())) |>
#'  dplyr::summarize(n = dplyr::n())
#' xt
#' stopifnot(all(xt$n == 4))
#'
table_factors <- function(pdata, configuration) {
  factorsTab <- pdata |>
    dplyr::select(c(configuration$fileName, configuration$sampleName, configuration$factor_keys())) |>
    dplyr::distinct() |>
    arrange(!!sym(configuration$sampleName))
  return(factorsTab)
}

#' Table of distinct factors (sample annotation)
#' @param pdata data.frame
#' @param configuration AnalysisConfiguration
#'
#' @export
#' @keywords internal
#' @family configuration
#' @examples
#'
#'
#' istar <- sim_lfq_data_peptide_config()
#'
#'
#' xx <- table_factors_size(istar$data,istar$config )
#' stopifnot(all(xx$n == 4))
#'
table_factors_size <- function(pdata, configuration) {
  xx <- table_factors(pdata, configuration)
  xx <- xx |> dplyr::group_by(dplyr::across(configuration$factor_keys_depth())) |> dplyr::summarize(n = dplyr::n())
  return(xx)
}

# Functions - summarize hierarchies

#' Count distinct elements for each level of hierarchy and istope
#'
#' E.g. number of proteins, peptides, precursors in the dataset
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @export
#' @keywords internal
#' @family summary
#' @examples
#'
#' bb <- prolfqua::sim_lfq_data_peptide_config()
#'
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' x <- hierarchy_counts(data, config)
#' x$protein_Id
#' stopifnot(ncol(x) == length(config$hierarchy_keys()) + 1)
#' # select non existing protein
#' data0 <- data |> dplyr::filter( protein_Id == "XYZ")
#' tmp <- hierarchy_counts(data0, config)
#' stopifnot(nrow(tmp) == 0)
hierarchy_counts <- function(pdata, config) {
  hierarchy <- config$hierarchy_keys()
  res <- pdata |>
    dplyr::group_by(across(all_of(config$isotopeLabel))) |>
    dplyr::summarise(across(all_of(hierarchy), n_distinct))

  return(res)
}

#' Count distinct elements for each level of hierarchy per sample
#'
#' @export
#' @param pdata data.frame
#' @param configuration \code{\link{AnalysisConfiguration}}
#' @param nr_children integer, minimum number of children required
#'
#' @keywords internal
#' @family summary
#' @examples
#'
#' bb <- prolfqua::sim_lfq_data_peptide_config()
#'
#' config <- bb$config
#' data <- bb$data
#' res <- hierarchy_counts_sample(data, config, nr_children = 1)
#' x <- res("long")
#' # filters on peptide level
#' res <- hierarchy_counts_sample(data, config, nr_children = 2)
#' x2 <- res("long")
#' # filters on protein level based on peptide count
#' bb <- prolfqua::sim_lfq_data_protein_config()
#' res <- hierarchy_counts_sample(bb$data, bb$config, nr_children = 2)
#' x1 <- res()
#' res <- hierarchy_counts_sample(bb$data, bb$config, nr_children = 1)
#' x2 <- res()
#' x1$nr_children <- 2
#' x2$nr_children <- 1
#' xl <- dplyr::bind_rows(x1, x2)
#'
#' xl$nr_children |> table()
#' nudgeval <-  -mean(xl$protein_Id) * 0.05
#' ggplot2::ggplot(xl,
#'   ggplot2::aes(x = sampleName, y = protein_Id, fill = as.character(nr_children))) +
#'  ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge())
#'
hierarchy_counts_sample <- function(
  pdata,
  configuration,
  nr_children = 1
) {
  hierarchy <- configuration$hierarchy_keys()
  summary <- pdata |>
    dplyr::filter(
      !is.na(!!rlang::sym(configuration$get_response())),
      !!rlang::sym(configuration$nr_children) >= .env$nr_children
    ) |>
    dplyr::group_by(across(all_of(c(configuration$isotopeLabel, configuration$sampleName)))) |>
    dplyr::summarise(across(all_of(hierarchy), dplyr::n_distinct))

  res <- function(value = c("wide", "long", "plot")) {
    value <- match.arg(value)
    if (value == "wide") {
      return(summary)
    } else {
      long <- summary |>
        tidyr::pivot_longer(
          cols = -dplyr::all_of(c(configuration$isotopeLabel, configuration$sampleName)),
          names_to = "key",
          values_to = "nr"
        )
      if (value == "long") {
        return(long)
      } else if (value == "plot" && nrow(long) > 0) {
        nudgeval <- -mean(long$nr) * 0.05
        # TODO(WEW) check potential problem with sampleName
        ggplot2::ggplot(long, ggplot2::aes(x = !!rlang::sym(configuration$sampleName), y = .data$nr)) +
          ggplot2::geom_bar(stat = "identity", position = "dodge", colour = "black", fill = "white") +
          ggplot2::facet_wrap(~key, scales = "free_y", ncol = 1) +
          ggplot2::geom_text(ggplot2::aes(label = .data$nr), nudge_y = nudgeval, angle = 65) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
      }
    }
  }
  return(res)
}


#' Summarize hierarchy counts
#'
#' E.g compute number of peptides for each protein
#'
#' @export
#' @keywords internal
#' @family summary
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param hierarchy for which hierarchy level (default up to hierarchy depth)
#' @param factors which factors to include
#'
#' @examples
#'
#'
#'
#' bb <- sim_lfq_data_peptide_config()
#' data <- bb$data
#' configur <- bb$config
#' summarize_hierarchy(data, configur)
#' summarize_hierarchy(data, configur, factors = character())
#'
#' summarize_hierarchy(data, configur,
#'  hierarchy = configur$hierarchy_keys_depth() )
#' summarize_hierarchy(data, configur,
#'  hierarchy = NULL, factors = configur$factor_keys_depth() )
#' configur$hierarchyDepth = 1
#' summarize_hierarchy(data, configur,
#'  factors = configur$factor_keys_depth())
#' configur$hierarchyDepth = 2
#' summarize_hierarchy(data, configur)
#' configur$hierarchyDepth = 3
#' summarize_hierarchy(data, configur )
#' configur$hierarchyDepth = 4
#' summarize_hierarchy(data, configur )
#'
summarize_hierarchy <- function(pdata, config, hierarchy = config$hierarchy_keys_depth(), factors = character()) {
  all_hierarchy <- c(config$isotopeLabel, config$hierarchy_keys())

  precursor <- pdata |> dplyr::select(dplyr::all_of(c(factors, all_hierarchy))) |> dplyr::distinct()
  x3 <- precursor |>
    dplyr::group_by(across(all_of(c(factors, hierarchy)))) |>
    dplyr::summarize(across(all_of(base::setdiff(all_hierarchy, hierarchy)), list(n = dplyr::n_distinct)))
  return(x3)
}

# Functions - Handling isotopes ----

# Computing protein Intensity summaries ---
