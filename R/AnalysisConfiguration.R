# AnalysisConfiguration ----
#'
#' Analysis Configuration
#' @description
#' Analysis Configuration - holds all table annotations, hierarchy
#' definitions, factor definitions, and analysis parameters in a single
#' flat object.
#'
#' @family configuration
#' @return An R6 class generator.
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
    #' @field file_name column name of column containing raw file names
    file_name = NULL,
    #' @field sample_name (will be generated from factors or file_name)
    sample_name = "sampleName",
    #' @field norm_value optional column with normalization values (e.g., Creatinine)
    norm_value = NULL,

    #' @field identification_confidence column with the identification confidence of each feature in each sample,
    #'   in \[0, 1\], larger is better (e.g. 1 - PEP); used as model weight
    identification_confidence = "identification_confidence",
    #' @field identification_confidence_feature column with one identification confidence per feature, the maximum
    #'   of \code{identification_confidence} over samples
    identification_confidence_feature = "identification_confidence_feature",
    #' @field ident_score column with identification score (larger better)
    ident_score = character(),

    #' @field opt_rt optional column with rt information
    opt_rt = character(),
    #' @field opt_mz optional column with mz information
    opt_mz = character(),
    #' @field opt_se optional column with standard errors (e.g. from limpa aggregation)
    opt_se = character(),
    #' @field nr_children optional column containing for instance the number of peptides
    nr_children = "nr_children",
    #' @field nr_children_feature column with one child count per feature, the maximum of \code{nr_children} over
    #'   samples
    nr_children_feature = "nr_children_feature",

    #' @field work_intensity column which contains the intensities
    work_intensity = NULL,
    #' @field is_response_transformed are the intensities transformed for constant variance
    is_response_transformed = FALSE,
    #' @field bin_resp column with encoded missing information
    bin_resp = character(),

    #' @field pattern_decoys optional regex identifying decoy / reverse-database
    #'   ids, unioned with the built-in default prefixes (see
    #'   \code{\link{is_decoy}}). \code{NULL} (default) leaves decoy handling off;
    #'   set it (even to \code{""} for defaults only) to opt into flagging and the
    #'   targets-only model fit.
    pattern_decoys = NULL,
    #' @field pattern_contaminants optional regex identifying contaminant ids
    #'   (keratin, trypsin, BSA, ...). \code{NULL} (default) means contaminants are
    #'   not flagged; set it to enable contaminant flagging / removal.
    pattern_contaminants = NULL,

    #' @field factors Names of columns containing factors (annotations)
    factors = list(),
    #' @field factor_depth number of relevant factors (used by plotting functions etc)
    factor_depth = 1,

    #' @field hierarchy list with columns describing the measurement hierarchy (i.e. protein peptide precursor fragment)
    hierarchy = list(),
    #' @field hierarchy_depth At which depth do you want to model i.e. protein than 1
    hierarchy_depth = 1,

    # --- from AnalysisParameters ---
    #' @field min_peptides_protein minimum number of peptides per protein
    min_peptides_protein = 2,

    # --- methods from AnalysisTableAnnotation ---
    #' @description
    #' Add name of intensity column
    #' @param col_name name of intensity column
    set_response = function(col_name) {
      self$work_intensity <- c(self$work_intensity, col_name)
    },
    #' @description
    #' Get name of working intensity column
    get_response = function() {
      tail(self$work_intensity, n = 1)
    },
    #' @description
    #' Remove last name in array of working intensity column names
    pop_response = function() {
      res <- self$work_intensity[length(self$work_intensity)]
      self$work_intensity <- self$work_intensity[-length(self$work_intensity)]
      res
    },

    #' @description
    #' Get factor keys
    #' @return array with keys
    factor_keys = function() {
      names(self$factors)
    },
    #' @description
    #' Get the first `factor_depth` factor keys
    #' @return character vector of factor column names
    factor_keys_depth = function() {
      res <- head(self$factors, n = self$factor_depth)
      names(res)
    },

    #' @description
    #' get hierarchy keys
    #' @return array of column names
    hierarchy_keys = function() {
      names(self$hierarchy)
    },

    #' @description
    #' get hierarchy keys up to depth
    #' @param names if TRUE names only if FALSE key value pairs
    #' @return array of column names
    hierarchy_keys_depth = function(names = TRUE) {
      res <- head(self$hierarchy, n = self$hierarchy_depth)
      if (names) names(res) else res
    },

    #' @description
    #' get names of columns annotating values (e.g. intensities)
    #' @return character array
    id_vars = function() {
      id_vars <- c(
        self$file_name,
        names(self$factors),
        names(self$hierarchy),
        self$sample_name,
        self$norm_value
      )
      id_vars
    },
    #' @description
    #' get names of columns containing observations e.g. (intensity, identification confidence, mz or rt)
    value_vars = function() {
      value_vars <- c(
        self$get_response(),
        self$identification_confidence,
        self$identification_confidence_feature,
        self$ident_score,
        self$opt_mz,
        self$opt_rt,
        self$opt_se,
        self$nr_children,
        self$nr_children_feature
      )
      value_vars
    },
    #' @description
    #' get names of columns with sample annotations
    #'
    annotation_vars = function() {
      annotation_vars <- c(self$file_name, self$sample_name, self$factor_keys(), self$norm_value)
      annotation_vars
    }
  )
)


#' Reconstruct a configuration from minimal YAML
#'
#' The inverse of [r6_extract_values()]: rebuilds an [AnalysisConfiguration] from the list read from a YAML file,
#' in the current flat format or the older nested one. Unknown keys raise an error.
#' @param dd list with table and parameter elements as produced by r6_extract_values
#' @return an [AnalysisConfiguration] with the fields set from `dd`.
#' @export
#' @examples
#'
#' DEAconfig <- AnalysisConfiguration$new()
#' DEAconfig$file_name <- "Replicate.Name"
#' DEAconfig$hierarchy[["protein_Id"]] <- "Protein.Name"
#' DEAconfig$hierarchy[["peptide_Id"]] <- "Peptide.Sequence"
#' DEAconfig$hierarchy[["precursor_Id"]] <- c("Peptide.Sequence", "Precursor.Charge")
#' DEAconfig$hierarchy[["fragment_Id"]] <- c(
#'   "Peptide.Sequence", "Precursor.Charge",
#'   "Fragment.Ion", "Product.Charge"
#' )
#' DEAconfig$identification_confidence <- "annotation_confidence"
#' DEAconfig$set_response("Area")
#' configList <- prolfqua::r6_extract_values(DEAconfig)
#' stopifnot(class(configList) == "list")
#' config <- list_to_analysis_configuration(configList)
#' all.equal(prolfqua::r6_extract_values(config), configList)
list_to_analysis_configuration <- function(dd) {
  config <- AnalysisConfiguration$new()
  fields <- names(Filter(function(x) !is.function(x) && !is.environment(x), as.list(config)))
  if (!is.null(dd$table)) {
    .set_nested_configuration(config, dd, fields)
  } else {
    # New flat format
    dd$isotope_label <- NULL
    .assert_configuration_keys(setdiff(names(dd), "sep"), fields)
    for (i in names(dd)) {
      if (i != "sep") {
        config[[i]] <- dd[[i]]
      }
    }
  }
  if (!is.null(dd$sep)) {
    config$sep <- dd$sep
  }
  config
}

#' Set a configuration from the old nested format
#'
#' Keys of the `table` block must be configuration fields; keys of the `parameter` block that are not fields
#' (for example method names) are skipped, as before.
#' @param config AnalysisConfiguration to set
#' @param dd list with `table` and `parameter` elements
#' @param fields configuration fields
#' @param call calling environment named in errors
#' @return `config`, invisibly
#' @noRd
.set_nested_configuration <- function(config, dd, fields, call = rlang::caller_env()) {
  # configurations stored before prolfqua became label-free carry the removed `isotope_label`
  dd$table$isotope_label <- NULL
  .assert_configuration_keys(names(dd$table), fields, call = call)
  for (i in names(dd$table)) {
    config[[i]] <- dd$table[[i]]
  }
  for (i in intersect(names(dd$parameter), fields)) {
    config[[i]] <- dd$parameter[[i]]
  }
  invisible(config)
}

#' Abort unless every key names a field of AnalysisConfiguration
#' @param keys list names to set on the configuration
#' @param fields configuration fields
#' @param call calling environment named in the error
#' @noRd
.assert_configuration_keys <- function(keys, fields, call = rlang::caller_env()) {
  unknown <- setdiff(keys, fields)
  if (length(unknown) > 0) {
    abort_bad_argument(
      "dd",
      "name only AnalysisConfiguration fields",
      not = paste0("unknown keys ", toString(unknown)),
      call = call
    )
  }
}


#' Make reduced hierarchy configuration
#' @param config AnalysisConfiguration
#' @param work_intensity work intensity column
#' @param hierarchy new reduced hierarchy
#' @return AnalysisConfiguration with reduced hierarchy
#' @noRd
.make_reduced_hierarchy_config <- function(config, work_intensity, hierarchy) {
  new_config <- config$clone(deep = TRUE)
  new_config$hierarchy <- hierarchy
  new_config$work_intensity <- work_intensity
  new_config
}


# Functions - Configuration ----
#' Extract all value slots in an R6 object
#'
#' Turns an R6 object, such as an [AnalysisConfiguration], into a plain list that can be written to YAML;
#' [list_to_analysis_configuration()] reads it back.
#' @param r6class r6 class
#' @family configuration
#' @return a named list with the value fields of `r6class`; fields holding other R6 objects become nested
#'   lists, functions and environments are dropped.
#' @export
#' @examples
#' bb <- sim_lfq_data_peptide_config()
#' values <- r6_extract_values(bb$config)
#' names(values)
r6_extract_values <- function(r6class) {
  tmp <- vapply(r6class, function(x) class(x)[1], character(1))
  slots <- tmp[!tmp %in% c("environment", "function")]
  res <- list()
  for (i in names(slots)) {
    val <- r6class[[i]]
    if (inherits(val, "R6")) {
      if (identical(val, r6class)) {
        next
      }
      res[[i]] <- r6_extract_values(val)
    } else {
      res[[i]] <- val
    }
  }
  res
}
