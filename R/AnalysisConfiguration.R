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
    #' @field file_name column name of column containing raw file names
    file_name = NULL,
    #' @field sample_name (will be generated from factors or file_name)
    sample_name = "sampleName",
    #' @field norm_value optional column with normalization values (e.g., Creatinine)
    norm_value = NULL,

    #' @field isotope_label which column contains the isotope label (e.g. heavy or light), or light only if LFQ.
    isotope_label = "isotopeLabel",
    #' @field ident_q_value column name with identification QValues (smaller better)
    ident_q_value = "qValue",
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

    #' @field work_intensity column which contains the intensities
    work_intensity = NULL,
    #' @field is_response_transformed are the intensities transformed for constant variance
    is_response_transformed = FALSE,
    #' @field bin_resp column with encoded missing information
    bin_resp = character(),

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

    #' @description
    #' create AnalysisConfiguration
    initialize = function() {},

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
      return(tail(self$work_intensity, n = 1))
    },
    #' @description
    #' Remove last name in array of working intensity column names
    pop_response = function() {
      res <- self$work_intensity[length(self$work_intensity)]
      self$work_intensity <- self$work_intensity[-length(self$work_intensity)]
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
      res <- head(self$factors, n = self$factor_depth)
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
    #' get hierarchy keys up to depth
    #' @param names if TRUE names only if FALSE key value pairs
    #' @return array of column names
    hierarchy_keys_depth = function(names = TRUE) {
      res <- head(self$hierarchy, n = self$hierarchy_depth)
      res <- if (names) {
        names(res)
      } else {
        res
      }
      return(res)
    },

    #' @description
    #' Id Columns which must be in the input data frame
    #' @return character array
    id_required = function() {
      id_vars <- c(
        self$file_name,
        unlist(self$factors),
        unlist(self$hierarchy),
        self$isotope_label
      )
      return(id_vars)
    },
    #' @description
    #' get names of columns annotating values (e.g. intensities)
    #' @return character array
    id_vars = function() {
      "Id Columns which must be in the output data frame"
      id_vars <- c(
        self$file_name,
        names(self$factors),
        names(self$hierarchy),
        self$isotope_label,
        self$sample_name,
        self$norm_value
      )
      return(id_vars)
    },
    #' @description
    #' get names of columns containing observations e.g. (intensity, qValue, mz or rt)
    value_vars = function() {
      "Columns containing values"
      value_vars <- c(
        self$get_response(),
        self$ident_q_value,
        self$ident_score,
        self$opt_mz,
        self$opt_rt,
        self$opt_se,
        self$nr_children
      )
      return(value_vars)
    },
    #' @description
    #' get names of columns with sample annotations
    #'
    annotation_vars = function() {
      annotation_vars <- c(self$file_name, self$sample_name, self$factor_keys(), self$norm_value)
      return(annotation_vars)
    }
  )
)


#' read minimal yaml to reconstruct configuration
#' @param dd list with table and parameter elements as produced by R6_extract_values
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
#' DEAconfig$ident_q_value <- "annotation_QValue"
#' DEAconfig$set_response("Area")
#' DEAconfig$isotope_label <- "Isotope.Label"
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
#' @param work_intensity work intensity column
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
make_reduced_hierarchy_config <- function(config, work_intensity, hierarchy) {
  new_config <- config$clone(deep = TRUE)
  new_config$hierarchy <- hierarchy
  new_config$work_intensity <- work_intensity
  return(new_config)
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
