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
    #' @field opt_se optional column with standard errors (e.g. from limpa aggregation)
    opt_se = character(),
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
          "opt_se",
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
        self$opt_se,
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
