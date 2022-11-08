# AnalysisTableAnnotation ----
#'
#' Create Annotation
#' @description
#' Annotates Data Table
#' @family configuration
#' @export
#' @examples
#'
#' ata <- AnalysisTableAnnotation$new()
#' ata$fileName <- "rawfile.column"
#' ata$hierarchy[["protein"]] <- "protein.column"
#' ata$factors[["explanatory"]] <- "explanatory.column"
#' ata$set_response("abundance")
#' ata$id_required()
#' ata$id_vars()
#' ata$value_vars()
#' ata$annotation_vars()
#' ac <- AnalysisConfiguration$new(ata)
#'
AnalysisTableAnnotation <- R6::R6Class(
  "AnalysisTableAnnotation",
  public = list(
    #' @description
    #' create a new  AnalysisTableAnnotationG
    initialize = function() {
    },

    #' @field fileName column name of column containing raw file names
    fileName = NULL,
    #' @field sampleName (will be generated from factors)
    sampleName = "sampleName",

    #' @field isotopeLabel which column contains the isotope label (e.g. heavy or light), Gor light only if LFQ.
    isotopeLabel = "isotopeLabel",
    # do you want to model charge sequence etc?
    #' @field ident_qValue column name with identification QValues (smaller better)
    ident_qValue = character(),
    #' @field ident_Score column with identification score (lager better)
    ident_Score = character(),

    #' @field opt_rt optional column with rt information
    opt_rt = character(),
    #' @field opt_mz optional column with mz information
    opt_mz = character(),


    #' @field workIntensity column which contains the intensities
    workIntensity = NULL, # could be list with names and functions
    #' @field is_response_transformed are the intensities transformed for constant variance
    is_response_transformed = FALSE,
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


    #' @field factors Names of columns containing factors (annotations)
    factors = list(), # ordering is important - first is considered the main
    #' @field factorDepth number of relevant factors (used by plotting functions etc)
    factorDepth = 1,
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

    #' @field hierarchy list with columns describing the measurement hierarchy (i.e. protein peptide precursor fragment)
    hierarchy = list(),
    #' @field hierarchyDepth At which depth do you want to model i.e. i.e. protein than 1
    hierarchyDepth = 1,
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
    #' get hierarchy keys
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
    #' get hierarchy keys up to depth
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
        self$sampleName
      )
      return(id_vars)
    },
    #' @description
    #' get names of columns containing observations e.g. (intensity, qValue, mz or rt)
    value_vars = function() {
      "Columns containing values"
      valueVars <- c(self$get_response(), self$ident_qValue, self$ident_Score, self$opt_mz, self$opt_rt)
      return(valueVars)
    },
    #' @description
    #' get names of columns with sample annotations
    #'
    annotation_vars = function() {
      annotationVars <- c(self$fileName, self$sampleName, self$factor_keys())
      return(annotationVars)
    }
  )
)
