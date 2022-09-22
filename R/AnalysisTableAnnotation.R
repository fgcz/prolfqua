# AnalysisTableAnnotation ----
#'
#' Create Annotation
#' @description
#' Annotates Data Table
#' @keywords internal
#' @family configuration
#' @export
#'
AnalysisTableAnnotation <- R6::R6Class(
  "AnalysisTableAnnotation",
  public = list(
    #' @description
    #' create a new  AnalysisTableAnnotationG
    initialize = function(){
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
    opt_rt  = character(),
    #' @field opt_mz optional column with mz information
    opt_mz = character(),


    #' @field workIntensity column which contains the intensities
    workIntensity = NULL, # could be list with names and functions
    #' @field is_intensity_transformed are the intensities transformed for constant variance
    is_intensity_transformed = FALSE,
    #' @description
    #' add name of intensity column
    #' @param colName name of intensity column
    setWorkIntensity = function(colName){
      self$workIntensity <- c(self$workIntensity, colName)
    },
    #' @description
    #' get name of working intensity column
    getWorkIntensity = function(){
      return(tail(self$workIntensity, n = 1))
    },
    #' @description
    #' remove last name in array of working intensity column names
    popWorkIntensity = function(){
      res <- self$workIntensity[length(self$workIntensity)]
      self$workIntensity <- self$workIntensity[-length(self$workIntensity)]
      return(res)
    },


    #' @field factors names of columns containing factors (annotions)
    factors = list(), # ordering is important - first is considered the main
    #' @field factorDepth facet plot according to the first or the first two factors (factorDepth can be 1 or 2)
    factorDepth = 1,
    #' @description
    #' get factor keys
    #' @return array with keys
    factorKeys = function(){
      return(names(self$factors))
    },
    #' @description
    #' get factor keys till factorDepth
    fkeysDepth = function(){
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
    hierarchyKeys = function(rev = FALSE){
      if (rev) {
        return(rev(names(self$hierarchy)))
      }else{
        return(names(self$hierarchy))
      }
    },
    #' @description
    #' get hierarchy keys up to depth
    #' @param names if TRUE names only if FALSE key value pairs
    #' @return array of column names
    hkeysDepth = function(names = TRUE){
      res <- head( self$hierarchy,n = self$hierarchyDepth)
      res <- if (names) {
        names(res)
      }else{
        res
      }
      return(res)
    },

    #' @description
    #' Id Columns which must be in the input data frame
    #' @return character array
    idRequired = function(){
      idVars <- c(
        self$fileName,
        purrr::map_chr(self$factors,"colnames"),
        unlist(self$hierarchy),
        self$isotopeLabel
      )
      return(idVars)
    },
    #' @description
    #' get names of columns annotating values (e.g. intensities)
    #' @return character array
    idVars = function(){
      "Id Columns which must be in the output data frame"
      idVars <- c(
        self$fileName,
        names(self$factors),
        names(self$hierarchy),
        self$isotopeLabel,
        self$sampleName)
      return(idVars)
    },
    #' @description
    #' get names of columns containing observations e.g. (intensity, qValue, mz or rt)
    valueVars = function(){
      "Columns containing values"
      valueVars <- c( self$getWorkIntensity(), self$ident_qValue, self$ident_Score, self$opt_mz, self$opt_rt)
      return(valueVars)
    },
    #' @description
    #' get names of columns with sample annotations
    #'
    annotationVars = function(){
      annotationVars <- c(self$fileName, self$sampleName, self$factorKeys() )
      return(annotationVars)
    }
  )
)
