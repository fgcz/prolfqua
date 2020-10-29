library(R6)

# AnalysisParameters ----
#' Analysis parameters
#' @description Analysis parameters
#' @keywords internal
#' @family configuration
#' @export
AnalysisParameters <- R6::R6Class(
  "AnalysisParameters",
  public = list(
    #' @field qVal_individual_threshold qValue threshold for sample
    qVal_individual_threshold  = 0.05,
    #' @field qVal_experiment_threshold qValue threshold for dataset
    qVal_experiment_threshold = 0.01,
    #' @field qVal_minNumber_below_experiment_threshold how many samples need to meet qVal_experiment_threshold
    qVal_minNumber_below_experiment_threshold = 3,
    #' @field min_nr_of_notNA minimum number of not NA's in all samples default 1
    min_nr_of_notNA = 1, # how many values per transition total
    #' @field min_nr_of_notNA_condition minimum number of not NA's in interaction
    min_nr_of_notNA_condition = 0, # how many not missing in condition
    #' @field min_peptides_protein minimum number of peptides per protein
    min_peptides_protein = 2
  )
)

# AnalysisTableAnnotation ----



#'
#' Create Annotation
#' @description Annotates Data Table
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
    isotopeLabel = character(),
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

# AnalysisConfiguration ----
#' Analysis Configuration
#' @description Hello world
#' @keywords internal
#' @family configuration
#' @export
#'
AnalysisConfiguration <- R6::R6Class("AnalysisConfiguration",
                                     public = list(
                                       #' @field sep separator to use when uniting columns is necessary
                                       sep = "~",
                                       #' @field table AnalysisTableAnnotation
                                       table = NULL,
                                       #' @field parameter AnalysisParameter
                                       parameter = NULL,
                                       #' @description
                                       #' create AnalysisConfiguration
                                       #' @param analysisTableAnnotation instance of AnalysisTableAnnotation
                                       #' @param analysisParameter instance of AnalysisParameter
                                       initialize = function(analysisTableAnnotation, analysisParameter){
                                         self$table <- analysisTableAnnotation
                                         self$parameter <- analysisParameter
                                       }
                                     )
)

#' Make reduced hierarchy configuration
#' @export
#' @keywords internal
#' @param config AnalysisConfiguration
#' @param workIntensity work intensity column
#' @param hierarchy new reduced hierarchy
#' @family configuration
#' @return AnalysisConfiguration with reduced hieararchy
#' @examples
#'
#' bb <- LFQService::ionstar$filtered()
#' conf <- bb$config$clone(deep=TRUE)
#' analysis <- bb$data
#'
#' make_reduced_hierarchy_config(conf,
#'  "testintensity",
#'  conf$table$hierarchy[1:2])
#'
make_reduced_hierarchy_config <- function(config, workIntensity , hierarchy ){
  newConfig <- config$clone(deep = TRUE)
  newConfig$table$hierarchy = hierarchy
  newConfig$table$workIntensity = workIntensity
  return(newConfig)
}

#' create interaction column from factors
#' @export
#' @keywords internal
#' @family configuration
#' @examples
#' xx <- data.frame(A = c("a","a","a"), B = c("d","d","e"))
#' # debug(make_interaction_column)
#' x <- make_interaction_column(xx, c("B","A"))
#' x <- make_interaction_column(xx, c("A"))
#'
#' bb <- LFQService::ionstar$filtered()
#' config <- bb$config$clone(deep=TRUE)
#' analysis <- bb$data
#'
#' config$table$factorKeys()
#' config$table$factorDepth <- 1
#' make_interaction_column(analysis,
#'    config$table$fkeysDepth())
#'
make_interaction_column <- function(data, columns, sep="."){
  intr <- dplyr::select(data, columns)
  intr <- purrr::map_dfc(intr, factor)
  names(columns) <- columns
  newlev <- purrr::map2(columns, intr, function(x,y){paste0(x,levels(y))})
  intr <- purrr::map2_dfc(columns, intr, paste0)
  intr <- purrr::map2_dfc(intr , newlev, forcats::fct_relevel)

  colnames(intr) <- paste0("interaction_",columns)
  colname <- "interaction"
  data <- data %>% dplyr::mutate(!!colname := interaction(intr, sep = sep))
  return(data)
}


# Functions - Configuration ----
#' Extract all value slots in an R6 object
#' @param r6class r6 class
#' @keywords internal
#' @family configuration
#' @export
R6extractValues <- function(r6class){
  tmp <- sapply(r6class, class)
  slots <- tmp[ !tmp %in% c("environment", "function")]
  res <- list()
  for (i in names(slots)) {
    if ("R6" %in% class(r6class[[i]])) {
      res[[i]]  <- R6extractValues(r6class[[i]])
    }else{
      res[[i]] <- r6class[[i]]
    }
  }
  return(res)
}



#' Extracts columns relevant for a configuration from a data frame
#' and create new columns i sampleName column etc.
#' @param data data.frame
#' @param config AnlalysisConfiguration
#' @param cc complete cases default TRUE
#' @export
#' @keywords internal
#' @family configuration
#' @examples
#'
#'
#' skylineconfig <- create_config_Skyline(isotopeLabel = "Isotope.Label.Type",
#'  ident_qValue = "Detection.Q.Value")
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#'
#' sample_analysis <- setup_analysis(LFQServiceData::skylinePRMSampleData_A$data, skylineconfig)
#'
setup_analysis <- function(data, configuration, cc = TRUE ){
  table <- configuration$table
  for (i in 1:length(table$hierarchy))
  {
    data <- tidyr::unite(data, UQ(sym(table$hierarchyKeys()[i])), table$hierarchy[[i]],remove = FALSE, sep = configuration$sep)
  }
  data <- dplyr::select(data , -dplyr::one_of(dplyr::setdiff(unlist(table$hierarchy), table$hierarchyKeys() )))

  for (i in 1:length(table$factors))
  {
    if ( length(table$factors[[i]]) > 1) {
      data <- tidyr::unite(data, UQ(sym(table$factorKeys()[i])), table$factors[[i]],remove = FALSE, sep = configuration$sep)
    } else {
      newname <- table$factorKeys()[i]
      data <- dplyr::mutate(data, !!newname := !!sym(table$factors[[i]]))
    }
  }

  sampleName <- table$sampleName

  if (!sampleName  %in% names(data)) {
    message("creating sampleName")

    data <- data %>%  tidyr::unite( UQ(sym( sampleName)) , unique(unlist(table$factors)), remove = TRUE , sep = configuration$sep) %>%
      dplyr::select(sampleName, table$fileName) %>% dplyr::distinct() %>%
      dplyr::mutate_at(sampleName, function(x){ x <- make.unique( x, sep = configuration$sep )}) %>%
      dplyr::inner_join(data, by = table$fileName)
  } else {
    warning(sampleName, " already exists")
  }

  data <- data %>% dplyr::select(-dplyr::one_of(dplyr::setdiff(unlist(table$factors), table$factorKeys())))

  # Make implicit NA's explicit
  data <- data %>% dplyr::select(c(configuration$table$idVars(),configuration$table$valueVars()))

  if (cc) {
    data <- complete_cases( data , configuration)
  }
  return( data )
}

#' separates hierarchy columns into starting columns
#'
#' @export
#' @param data data.frame
#' @param config AnlalysisConfiguration
#' @family configuration
#' @keywords internal
separate_hierarchy <- function(data, config){
  for (hkey in config$table$hkeysDepth()) {
    data <- data %>% tidyr::separate( hkey, config$table$hierarchy[[hkey]], sep = config$sep, remove = FALSE)
  }
  return(data)
}

#' separates factor columns into starting columns
#' @param data data.frame
#' @param config AnlalysisConfiguration
#' @export
#'
#' @keywords internal
#' @family configuration
separate_factors <- function(data, config) {
  for (fkey in config$table$factorKeys()) {
    data <- data %>% tidyr::separate( fkey, config$table$factors[[fkey]], sep = config$sep, remove = FALSE)
  }
  return(data)
}


#' Complete cases
#' @export
#' @param pdata data.frame
#' @param config AnlalysisConfiguration
#' @keywords internal
#' @family configuration
#' @examples
#'
#' bb <- LFQService::ionstar$filtered()
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' xx <- complete_cases(data, config)
#' stopifnot(nrow(data) <= nrow(xx))
#'
complete_cases <- function(pdata, config) {
  message("completing cases")
  fkeys <- c(config$table$fileName,config$table$sampleName, config$table$factorKeys())
  hkeys <- c(config$table$isotopeLabel, config$table$hierarchyKeys())
  res <- tidyr::complete(
    pdata,
    tidyr::nesting(!!!syms(fkeys)),
    tidyr::nesting(!!!syms(hkeys))
  )
  return(res)
}


#' sample subset of data
#' @export
#' @keywords internal
#' @family configuration
#'
sample_subset <- function(size, pdata, config){
  hk <- config$table$hkeysDepth()
  message("Sampling ", size, paste(hk, collapse = "," ) )
  hkdf <- pdata %>% select(all_of(hk)) %>% distinct() %>% sample_n(size = size)
  sdata <- inner_join(hkdf, pdata)
  return(sdata)
}

# Functions - summary ----

# Functions - summarize factors ----

#' table of distinct factors (sample annotation)
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#'
#' @export
#' @keywords internal
#' @family configuration
#' @examples
#'
#' library(tidyverse)
#' bb <- LFQService::ionstar$filtered()
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' xx <- table_factors(data,config )
#' xx
#' xx %>% dplyr::group_by(!!!syms(config$table$factorKeys())) %>%
#'  dplyr::summarize(n = n())
#'
table_factors <- function(pdata, configuration){
  factorsTab <- pdata %>% dplyr::select(c(configuration$table$fileName, configuration$table$sampleName, configuration$table$factorKeys())) %>%
    dplyr::distinct() %>%
    arrange(!!sym(configuration$table$sampleName))
  return(factorsTab)
}


# Functions - summarize hierarchies

#' Count distinct elements for each level of hierarchy
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#'
#'
#' @export
#' @keywords internal
#' @family summary
#' @examples
#' library(tidyverse)
#' library(LFQService)
#' bb <- LFQService::ionstar$filtered()
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' hierarchy_counts(data, config)
hierarchy_counts <- function(pdata, config){
  hierarchy <- config$table$hierarchyKeys()
  res <- pdata %>%
    dplyr::group_by_at(config$table$isotopeLabel) %>%
    dplyr::summarise_at( hierarchy, n_distinct )
  return(res)
}

#' Count distinct elements for each level of hierarchy per sample
#' @export
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#'
#' @keywords internal
#' @family summary
#' @examples
#'
#' library(LFQService)
#' library(tidyverse)
#'
#' bb <- LFQService::ionstar$filtered()
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' res <- hierarchy_counts_sample(data, config)
#' res()
#' res("long")
#' res("plot")
hierarchy_counts_sample <- function(pdata,
                                    configuration)
{
  hierarchy <- configuration$table$hierarchyKeys()
  summary <- pdata %>% dplyr::filter(!is.na(!!sym(configuration$table$getWorkIntensity() ))) %>%
    dplyr::group_by_at(c(configuration$table$isotopeLabel, configuration$table$sampleName)) %>%
    dplyr::summarise_at( hierarchy, n_distinct )

  res <- function(value = c("wide", "long", "plot")){
    value <- match.arg(value)
    if (value == "wide") {
      return(summary)
    }else{
      long <- summary %>% tidyr::gather("key",
                                        "nr",-dplyr::one_of(configuration$table$isotopeLabel,
                                                            configuration$table$sampleName))
      if (value == "long") {
        return(long)
      }else if (value == "plot") {
        nudgeval <- mean(long$nr) * 0.05
        ggplot2::ggplot(long, aes(x = sampleName, y = nr)) +
          geom_bar(stat = "identity", position = "dodge", colour = "black", fill = "white") +
          facet_wrap( ~ key, scales = "free_y", ncol = 1) +
          geom_text(aes(label = nr), nudge_y = nudgeval, angle = 45) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
      }
    }
  }
  return(res)
}



#' Summarize hierarchy counts
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
#' library(LFQService)
#' library(tidyverse)
#'
#' bb <- LFQService::ionstar$filtered()
#' configur <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' summarize_hierarchy(data, configur)
#' summarize_hierarchy(data, configur, factors = character())
#'
#' summarize_hierarchy(data, configur,
#'  hierarchy = configur$table$hkeysDepth() )
#' summarize_hierarchy(data, configur,
#'  hierarchy = NULL, factors = configur$table$fkeysDepth() )
#' configur$table$hierarchyDepth = 1
#' summarize_hierarchy(data, configur,
#'  factors = configur$table$fkeysDepth())
#' configur$table$hierarchyDepth = 2
#' summarize_hierarchy(data, configur)
#' configur$table$hierarchyDepth = 3
#' summarize_hierarchy(data, configur )
#' configur$table$hierarchyDepth = 4
#' summarize_hierarchy(data, configur )
#' #summarize_hierarchy(LFQServiceData::dataIonstarFilteredPep$data,
#' #LFQServiceData::dataIonstarFilteredPep$config)
#'
summarize_hierarchy <- function(pdata,
                                config,
                                hierarchy = config$table$hkeysDepth(),
                                factors = character())
{
  all_hierarchy <- c(config$table$isotopeLabel, config$table$hierarchyKeys() )

  precursor <- pdata %>% dplyr::select(factors, all_hierarchy) %>% dplyr::distinct()
  x3 <- precursor %>% dplyr::group_by_at(c(factors, hierarchy)) %>%
    dplyr::summarize_at( setdiff(all_hierarchy, hierarchy),
                         list( n = dplyr::n_distinct))
  return(x3)
}


#' Summarize Protein counts
#' Works for one isotope label level only
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @export
#' @keywords internal
#' @family summary
#' @examples
#'
#' library(LFQService)
#'
#' bb <- LFQService::ionstar$filtered()
#' configur <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' LFQService:::summarize_protein(data, configur)
#'
summarize_protein <- function(pdata, config ){
  warning("DEPRECATED use summarize_hierarchy instead")
  rev_hierarchy <- config$table$hierarchyKeys(TRUE)
  if (length(rev_hierarchy) < 4) {
    warning("only usable for 4 level hierarchies : protein_Id, peptide_Id, precursor_Id, fragment_Id")
    warning("your data has only: ", paste(config$table$hierarchyKeys(), collapse = ", ") )
  }

  precursorSum <- pdata %>% dplyr::select(rev_hierarchy) %>% dplyr::distinct() %>%
    group_by_at(rev_hierarchy[-1]) %>%
    dplyr::summarize(nrFragments = n())

  peptideSum <- precursorSum %>% group_by_at(rev_hierarchy[-(1:2)]) %>%
    dplyr::summarize(nrPrecursors = n(),
                     minNrFragments = min(nrFragments),
                     maxNrFragments = max(nrFragments))


  proteinSum <- peptideSum %>% group_by_at(rev_hierarchy[-(1:3)])  %>%
    dplyr::summarize(nrpeptides = n(),
                     minNrPrecursors = min(nrPrecursors),
                     maxNrPrecursors = max(nrPrecursors),
                     minNrFragments = min(minNrFragments),
                     maxNrFragments = max(maxNrFragments)
    )
  proteinPeptide <- proteinSum %>% tidyr::unite(Precursors ,minNrPrecursors , maxNrPrecursors, sep="-", remove=FALSE)
  proteinPeptide <- proteinPeptide %>% tidyr::unite(Fragments ,minNrFragments , maxNrFragments, sep="-", remove=FALSE)
  return(proteinPeptide)
}


# Functions - Handling isotopes ----

#' spreads isotope label heavy light into two columns
#' @export
#' @keywords internal
#' @family configuration
#' @examples
#' library(LFQService)
#'
#' bb <- LFQService::ionstar$filtered()
#' configur <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' x<-spreadValueVarsIsotopeLabel(data,configur)
#' head(x)
#'
#' bb <- LFQService::skylineSRM_HL_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#'
#' bb <- LFQService::skylineSRM_HL_A
#' conf <- bb$config_f()
#' analysis <- bb$analysis(bb$data, bb$config_f())
#' x <- spreadValueVarsIsotopeLabel(analysis, conf)
#' head(x[,5:ncol(x)])
#'
spreadValueVarsIsotopeLabel <- function(resData, config){
  table <- config$table
  idVars <- table$idVars()
  resData2 <- resData %>% dplyr::select(c(table$idVars(), table$valueVars()) )
  resData2 <- resData2 %>% tidyr::gather(variable, value, - idVars  )
  resData2 <- resData2 %>%  tidyr::unite(temp, table$isotopeLabel, variable )
  HLData <- resData2 %>% tidyr::spread(temp,value)
  invisible(HLData)
}

# Computing protein Intensity summaries ---


