library(R6)

# AnalysisParameters ----
#' Analysis parameters
#' @description Analysis parameters
#' @keywords internal
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
#'
#' @return AnalysisConfiguration with reduced hieararchy
#' @examples
#'
#' conf <- LFQServiceData::skylinePRMSampleData_A$config_f()
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
#' @examples
#' xx <- data.frame(A = c("a","a","a"), B = c("d","d","e"))
#' # debug(make_interaction_column)
#' x <- make_interaction_column(xx, c("B","A"))
#' x <- make_interaction_column(xx, c("A"))
#'
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' config <- bb$config_f()
#' analysis <- bb$analysis(bb$data, bb$config_f())
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
#' @examples
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
#' @export
#' @param data data.frame
#' @param config AnlalysisConfiguration
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
#' @keywords internal
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
#' @examples
#'
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' config <- bb$config_f()
#' data <- bb$analysis(bb$data, config)
#' dim(data)
#' xx <- complete_cases(data, config)
#' stopifnot(nrow(data) < nrow(xx))
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


# Functions - summary ----

# Functions - summarize factors ----

#' table of distinct factors (sample annotation)
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#'
#' @export
#' @keywords internal
#' @examples
#' library(tidyverse)
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' config <- bb$config_f()
#' data <- bb$analysis(bb$data, config)
#' xx <- table_factors(data,config )
#' xx %>% dplyr::group_by(!!sym(config$table$factorKeys())) %>%
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
#' @examples
#' library(tidyverse)
#' library(LFQService)
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' config <- bb$config_f()
#' data <- bb$analysis(bb$data, config)
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
#' @examples
#' library(LFQService)
#' library(tidyverse)
#'
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' config <- bb$config_f()
#' data <- bb$analysis(bb$data, config)
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
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param hierarchy for which hierarchy level (default up to hierarchy depth)
#' @param factors which factors to include
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
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
#' @examples
#'
#' library(LFQService)
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
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

# Functions - Missigness ----

#' compute missingness statistics per hierarchy and factor level
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param factors factor to include (default up to factor depth)
#' @param hierarchy hierarchy to include (default up to hierarchy depth)
#' @param workIntensity work intensity column
#' @export
#' @keywords internal
#' @examples
#' library(tidyverse)
#' library(LFQService)
#'
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#'
#' configur$parameter$qVal_individual_threshold <- 0.01
#' xx <- LFQService::removeLarge_Q_Values(data,
#'    configur)
#' xx <- complete_cases(xx, configur)
#' nrow(xx)
#' x <- interaction_missing_stats(xx, configur)$data %>% arrange(desc(nrNAs))
#' stopifnot(nrow(x) == 7416)
#' stopifnot(sum(is.na(x$meanArea)) == 249)
#' stopifnot(length(unique(x$protein_Id)) == 37)
#'
#' tmp <- interaction_missing_stats(xx, configur,
#'  factors= character(),
#'   hierarchy = configur$table$hierarchyKeys()[1])$data
#' stopifnot(nrow(tmp) == 37)
#'
#' tmp <- interaction_missing_stats(xx, configur,
#'   hierarchy = configur$table$hierarchyKeys()[1])$data
#' stopifnot(sum(is.na(tmp$nrMeasured))==0)
#'
interaction_missing_stats <- function(pdata,
                                      config,
                                      factors = config$table$fkeysDepth(),
                                      hierarchy = config$table$hierarchyKeys(),
                                      workIntensity = config$table$getWorkIntensity())
{
  pdata <- complete_cases(pdata, config)
  table <- config$table
  missingPrec <- pdata %>% group_by_at(c(factors,
                                     hierarchy,
                                     table$isotopeLabel
  ))
  missingPrec <- missingPrec %>%
    dplyr::summarize(nrReplicates = n(),
                     nrNAs = sum(is.na(!!sym(workIntensity))),
                     meanArea = mean(!!sym(workIntensity), na.rm = TRUE)) %>%
    mutate(nrMeasured = nrReplicates - nrNAs) %>% dplyr::ungroup()
  return(list(data = missingPrec,
              summaries = c("nrReplicates","nrNAs","nrMeasured","meanArea")))
}

#' Compute interaction averages and
#' impute data using mean of lowest 0.1 (default)
#'
#' used in Acetylation project p2916
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param factors factor to include (default up to factor depth)
#' @param probs quantile to take average from (default 0.1)
#' @return function with parameter `value`
#' `c("long", "nrReplicates", "nrMeasured", "meanArea", "imputed", "allWide", "all")`
#' @export
#' @keywords internal
#' @return function
#' @examples
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#' configur$parameter$qVal_individual_threshold <- 0.01
#'
#' xx <- LFQService::removeLarge_Q_Values(data, configur)
#' xx <- complete_cases(xx, configur)
#'
#' tmp <- interaction_missing_stats(xx, configur)
#' head(tmp)
#' fun <- missigness_impute_interactions(xx, configur)
#' dd <- fun("long")
#' xx <- fun(DEBUG=TRUE)
#' sum(is.na(xx$long$nrReplicates))
#' xxx <- (fun("nrReplicates"))
#' head(xxx)
#' xxx <- fun("all")
#' head(xxx)
#'
missigness_impute_interactions <- function(pdata,
                                           config,
                                           factors = config$table$fkeysDepth(),
                                           probs = 0.1){
  mstats <- interaction_missing_stats(pdata, config, factors = factors)
  x_summaries <- mstats$summaries
  xx <- mstats$data
  xx <- make_interaction_column(xx, factors, sep = ":")


  lowerMean <- function(meanArea, probs = probs){
    meanAreaNotNA <- na.omit(meanArea)
    small10 <- meanAreaNotNA[meanAreaNotNA < quantile(meanAreaNotNA, probs = probs)]
    meanArea[is.na(meanArea)] <- mean(small10)
    return(meanArea)
  }

  xx <- xx %>%
    group_by(interaction) %>%
    mutate(imputed = lowerMean(meanArea,probs = 0.2))

  res_fun <- function(value = c("long",
                                "nrReplicates",
                                "nrMeasured",
                                "meanArea",
                                "imputed",
                                "allWide",
                                "all" ),
                      add.prefix = TRUE,
                      DEBUG = FALSE){
    value <- match.arg(value)
    if (DEBUG) {
      return(list(value = value, long = xx , config = config ))
    }

    if (value == "long") {
      return(xx)
    }else{
      xx <- xx %>% dplyr::select(-one_of(factors))

      pid <- config$table$hkeysDepth()
      nrReplicates <- xx %>%
        dplyr::select( -one_of(c(setdiff(x_summaries,"nrReplicates"),"imputed") )) %>%
        tidyr::spread(interaction, nrReplicates, sep = ".nrReplicates.") %>%
        arrange(!!!syms(pid)) %>%
        dplyr::ungroup()
      nrMeasured <- xx %>% dplyr::select(-one_of(c(setdiff(x_summaries,"nrMeasured"),"imputed" ) )) %>%
        tidyr::spread(interaction, nrMeasured, sep = ".nrMeasured.") %>%
        arrange(!!!syms(pid)) %>% dplyr::ungroup()

      meanArea <- xx %>% dplyr::select(-one_of(c(setdiff(x_summaries,"meanArea"),"imputed" ) )) %>%
        tidyr::spread(interaction, meanArea, sep = ".meanArea.") %>%
        arrange(!!!syms(pid)) %>% dplyr::ungroup()

      meanAreaImputed <- xx %>% dplyr::select(-one_of(setdiff(x_summaries,"imputed" ) )) %>%
        tidyr::spread(interaction, imputed, sep = ".imputed.") %>%
        arrange(!!!syms(pid)) %>% dplyr::ungroup()

      allTables <- list(meanArea = meanArea,
                        nrMeasured = nrMeasured,
                        nrReplicates = nrReplicates,
                        meanAreaImputed = meanAreaImputed)

      if (value == "all") {
        allTables[["long"]] <- xx
        return(allTables)
      }else if (value == "allWide") {
        return(purrr::reduce(allTables, inner_join))
      }else if (value == "nrReplicates") {
        srepl <- if (add.prefix) {"nrRep."}else{""}
        colnames(nrReplicates) <- gsub("interaction.nrReplicates.", srepl ,colnames(nrReplicates))
        nrReplicates <- tibble::add_column( nrReplicates, "value" = value, .before = 1)
        return(nrReplicates)
      }else if (value == "nrMeasured") {
        srepl <- if (add.prefix) {"nrMeas."}else{""}
        colnames(nrMeasured) <- gsub("interaction.nrMeasured.", srepl ,colnames(nrMeasured))
        nrMeasured <- tibble::add_column( nrMeasured, "value" = value, .before = 1)
        return(nrMeasured)
      }else if (value == "meanArea") {
        srepl <- if (add.prefix) {"mean."}else{""}
        colnames(meanArea) <- gsub("interaction.meanArea.", srepl ,colnames(meanArea))
        meanArea <- tibble::add_column( meanArea, "value" = value, .before = 1)
        return(meanArea)
      }else if (value == "imputed") {
        srepl <- if (add.prefix) {"mean.imp."}else{""}
        colnames(meanAreaImputed) <- gsub("interaction.imputed.", srepl ,colnames(meanAreaImputed))
        meanAreaImputed <- tibble::add_column( meanAreaImputed, "value" = value, .before = 1)
        return(meanAreaImputed)
      }
    }
  }

  #  nrMeasured %>% dplyr::select(starts_with("interaction")) -> nrMeasuredM
  #  nrReplicates %>% dplyr::select(starts_with("interaction")) -> nrReplicatesM
  return(res_fun)
}


#' compute per group averages and impute values
#' should generalize at some stage
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param probs quantile to take average from (default 0.1)
#' @param value use default
#' @param add.prefix use default
#'
#' @export
#' @keywords internal
#' @examples
#'
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#' configur$parameter$qVal_individual_threshold <- 0.01
#' xx <- LFQService::removeLarge_Q_Values(data, configur)
#' xx <- complete_cases(xx, configur)
#' #undebug(missigness_impute_factors_interactions)
#' fun <- missigness_impute_factors_interactions(xx, configur)
#' head(fun)
#' fun <- missigness_impute_factors_interactions(xx, configur, value = "imputed")
#' head(fun)
#' dim(fun)
#' dim(dplyr::distinct(fun[,1:6]))
#' fun <- missigness_impute_factors_interactions(xx, configur, value = "nrMeasured")
missigness_impute_factors_interactions <-
  function(pdata,
           config,
           probs = 0.1,
           value = c("nrReplicates", "nrMeasured", "meanArea", "imputed"),
           add.prefix = FALSE){
    value <- match.arg(value)
    fac_fun <- list()
    fac_fun[["interaction"]] <- missigness_impute_interactions(pdata,
                                                               config,
                                                               probs = probs)
    if (config$table$factorDepth > 1 ) { # if 1 only then done
      for (factor in config$table$fkeysDepth()) {
        fac_fun[[factor]] <- missigness_impute_interactions(pdata,
                                                            config,
                                                            factors = factor,
                                                            probs = probs)
      }
    }
    fac_res <- list()
    for (fun_name in names(fac_fun)) {
      fac_res[[fun_name]] <- fac_fun[[fun_name]](value, add.prefix = add.prefix)
    }
    intfact <- purrr::reduce(fac_res,
                             dplyr::inner_join,
                             by = c(config$table$hierarchyKeys(),
                                    config$table$isotopeLabel, "value"))
    return(intfact)
  }



#' Compute fold changes given Contrasts
#' @keywords internal
#' @export
#'
#' @examples
#'
#' library(LFQService)
#' library(tidyverse)
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#' configur$parameter$qVal_individual_threshold <- 0.01
#' data <- LFQService::removeLarge_Q_Values(data, configur)
#' data <- complete_cases(data, configur)
#'
#' Contrasts <- c("TimeT168vsT2" = "TimeT168 - TimeT2","TimeT168vsT24" = "TimeT168 - TimeT24" )
#' message("missigness_impute_factors_interactions : imputed")
#' xx <- missigness_impute_factors_interactions(data, configur, value = "nrMeasured" )
#' imputed <- missigness_impute_contrasts(xx, configur, Contrasts)
#'
#' xx <- missigness_impute_factors_interactions(data, configur, value = "imputed" )
#' imputed <- missigness_impute_contrasts(xx, configur, Contrasts)
#' xx <- missigness_impute_factors_interactions(data, configur, value = "meanArea" )
#' message("missigness_impute_factors_interactions : meanArea")
#' mean <- missigness_impute_contrasts(xx, configur, Contrasts)
#'
missigness_impute_contrasts <- function(data,
                                        config,
                                        contrasts,
                                        agg_fun = function(x){median(x, na.rm = TRUE)})
{
  for (i in 1:length(contrasts)) {
    message(names(contrasts)[i], "=", contrasts[i],"\n")
    data <- dplyr::mutate(data, !!names(contrasts)[i] := !!rlang::parse_expr(contrasts[i]))
  }

  if (!is.null(agg_fun)) {
    data <- data %>% group_by_at(c("value" , config$table$hkeysDepth())) %>%
      summarise_if(is.numeric, agg_fun)
  }
  return(data)
}

#' Compute fold changes given Contrasts 2
#' @keywords internal
#' @export
#'
workflow_missigness_impute_contrasts <- function(data,
                                                 config,
                                                 contrasts){

  xx <- missigness_impute_factors_interactions(data, config, value = "imputed" )
  message("missigness_impute_factors_interactions : imputed")
  imputed <- missigness_impute_contrasts(xx, config, contrasts)
  xx <- missigness_impute_factors_interactions(data, config, value = "meanArea" )
  message("missigness_impute_factors_interactions : meanArea")
  mean <- missigness_impute_contrasts(xx, config, contrasts)

  dd <- dplyr::bind_rows(imputed, mean)
  dd_long <- dd %>% tidyr::gather("contrast","int_val",
                                  colnames(dd)[sapply(dd, is.numeric)])

  res_fun <- function(value = c("long", "wide","raw"),
                      what = c("contrasts", "factors", "all"),
                      DEBUG = FALSE){
    value <- match.arg( value )
    what  <- match.arg( what  )
    if (DEBUG) {
      return(list(value = value,
                  what = what,
                  dd_long = dd_long,
                  contrasts = contrasts,
                  config = config ))
    }

    if (what == "contrasts") {
      dd_long <- dplyr::filter(dd_long, contrast %in% names(contrasts))
    }else if (what == "factors") {
      dd_long <- dplyr::filter(dd_long, !contrast %in% names(contrasts))
    }else if (what == "all") {
    }

    if (value == "long") {
      long_xxxx <- dd_long %>% spread(value, int_val)
      return(long_xxxx)
    }else if (value == "wide") {
      dd <- dd_long %>% unite(contrast.v , value, contrast, sep="~") %>% spread(contrast.v, int_val)
      xxx_imputed <- inner_join(LFQService::summarize_hierarchy(data,config),dd)
      return(xxx_imputed)
    }else if (value == "raw") {
      return(dd_long)
    }
  }
  return(res_fun)
}
#' Histogram summarizing missigness
#' @export
#' @keywords internal
#' @family plotting
#' @examples
#' library(tidyverse)
#' library(LFQService)
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#'
#' xx <- complete_cases(data, configur)
#' configur$parameter$qVal_individual_threshold <- 0.01
#' xx <- LFQService::removeLarge_Q_Values(data, configur)
#' xx <- complete_cases(xx, configur)
#' missigness_histogram(xx, configur)
#'
#' missingPrec <- interaction_missing_stats(xx, configur)
#'
#' setNa <- function(x){ifelse(x < 100, NA, x)}
#' data %>% dplyr::mutate(Area = setNa(Area)) -> data
#' missigness_histogram(data, configur)
#'
missigness_histogram <- function(x, config, showempty = TRUE, factors = config$table$fkeysDepth()){
  table <- config$table
  missingPrec <- interaction_missing_stats(x, config , factors)$data
  missingPrec <- missingPrec %>%  dplyr::ungroup() %>% dplyr::mutate(nrNAs = as.factor(nrNAs))

  if (showempty) {
    if (config$table$is_intensity_transformed) {
      missingPrec <- missingPrec %>% dplyr::mutate(meanArea = ifelse(is.na(meanArea),1,meanArea))
    }else{
      missingPrec <- missingPrec %>% dplyr::mutate(meanArea = ifelse(is.na(meanArea),-20,meanArea))
    }

  }

  factors <- table$fkeysDepth()
  formula <- paste(table$isotopeLabel, "~", paste(factors, collapse = "+"))
  message(formula)

  p <- ggplot(missingPrec, aes(x = meanArea, fill = nrNAs, colour = nrNAs)) +
    geom_histogram(alpha = 0.2, position = "identity") +
    facet_grid(as.formula(formula)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  if (!config$table$is_intensity_transformed) {
    p <- p + scale_x_log10()
  }
  p
}

#' cumulative sums of missing
#' @export
#' @keywords internal
#' @family plotting
#' @examples
#'
#' setNa <- function(x){ifelse(x < 100, NA, x)}
#'
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#'
#' data %>% dplyr::mutate(Area = setNa(Area)) -> data
#' res <- missingness_per_condition_cumsum(data,configur)
#' stopifnot("ggplot" %in% class(res$figure))
#' print(res$figure)
#' res$data
missingness_per_condition_cumsum <- function(x,
                                             config,
                                             factors = config$table$fkeysDepth()){
  table <- config$table
  missingPrec <- interaction_missing_stats(x, config,factors)$data

  xx <- missingPrec %>% group_by_at(c(table$isotopeLabel, factors,"nrNAs","nrReplicates")) %>%
    dplyr::summarize(nrTransitions = n())

  xxcs <- xx %>% group_by_at( c(table$isotopeLabel,factors)) %>% arrange(nrNAs) %>%
    dplyr::mutate(cumulative_sum = cumsum(nrTransitions))
  res <- xxcs  %>% dplyr::select(-nrTransitions)

  formula <- paste(table$isotopeLabel, "~", paste(factors, collapse = "+"))
  message(formula)

  nudgeval = mean(res$cumulative_sum) * 0.05
  p <- ggplot(res, aes(x = nrNAs, y = cumulative_sum)) +
    geom_bar(stat = "identity", color = "black", fill = "white") +
    geom_text(aes(label = cumulative_sum), nudge_y = nudgeval, angle = -45) +
    facet_grid(as.formula(formula))

  res <- res %>% tidyr::spread("nrNAs","cumulative_sum")
  return(list(data = res, figure = p))
}

#' Summarize missing in condtion as barplot
#' @export
#' @keywords internal
#' @family plotting
#' @examples
#' setNa <- function(x){ifelse(x < 100, NA, x)}
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#' data %>% dplyr::mutate(Area = setNa(Area)) -> data
#' res <- missingness_per_condition(data, configur)
#' names(res)
#' stopifnot(c(8,7) == dim(res$data))
#' stopifnot("ggplot" %in% class(res$figure))
#' print(res$figure)
#'
missingness_per_condition <- function(x, config, factors = config$table$fkeysDepth()){
  table <- config$table
  missingPrec <- interaction_missing_stats(x, config, factors)$data
  hierarchyKey <- tail(config$table$hierarchyKeys(),1)
  hierarchyKey <- paste0("nr_",hierarchyKey)
  xx <- missingPrec %>% group_by_at(c(table$isotopeLabel,
                                      factors,"nrNAs","nrReplicates")) %>%
    dplyr::summarize( !!sym(hierarchyKey) := n())

  formula <- paste(table$isotopeLabel, "~", paste(factors, collapse = "+"))
  #message(formula)

  nudgeval = max(xx[[hierarchyKey]]) * 0.05

  p <- ggplot(xx, aes_string(x = "nrNAs", y = hierarchyKey)) +
    geom_bar(stat = "identity", color = "black", fill = "white") +
    geom_text(aes(label = !!sym(hierarchyKey)), nudge_y = nudgeval, angle = 45) +
    facet_grid(as.formula(formula))
  xx <- xx %>% tidyr::spread("nrNAs",hierarchyKey)

  return(list(data = xx ,figure = p))
}

# Functions - Handling isotopes ----

#' spreads isotope label heavy light into two columns
#' @export
#' @keywords internal
#' @examples
#' library(LFQService)
#' setNa <- function(x){ifelse(x < 100, NA, x)}
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#'
#' data %>% dplyr::mutate(Area = setNa(Area)) -> data
#' x<-spreadValueVarsIsotopeLabel(data,configur)
#' head(x)
#'
#' bb <- LFQServiceData::skylineSRM_HL_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#' conf <- LFQServiceData::skylineconfig_HL$clone(deep=TRUE)
#' x<-spreadValueVarsIsotopeLabel(LFQServiceData::sample_analysis_HL, conf)
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


