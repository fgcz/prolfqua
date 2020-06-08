library(R6)

# AnalysisParameters ----
#' Analysis parameters
#' @description Analysis parameters
#' @keywords internal
#' @export
AnalysisParameters <- R6::R6Class("AnalysisParameters",
                                  public = list(
                                    #' @field qVal_individual_threshold qValue threshold for sample
                                    qVal_individual_threshold  = 0.05,
                                    #' @field qVal_experiment_threshold qValue threshold for dataset
                                    qVal_experiment_threshold = 0.01,
                                    #' @field qVal_minNumber_below_experiment_threshold how many samples need to meet qVal_experiment_threshold
                                    qVal_minNumber_below_experiment_threshold = 3,
                                    #' @field is_intensity_transformed is the intensity transformed (typically log2)
                                    is_intensity_transformed = FALSE, # important for some plotting functions
                                    #' @field min_nr_of_notNA minimum number of not NA's in all samples default 1
                                    min_nr_of_notNA = 1, # how many values per transition total
                                    #' @field min_nr_of_notNA_condition minimum number of not NA's in interaction
                                    min_nr_of_notNA_condition = 0, # how many not missing in condition
                                    #' @field min_peptides_protein minimum number of peptides per protein
                                    min_peptides_protein = 2
                                  )
)

# AnalysisTableAnnotation ----



#' Create Annotation
#' @description Annotates Data Table
#' @keywords internal
#' @export
AnalysisTableAnnotation <- R6::R6Class("AnalysisTableAnnotation",
                                       public = list(
                                         #' @field fileName some funny name
                                         fileName = NULL,
                                         factors = list(), # ordering is important - first is considered the main
                                         factorDepth = 1,

                                         sampleName = "sampleName",
                                         # measurement levels
                                         hierarchy = list(),
                                         hierarchyDepth = 1,
                                         isotopeLabel = character(),
                                         # do you want to model charge sequence etc?

                                         ident_qValue = character(), # smaller better
                                         ident_Score = character(), # larger better
                                         workIntensity = NULL, # could be list with names and functions

                                         opt_rt  = character(),
                                         opt_mz = character(),
                                         is_intensity_transformed = FALSE,

                                         #' create a new  AnalysisTableAnnotation
                                         initialize = function(){
                                         },
                                         #' get number of factor Levels
                                         getfactorDepth = function(){
                                           if (length(self$factorDepth) == 0) {
                                             return(length(self$factors))
                                           }else{
                                             return(self$factorDepth)
                                           }
                                         },
                                         #' set name of working intensity
                                         setWorkIntensity = function(colName){
                                           self$workIntensity <- c(self$workIntensity, colName)
                                         },
                                         #' get name of working intensity column
                                         getWorkIntensity = function(){
                                           return(tail(self$workIntensity, n = 1))
                                         },
                                         #' remove working intensity column
                                         popWorkIntensity = function(){
                                           res <- self$workIntensity[length(self$workIntensity)]
                                           self$workIntensity <- self$workIntensity[-length(self$workIntensity)]
                                           return(res)
                                         },
                                         #' Id Columns which must be in the input data frame
                                         idRequired = function(){
                                           idVars <- c(
                                             self$fileName,
                                             purrr::map_chr(self$factors,"colnames"),
                                             unlist(self$hierarchy),
                                             self$isotopeLabel
                                           )
                                           return(idVars)
                                         },
                                         hierarchyKeys = function(rev = FALSE){
                                           if (rev) {
                                             return(rev(names(self$hierarchy)))
                                           }else{
                                             return(names(self$hierarchy))
                                           }
                                         },
                                         hkeysDepth = function(names = TRUE){
                                           res <- head( self$hierarchy,n = self$hierarchyDepth)
                                           res <- if (names) {
                                             names(res)
                                           }else{
                                             res
                                           }
                                           return(res)
                                         },
                                         factorKeys = function(){
                                           return(names(self$factors))
                                         },
                                         fkeysDepth = function(){
                                           res <- head(self$factors, n = self$factorDepth)
                                           return(names(res))
                                         },
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
                                         valueVars = function(){
                                           "Columns containing values"
                                           valueVars <- c( self$getWorkIntensity(), self$ident_qValue, self$ident_Score, self$opt_mz, self$opt_rt)
                                           return(valueVars)
                                         },
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
                                       #' @field project_Id the project Id
                                       project_Id = "",
                                       order_Id = "",
                                       workunit_Id = "",
                                       sep = "~",
                                       table = NULL,
                                       parameter = NULL,
                                       initialize = function(analysisTableAnnotation, analysisParameter){
                                         self$table <- analysisTableAnnotation
                                         self$parameter <- analysisParameter
                                       }
                                     )
)

#' Make reduced hierarchy configuration
#' @export
#' @keywords internal
#' @examples
#' make_reduced_hierarchy_config(LFQServiceData::skylineconfig, "testintensity", LFQServiceData::skylineconfig$table$hierarchy[1:2])
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
#' x <- make_interaction_column(xx, c("B","A"))
#' x <- make_interaction_column(xx, c("A"))
make_interaction_column <- function(data, columns, sep="."){
  intr <- dplyr::select(data, columns)
  intr <- purrr::map_dfc(intr, factor)

  newlev <- purrr::map2(columns, intr, function(x,y){paste0(x,levels(y))})
  intr <- purrr::map2_dfc(columns, intr, paste0)
  intr <- purrr::map2_dfc(intr , newlev, fct_relevel)

  colnames(intr) <- paste0("interaction_",columns)
  colname <- "interaction"
  data <- data %>% dplyr::mutate(!!colname := interaction(intr, sep = sep))
  return(data)
}


#' create interaction column from factors
#' @export
#' @keywords internal
#' @examples
#' skylineconfig <- LFQServiceData::skylineconfig
#' skylineconfig$table$factorKeys()
#' skylineconfig$table$factorDepth <- 1
#' make_interaction_column_config(LFQServiceData::sample_analysis,skylineconfig)
make_interaction_column_config <- function(data, config, sep="."){
  columns <- config$table$fkeysDepth()
  data <- make_interaction_column(data, columns, sep = sep)
  return(data)
}

# Functions - Configuration ----
#' Helper function to extract all value slots in an R6 object
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
#' @export
#' @keywords internal
#' @examples
#'
#' skylineconfig <- createSkylineConfiguration(isotopeLabel = "Isotope.Label.Type",
#'  ident_qValue = "Detection.Q.Value")
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#' skylinePRMSampleData <- LFQServiceData::skylinePRMSampleData
#'
#' sample_analysis <- setup_analysis(skylinePRMSampleData, LFQServiceData::skylineconfig)
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

#' separates hierarchies into starting columns
#' @export
#' @keywords internal
separate_hierarchy <- function(data, config){
  for (hkey in config$table$hkeysDepth()) {
    data <- data %>% tidyr::separate( hkey, config$table$hierarchy[[hkey]], sep = config$sep, remove = FALSE)
  }
  return(data)
}

#' separates hierarchies into starting columns
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
#' @keywords internal
#' @examples
#'
#'  config <- LFQServiceData::skylineconfig
#'  config$table$isotopeLabel <- "Isotope.Label.Type"
#'  data <- LFQServiceData::sample_analysis
#'  xx <- complete_cases(data, config)
#'
complete_cases <- function(data, config) {
  message("completing cases")
  fkeys <- c(config$table$fileName,config$table$sampleName, config$table$factorKeys())
  hkeys <- c(config$table$isotopeLabel, config$table$hierarchyKeys())
  res <- tidyr::complete(
    data,
    tidyr::nesting(!!!syms(fkeys)),
    tidyr::nesting(!!!syms(hkeys))
  )
  return(res)
}


# Functions - Plotting ----
#' Plot peptide and fragments
plot_hierarchies_line_default <- function(data,
                                          proteinName,
                                          sample,
                                          intensity,
                                          peptide,
                                          fragment,
                                          factor,
                                          isotopeLabel,
                                          separate = FALSE,
                                          log_y = FALSE
) {
  if (length(isotopeLabel)) {
    if (separate) {
      formula <- paste(paste( isotopeLabel, collapse = "+"), "~", paste(factor , collapse = "+"))
      p <- ggplot(data, aes_string(x = sample,
                                   y = intensity,
                                   group = fragment,
                                   color = peptide
      ))
    }else{
      formula <- sprintf("~%s",paste(factor, collapse = " + "))
      data <- tidyr::unite(data, "fragment_label", fragment, isotopeLabel, remove = FALSE)
      p <- ggplot(data, aes_string(x = sample,
                                   y = intensity,
                                   group = "fragment_label",
                                   color = peptide
      ))
    }
    p <- p +  geom_point(aes_string(shape = isotopeLabel)) +
      geom_line(aes_string(linetype = isotopeLabel))
  }else{
    formula <- sprintf("~%s", paste(factor, collapse = " + "))
    p <- ggplot(data, aes_string(x = sample, y = intensity, group = fragment,  color = peptide))
    p <- p +  geom_point() + geom_line()
  }

  #p <- ggplot(data, aes_string(x = sample, y = intensity, group = fragment,  color= peptide, linetype = isotopeLabel))
  p <- p + facet_grid(as.formula(formula), scales = "free_x"   )
  p <- p + ggtitle(proteinName) + theme_classic()
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top")
  if (log_y) {
    p <- p + scale_y_continuous(trans = 'log10')
  }
  return(p)
}

#' extracts the relevant information from the configuration to make the plot.
#' @export
#' @keywords internal
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' conf <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#' xnested <- LFQServiceData::sample_analysis %>%
#'  group_by_at(conf$table$hkeysDepth()) %>% tidyr::nest()
#'
#' LFQService::plot_hierarchies_line(xnested$data[[1]], xnested$protein_Id[[1]],conf )
#'
plot_hierarchies_line <- function(res, proteinName,
                                  configuration,
                                  factor_level = 1,
                                  separate = FALSE){

  rev_hnames <- configuration$table$hierarchyKeys(TRUE)
  fragment <- rev_hnames[1]
  peptide <- rev_hnames[1]

  if(length(rev_hnames) > 2){
    peptide <- rev_hnames[2]
  }
  res <- LFQService:::plot_hierarchies_line_default(res,
                                                    proteinName = proteinName,
                                                    sample = configuration$table$sampleName,
                                                    intensity = configuration$table$getWorkIntensity(),
                                                    peptide = peptide,
                                                    fragment = fragment,
                                                    factor = configuration$table$factorKeys()[1:factor_level],
                                                    isotopeLabel = configuration$table$isotopeLabel,
                                                    separate = separate,
                                                    log_y = !configuration$parameter$is_intensity_transformed
  )
  return(res)
}


#' generates peptide level plots for all Proteins
#' @export
#' @keywords internal
#' @examples
#' library(tidyverse)
#' resDataStart <- LFQServiceData::testDataStart2954$resDataStart
#' config <-  LFQServiceData::testDataStart2954$config
#' res <- plot_hierarchies_line_df(resDataStart, config)
#' res[[1]]
#' config <- config$clone(deep = TRUE)
#' #TODO make it work for other hiearachy levels.
#' #config$table$hierarchyDepth = 2
#' #res <- plot_hierarchies_line_df(resDataStart, config)
#filteredPep <- resDataStart
plot_hierarchies_line_df <- function(filteredPep, config){
  factor_level <- config$table$factorDepth

  hierarchy_ID <- "hierarchy_ID"
  filteredPep <- filteredPep %>% tidyr::unite(hierarchy_ID , !!!syms(config$table$hkeysDepth()), remove = FALSE)

  xnested <- filteredPep %>% dplyr::group_by_at(hierarchy_ID) %>% tidyr::nest()

  figs <- xnested %>%
    dplyr::mutate(plot = map2(data, !!sym(hierarchy_ID) ,
                              plot_hierarchies_line,
                              factor_level = factor_level, config ))
  return(figs$plot)
}



#' add quantline to plot
#' @export
#' @keywords internal
#' @examples
#'
plot_hierarchies_add_quantline <- function(p, data, aes_y,  configuration){
  table <- configuration$table
  p + geom_line(data = data,
                aes_string(x = table$sampleName , y = aes_y, group = 1),
                size = 1.3,
                color = "black",
                linetype = "dashed") +
    geom_point(data = data,
               aes_string(x = table$sampleName , y = aes_y, group = 1), color = "black", shape = 10)
}

#' plot peptides by factors and its levels.
#'
#' @export
#' @keywords internal
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' conf <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#' xnested <- LFQServiceData::sample_analysis %>%
#'  group_by_at(conf$table$hkeysDepth()) %>% tidyr::nest()
#'
#' p <- plot_hierarchies_boxplot(xnested$data[[3]],
#'  xnested$protein_Id[[3]],
#'   conf ,
#'    hierarchy_level = NULL)
#' p
#' p <- plot_hierarchies_boxplot(xnested$data[[3]],
#'  xnested$protein_Id[[3]],
#'   conf )
#' p
#' p <- plot_hierarchies_boxplot(xnested$data[[3]],
#'  xnested$protein_Id[[3]],
#'   conf,
#'   hierarchy_level =  conf$table$hierarchyKeys()[3])
#' p
#' #ddd <- xnested$data[[3]]
#' #proteinName <- xnested$protein_Id[[3]]
#' #config <- conf
#' #boxplot = TRUE
#' #factor_level = 1
#' plot_hierarchies_boxplot(xnested$data[[3]], xnested$protein_Id[[3]],conf, beeswarm = FALSE )
plot_hierarchies_boxplot <- function(ddd,
                                     proteinName,
                                     config,
                                     hierarchy_level = tail(config$table$hierarchyKeys(),1),
                                     beeswarm = TRUE){

  isotopeLabel <- config$table$isotopeLabel
  ddd <- LFQService::make_interaction_column( ddd , config$table$fkeysDepth())
  p <- ggplot(ddd, aes_string(x = "interaction",
                              y = config$table$getWorkIntensity()
  )) + theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    ggtitle(proteinName)

  if (!config$parameter$is_intensity_transformed) {
    p <- p + scale_y_continuous(trans = "log10")
  }
  p <- p + geom_boxplot()

  if ( beeswarm ) {
    p <- p + ggbeeswarm::geom_quasirandom(dodge.width = 0.7, grouponX = FALSE)
  }
  if (!is.null( hierarchy_level ) && hierarchy_level %in% colnames(ddd)) {
    p <- p + facet_grid( formula(paste0("~", hierarchy_level ) ))
  }
  return(p)
}


#' generates peptide level plots for all Proteins
#' @export
#' @keywords internal
#' @examples
#' resDataStart <- LFQServiceData::testDataStart2954$resDataStart
#' config <-  LFQServiceData::testDataStart2954$config
#' res <- plot_hierarchies_boxplot_df(resDataStart, config)
#' res$boxplot[[1]]
#' res <- plot_hierarchies_boxplot_df(resDataStart,
#'  config,
#'  hierarchy_level = "peptide_Id")
#' res$boxplot[[1]]
#' config <- config$clone(deep = TRUE)
#' #TODO make it work for other hiearachy levels.
#' config$table$hierarchyDepth = 2
#' res <- plot_hierarchies_boxplot_df(resDataStart, config)
#' res$boxplot[[1]]
plot_hierarchies_boxplot_df <- function(filteredPep,
                                        config,
                                        hierarchy = "protein_Id",
                                        hierarchy_level = NULL){

  xnested <- filteredPep %>% dplyr::group_by_at(hierarchy) %>% tidyr::nest()

  figs <- xnested %>%
    dplyr::mutate(boxplot = map2(data, !!sym(hierarchy),
                                 plot_hierarchies_boxplot,
                                 config ,
                                 hierarchy_level = hierarchy_level))
  return(dplyr::select_at(figs, c(hierarchy, "boxplot")))
}

# Functions - summary ----

# Functions - summarize factors ----

#' table factors
#' @export
#' @keywords internal
#' @examples
#' library(tidyverse)
#' conf <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#' configuration <- conf
#' data <- LFQServiceData::sample_analysis
#' xx <- table_factors(data,configuration )
#' xx %>% dplyr::group_by(!!sym(configuration$table$factorKeys())) %>% dplyr::summarize(n = n())
#'
table_factors <- function(data, configuration){
  factorsTab <- data %>% dplyr::select(c(configuration$table$fileName, configuration$table$sampleName, configuration$table$factorKeys())) %>%
    dplyr::distinct() %>%
    arrange(!!sym(configuration$table$sampleName))
  return(factorsTab)
}


# Functions - summarize hierarchies

#' Count distinct elements for each level of hierarchy
#'
#' @export
#' @keywords internal
#' @examples
#' library(tidyverse)
#' library(LFQService)
#' skylineconfig <- createSkylineConfiguration(isotopeLabel = "Isotope.Label.Type", ident_qValue = "Detection.Q.Value")
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#'
#' sample_analysis <- setup_analysis(LFQServiceData::skylinePRMSampleData, skylineconfig)
#' hierarchy_counts(sample_analysis, skylineconfig)
hierarchy_counts <- function(x, configuration){
  hierarchy <- names( configuration$table$hierarchy )
  res <- x %>% dplyr::group_by_at(configuration$table$isotopeLabel) %>%
    dplyr::summarise_at( hierarchy, n_distinct )
  return(res)
}

#' Count distinct elements for each level of hierarchy per sample
#' @export
#' @keywords internal
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' skylineconfig <- createSkylineConfiguration(isotopeLabel = "Isotope.Label.Type", ident_qValue = "Detection.Q.Value")
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#' skylinePRMSampleData <- LFQServiceData::skylinePRMSampleData
#' sample_analysis <- setup_analysis(skylinePRMSampleData, skylineconfig)
#' res <- hierarchy_counts_sample(sample_analysis, skylineconfig)
#' res()
#' res("long")
#' res("plot")
hierarchy_counts_sample <- function(data,
                                    configuration)
{
  hierarchy <- configuration$table$hierarchyKeys()
  summary <- data %>% dplyr::filter(!is.na(!!sym(configuration$table$getWorkIntensity() ))) %>%
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



#' Summarize peptide Counts
#' @export
#' @keywords internal
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' skylineconfig <- createSkylineConfiguration(isotopeLabel = "Isotope.Label.Type",
#'  ident_qValue = "Detection.Q.Value")
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#' skylinePRMSampleData <- LFQServiceData::skylinePRMSampleData
#' configuration <- skylineconfig$clone(deep=TRUE)
#' sample_analysis <- setup_analysis(skylinePRMSampleData, configuration)
#'
#'
#' summarize_hierarchy(sample_analysis, configuration)
#' summarize_hierarchy(sample_analysis, configuration, factors = character())
#'
#' summarize_hierarchy(sample_analysis, configuration,
#'  hierarchy = configuration$table$hkeysDepth() )
#' summarize_hierarchy(sample_analysis, configuration,
#'  hierarchy = NULL, factors = configuration$table$fkeysDepth() )
#' configuration$table$hierarchyDepth = 1
#' summarize_hierarchy(sample_analysis, configuration,
#'  factors = configuration$table$fkeysDepth())
#' configuration$table$hierarchyDepth = 2
#' summarize_hierarchy(sample_analysis, configuration)
#' configuration$table$hierarchyDepth = 3
#' summarize_hierarchy(sample_analysis, configuration )
#' configuration$table$hierarchyDepth = 4
#' summarize_hierarchy(sample_analysis, configuration )
#' summarize_hierarchy(testDataStart2954$resDataStart, testDataStart2954$config)
summarize_hierarchy <- function(x,
                                configuration,
                                hierarchy = configuration$table$hkeysDepth(),
                                factors = character())
{
  all_hierarchy <- c(configuration$table$isotopeLabel, configuration$table$hierarchyKeys() )
  #factors <- configuration$table$factorKeys()[ifelse(factor_level < 1, 0, 1): factor_level]

  precursor <- x %>% dplyr::select(factors,all_hierarchy) %>% dplyr::distinct()
  x3 <- precursor %>% dplyr::group_by_at(c(factors,hierarchy)) %>%
    dplyr::summarize_at( setdiff(all_hierarchy,hierarchy),
                         list( n = dplyr::n_distinct))
  return(x3)
}

#' Light only version.
#' Summarize Protein counts
#' @export
#' @keywords internal
#' @importFrom dplyr group_by_at
#' @examples
#' library(LFQService)
#' skylineconfig <- createSkylineConfiguration(isotopeLabel="Isotope.Label.Type",
#'   ident_qValue="Detection.Q.Value")
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#' skylinePRMSampleData <- LFQServiceData::skylinePRMSampleData
#' sample_analysis <- setup_analysis(skylinePRMSampleData, skylineconfig)
#' LFQService:::summarizeProteins(sample_analysis, skylineconfig)
#'configuration <- skylineconfig$clone(deep=TRUE)
#'summarize_hierarchy(testDataStart2954$resDataStart, testDataStart2954$config)
#'summarizeProteins(testDataStart2954$resDataStart, testDataStart2954$config)
summarizeProteins <- function(x, configuration ){
  warning("DEPRECATED use summarize_hierarchy instead")
  rev_hierarchy <- configuration$table$hierarchyKeys(TRUE)

  precursorSum <- x %>% dplyr::select(rev_hierarchy) %>% dplyr::distinct() %>%
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
                     minNrFragments= min(minNrFragments),
                     maxNrFragments = max(maxNrFragments)
    )
  proteinPeptide <- proteinSum %>% tidyr::unite(Precursors ,minNrPrecursors , maxNrPrecursors, sep="-", remove=FALSE)
  proteinPeptide <- proteinPeptide %>% tidyr::unite(Fragments ,minNrFragments , maxNrFragments, sep="-", remove=FALSE)
  return(proteinPeptide)
}

# Functions - Missigness ----
#' compute missing statistics
#' @export
#' @keywords internal
#' @examples
#' library(tidyverse)
#' library(LFQService)
#'
#' skylineconfig <- LFQServiceData::skylineconfig$clone(deep=TRUE)
#' skylineconfig$parameter$qVal_individual_threshold <- 0.01
#' xx <- LFQService::removeLarge_Q_Values(LFQServiceData::sample_analysis,
#'    skylineconfig)
#' xx <- complete_cases(xx, skylineconfig)
#' interaction_missing_stats(xx, skylineconfig)$data %>% arrange(desc(nrNAs))
#' tmp <- interaction_missing_stats(xx, skylineconfig,
#'  factors= character(),
#'   hierarchy = skylineconfig$table$hierarchyKeys()[1])$data
#' tmp %>% mutate(perc = nrNAs/nrReplicates )
#' summarize_hierarchy(xx , skylineconfig)
#'
interaction_missing_stats <- function(x,
                                      config,
                                      factors = config$table$fkeysDepth(),
                                      hierarchy = config$table$hierarchyKeys(),
                                      workIntensity = config$table$getWorkIntensity())
{
  x <- complete_cases(x, config)
  table <- config$table
  missingPrec <- x %>% group_by_at(c(factors,
                                     hierarchy,
                                     table$isotopeLabel
  ))
  missingPrec <- missingPrec %>%
    dplyr::summarize(nrReplicates = n(),
                     nrNAs = sum(is.na(!!sym(workIntensity))),
                     meanArea = mean(!!sym(workIntensity), na.rm = TRUE)) %>%
    mutate(nrMeasured = nrReplicates-nrNAs) %>% dplyr::ungroup()
  return(list(data = missingPrec,
              summaries = c("nrReplicates","nrNAs","nrMeasured","meanArea")))
}
#' Compute interaction averages and
#' impute data using mean of lowest 0.1 (default)
#'
#' used in Acetylation project p2916
#' @export
#' @keywords internal
#' @return function
#' @examples
#'
#' skylineconfig$parameter$qVal_individual_threshold <- 0.01
#' xx <- LFQService::removeLarge_Q_Values(sample_analysis, skylineconfig)
#' xx <- complete_cases(xx, skylineconfig)
#' tmp <- interaction_missing_stats(xx, skylineconfig)
#' fun <- missigness_impute_interactions(xx, skylineconfig)
#' dd <- fun("long")
#' head(dd)
#' undebug(fun)
#' xxx <- (fun("nrReplicates"))
#' xxx <- fun("all")
#' head(xxx)
missigness_impute_interactions <- function(mdataTrans,
                                           config,
                                           factors = config$table$fkeysDepth(),
                                           probs = 0.1){
  x <- interaction_missing_stats(mdataTrans, config, factors = factors)
  x_summaries <- x$summaries
  xx <- x$data
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
        dplyr::select( -one_of(setdiff(x_summaries,"nrReplicates" ))) %>%
        tidyr::spread(interaction, nrReplicates, sep = ".nrReplicates.") %>%
        arrange(!!!syms(pid)) %>%
        dplyr::ungroup()
      nrMeasured <- xx%>% dplyr::select(-one_of(setdiff(x_summaries,"nrMeasured" ) )) %>%
        tidyr::spread(interaction, nrMeasured, sep = ".nrMeasured.") %>%
        arrange(!!!syms(pid)) %>% dplyr::ungroup()

      meanArea <- xx %>% dplyr::select(-one_of(setdiff(x_summaries,"meanArea" ) )) %>%
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
#' @export
#' @keywords internal
#' @examples
#' skylineconfig$parameter$qVal_individual_threshold <- 0.01
#' xx <- LFQService::removeLarge_Q_Values(sample_analysis, skylineconfig)
#' xx <- complete_cases(xx, skylineconfig)
#' fun <- missigness_impute_factors_interactions(xx, skylineconfig)
#' fun <- missigness_impute_factors_interactions(xx, skylineconfig, value = "imputed")
#' head(fun)
#' dim(fun)
#' dim(distinct(fun[,1:6]))
#'
missigness_impute_factors_interactions <-
  function(data,
           config,
           value = c("nrReplicates", "nrMeasured", "meanArea", "imputed"),
           probs = 0.1,
           add.prefix = FALSE){
    value <- match.arg(value)
    fac_fun <- list()
    fac_fun[["interaction"]] <- missigness_impute_interactions(data,
                                                               config,
                                                               probs = probs)
    if (config$table$factorDepth > 1 ) { # if 1 only then done
      for (factor in config$table$fkeysDepth()) {
        fac_fun[[factor]] <- missigness_impute_interactions(data,
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

  xx <- missigness_impute_factors_interactions(data, config, "imputed" )
  message("missigness_impute_factors_interactions : imputed")
  imputed <- missigness_impute_contrasts(xx, config, contrasts)
  xx <- missigness_impute_factors_interactions(data, config, "meanArea" )
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
      dd_long <- dplyr::filter(dd_long, ! contrast %in% names(contrasts))
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
#' @examples
#' library(tidyverse)
#' library(LFQService)
#' sample_analysis <- LFQServiceData::sample_analysis
#' skylineconfig <- LFQServiceData::skylineconfig
#'
#' xx <- complete_cases(sample_analysis, skylineconfig)
#' skylineconfig$parameter$qVal_individual_threshold <- 0.01
#' xx <- LFQService::removeLarge_Q_Values(sample_analysis, skylineconfig)
#' xx <- complete_cases(xx, skylineconfig)
#' missigness_histogram(xx, skylineconfig)
#'
#' missingPrec <- interaction_missing_stats(xx, skylineconfig)
#'
#' setNa <- function(x){ifelse(x < 100, NA, x)}
#' sample_analysis %>% dplyr::mutate(Area = setNa(Area)) -> sample_analysis
#' missigness_histogram(sample_analysis, skylineconfig)
#'
missigness_histogram <- function(x, config, showempty = TRUE, factors = config$table$fkeysDepth()){
  table <- config$table
  missingPrec <- interaction_missing_stats(x, config , factors)$data
  missingPrec <- missingPrec %>%  dplyr::ungroup() %>% dplyr::mutate(nrNAs = as.factor(nrNAs))

  if (showempty) {
    if (config$parameter$is_intensity_transformed) {
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

  if (!config$parameter$is_intensity_transformed) {
    p <- p + scale_x_log10()
  }
  p
}

#' cumulative sums of missing
#' @export
#' @keywords internal
#' @examples
#'
#' setNa <- function(x){ifelse(x < 100, NA, x)}
#' sample_analysis %>% dplyr::mutate(Area = setNa(Area)) -> sample_analysis
#' res <- missingness_per_condition_cumsum(sample_analysis,skylineconfig)
#' names(res)
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
#' @examples
#' setNa <- function(x){ifelse(x < 100, NA, x)}
#' sample_analysis %>% dplyr::mutate(Area = setNa(Area)) -> sample_analysis
#' res <- missingness_per_condition(sample_analysis,skylineconfig)
#' names(res)
#' res$data
#' res$figure
#' print(res$figure)
#' config <- skylineconfig$clone(deep = TRUE)
#' x <- sample_analysis
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
#' setNa <- function(x){ifelse(x < 100, NA, x)}
#' sample_analysis <- LFQServiceData::sample_analysis
#' skylineconfig <- LFQServiceData::skylineconfig
#' sample_analysis %>% dplyr::mutate(Area = setNa(Area)) -> sample_analysis
#' x<-spreadValueVarsIsotopeLabel(sample_analysis,skylineconfig)
#' head(x)
#'
#' x<-spreadValueVarsIsotopeLabel(sample_analysis_HL,skylineconfig_HL)
#' head(x[,5:ncol(x)])
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

.reestablish_condition <- function(data,
                                   medpolishRes,
                                   configuration
){
  table <- configuration$table
  xx <- data %>%  dplyr::select(c(table$sampleName,
                                  table$factorKeys(),
                                  table$fileName,
                                  table$isotopeLabel)) %>% dplyr::distinct()
  res <- dplyr::inner_join(xx,medpolishRes, by = table$sampleName)
  res
}


#' compute tukeys median polish from peptide or precursor intensities
#' @family matrix manipulation
#' @param name if TRUE returns the name of the summary column
#' @export
#' @keywords internal
#'
#' @examples
#' medpolishPly(name = T)
#' gg <- matrix(runif(20),4,5)
#' rownames(gg) <- make.names(1:4)
#' colnames(gg) <- make.names(1:5)
#' mx <- medpolishPly(gg)
#'
#' # compare it with other methods of protein inference
#' dd <- tidyr::gather(as_tibble(gg))
#' #x <- robust::lmRob(value ~ key, data = dd )
#' #pred_lmRob <- c(coef(x)[1] , coef(x)[1] + coef(x)[-1])
#' xl <- lm(value ~ key , data = dd)
#' pred_lm <- c(coef(xl)[1] , coef(xl)[1] + coef(xl)[-1])
#' xr <- MASS::rlm(value ~ key , data = dd)
#' pred_rlm <- c(coef(xr)[1] , coef(xr)[1] + coef(xr)[-1])
#'
#' xx <- cbind(medpolish = mx$medpolish, pred_lm = pred_lm,pred_rlm = pred_rlm )
#' head(xx)
#' matplot(xx, type = "l")
#'
medpolishPly <- function(x, name = FALSE){
  if (name) {
    return("medpolish")
  }
  X <- medpolish(x,na.rm = TRUE, trace.iter = FALSE, maxiter = 10);
  res <- tibble("sampleName" = names(X$col) , medpolish = X$col + X$overall)
  res
}


#' Summarizes the intensities within hierarchy
#'
#' @param func - a function working on a matrix of intensities for each protein.
#' @return retuns function object
#' @keywords internal
#' @export
#' @importFrom purrr map
#' @examples
#'
#' library(LFQService)
#' library(tidyverse)
#' config <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#' data <- LFQServiceData::sample_analysis
#' x <- intensity_summary_by_hkeys(data, config, func = medpolishPly)
#'
#' res <- x("unnest")
#' x("unnest")$data %>% dplyr::select(config$table$hierarchyKeys()[1] , "medpolish") %>% tidyr::unnest()
#' config <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#' config$table$hierarchyDepth <- 1
#' x <- intensity_summary_by_hkeys(data, config, func = medpolishPly)
#'
#' x("unnest")$data
#' xnested<-x()
#' dd <- x("plot")
#' dd$medpolishPly[[1]]
#' dd$plot[[2]]
#' # example how to add peptide count information
#'
#' tmp <- summarize_hierarchy(data, config)
#' tmp <- inner_join(tmp, x("wide")$data, by = config$table$hkeysDepth())
#' head(tmp)
intensity_summary_by_hkeys <- function(data, config, func)
{
  x <- as.list( match.call() )
  makeName <- make.names(as.character(x$func))
  config <- config$clone(deep = TRUE)

  xnested <- data %>% group_by_at(config$table$hkeysDepth()) %>% nest()

  xnested <- xnested %>%
    dplyr::mutate(spreadMatrix = map(data, extractIntensities, config))
  xnested <- xnested %>%
    dplyr::mutate(!!makeName := map(spreadMatrix, func))
  xnested <- xnested %>%
    dplyr::mutate(!!makeName := map2(data,!!sym(makeName),.reestablish_condition, config ))


  res_fun <- function(value = c("nested", "unnest", "wide", "plot"), DEBUG = FALSE){

    value <- match.arg(value)
    if (DEBUG) {
      return(list(config = config, value = value, xnested = xnested  ))
    }

    newconfig <- make_reduced_hierarchy_config(config,
                                               workIntensity = func(name = TRUE),
                                               hierarchy = config$table$hkeysDepth(names = FALSE))

    if (value == "nested") {
      return(list(xnested = xnested, config = newconfig))
    }else if (value == "unnest" || value == "wide") {
      unnested <- xnested %>%
        dplyr::select(config$table$hkeysDepth(), makeName) %>%
        tidyr::unnest(cols = c(medpolishPly)) %>%
        dplyr::ungroup()

      if (value == "wide") {
        wide <- LFQService::toWideConfig(unnested, newconfig)
        wide$config <- newconfig
        return(wide)
      }
      return(list(data = unnested, config = newconfig))
    }else if (value == "plot") {
      hierarchy_ID <- "hierarchy_ID"
      xnested <- xnested %>% tidyr::unite(hierarchy_ID , !!!syms(config$table$hkeysDepth()))
      figs <- xnested %>%
        dplyr::mutate(plot = map2(data, !!sym(hierarchy_ID) ,
                                  plot_hierarchies_line,
                                  factor_level = config$table$factorDepth, config ))

      figs <- figs %>%
        dplyr::mutate(plot = map2(plot, !!sym(makeName) ,
                                  plot_hierarchies_add_quantline, func(name = TRUE), config ))
      return(figs)
    }
  }
  return(res_fun)
}


#' median polish from normalized peptide intensities
#' @export
#' @keywords internal
#' @examples
#' library(tidyverse)
#' library(LFQService)
#' data <- LFQServiceData::dataIonstarNormalizedPep
#'
#' data$data <- data$data %>% filter(protein_Id %in% sample(protein_Id, 100))
#' res <- medpolish_protein_quants(data$data,
#' data$config )
#'
#' head(res("unnest")$data)
#'
medpolish_protein_quants <- function(data, config){
  protintensity <- LFQService::intensity_summary_by_hkeys(data ,
                                                          config,
                                                          medpolishPly)
  return(protintensity)
}

.string.to.colors <- function(string, colors = NULL){
  if (is.factor(string)) {
    string = as.character(string)
  }
  if (!is.null(colors)) {
    if (length(colors) != length(unique(string))) {
      stop("The number of colors must be equal to the number of unique elements.")
    }
    else {
      conv = cbind(unique(string), colors)
    }
  } else {
    conv = cbind(unique(string), rainbow(length(unique(string))))
  }
  unlist(lapply(string, FUN = function(x){conv[which(conv[,1] == x),2]}))
}

.number.to.colors = function(value, colors = c("red", "blue"), num = 100){
  cols = colorRampPalette(colors)(num)
  cols = 	cols[findInterval(value, vec = seq(from = min(value), to = max(value), length.out = num))]
  cols
}


#' plot correlation heatmap with annotations
#'
#' @export
#' @keywords internal
#' @importFrom pheatmap pheatmap
#' @examples
#' library(tidyverse)
#' data <- sample_analysis
#' config <- skylineconfig$clone( deep = TRUE )
#' LFQService::plot_heatmap_cor( data, config )
#' plot_heatmap_cor( data, config, R2 = TRUE )
#'
plot_heatmap_cor <- function(data,
                             config,
                             R2 = FALSE,
                             color = colorRampPalette(c("white", "red"))(1024),
                             ...){

  res <-  toWideConfig(data, config , as.matrix = TRUE)
  annot <- res$annotation
  res <- res$data


  cres <- cor(res, use = "pa")
  if (R2) {
    cres <- cres^2
  }

  factors <- dplyr::select_at(annot, config$table$factorKeys())
  factors <- as.data.frame(factors)
  rownames(factors) <- annot$sampleName

  res <- pheatmap::pheatmap(cres,
                            scale = "none",
                            silent = TRUE)

  res <- pheatmap::pheatmap(cres[res$tree_row$order,],
                            scale = "none",
                            cluster_rows  = FALSE,
                            annotation_col = factors,
                            show_rownames = F,
                            border_color = NA,
                            main = ifelse(R2, "R^2", "correlation"),
                            silent = TRUE,
                            color = color,
                            ... = ...)
  invisible(res)
}
#' plot heatmap with annotations
#'
#' @export
#' @keywords internal
#' @importFrom pheatmap pheatmap
#' @examples
#' library(tidyverse)
#' library(LFQService)
#' data <- LFQServiceData::sample_analysis
#' config <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#' graphics.off()
#' .Device
#'  p  <- plot_heatmap(data, config)
#' .Device
#'  print(p)
#'  .Device
#'  plot(1)
#'
#'  print(p)

plot_heatmap <- function(data, config, na_fraction = 0.4, ...){
  res <-  toWideConfig(data, config , as.matrix = TRUE)
  annot <- res$annotation
  resdata <- res$data

  factors <- dplyr::select_at(annot, config$table$factorKeys())
  factors <- as.data.frame(factors)
  rownames(factors) <- annot$sampleName

  resdata <- LFQService::removeNArows(resdata,floor(ncol(resdata)*na_fraction))


  # not showing row dendrogram trick
  res <- pheatmap::pheatmap(resdata,
                            scale = "row",
                            silent = TRUE)

  res <- pheatmap::pheatmap(resdata[res$tree_row$order,],
                            cluster_rows  = FALSE,
                            scale = "row",
                            annotation_col = factors,
                            show_rownames = F,
                            border_color = NA,
                            silent = TRUE,
                            ... = ...)

  invisible(res)
}

#' plot heatmap of NA values
#' @export
#' @keywords internal
#' @importFrom pheatmap pheatmap
#' @importFrom heatmap3 heatmap3 showLegend
#' @examples
#'
#' library(tidyverse)
#' library(LFQService)
#' data <- sample_analysis
#' config <- skylineconfig$clone(deep = TRUE)
#' tmp <- plot_NA_heatmap(data, config)
#' print(tmp)
#' xx <- plot_NA_heatmap(data, config,distance = "euclidean")
#' print(xx)
#' dev.off()
#' print(xx)
#' names(xx)
#'
plot_NA_heatmap <- function(data,
                            config,
                            limitrows = 10000,
                            distance = "binary"){
  res <-  toWideConfig(data, config , as.matrix = TRUE)
  annot <- res$annotation
  res <- res$data
  stopifnot(annot$sampleName == colnames(res))

  factors <- dplyr::select_at(annot, config$table$factorKeys())
  factors <- as.data.frame(factors)
  rownames(factors) <- annot$sampleName

  res[!is.na(res)] <- 0
  res[is.na(res)] <- 1
  allrows <- nrow(res)
  res <- res[apply(res,1, sum) > 0,]

  message("rows with NA's: ", nrow(res), "; all rows :", allrows, "\n")

  if (nrow(res) > 0) {
    res <- if (nrow(res) > limitrows ) {
      message("limiting nr of rows to:", limitrows,"\n")
      res[sample( 1:nrow(res),limitrows),]
    }else{
      res
    }

    # not showing row dendrogram trick
    resclust <- pheatmap::pheatmap(res,
                                   scale = "none",
                                   silent = TRUE,
                                   clustering_distance_cols = distance,
                                   clustering_distance_rows = distance)

    resclust <- pheatmap::pheatmap(res[resclust$tree_row$order,],
                                   cluster_rows  = FALSE,
                                   clustering_distance_cols = distance,
                                   scale = "none",
                                   annotation_col = factors,
                                   color = c("white","black"),
                                   show_rownames = FALSE,
                                   border_color = NA,
                                   legend = FALSE,
                                   silent = TRUE
    )
    return(resclust)
  } else {
    return(NULL)
  }
}



#' plot PCA
#' @export
#' @keywords internal
#' @import ggfortify
#' @examples
#'
#' library(tidyverse)
#' library(LFQService)
#' data <- sample_analysis
#' config <- skylineconfig$clone(deep = TRUE)
#'
#' tmp <- plot_pca(data, config, add_txt= TRUE)
#' print(tmp)
#' tmp <- plot_pca(data, config, add_txt= FALSE)
#' print(tmp)
#' plotly::ggplotly(tmp, tooltip = config$table$sampleName)
#'
plot_pca <- function(data , config, add_txt = FALSE, plotly = FALSE){
  wide <- toWideConfig(data, config ,as.matrix = TRUE)
  ff <- na.omit(wide$data)
  ff <- t(ff)
  xx <- as_tibble(prcomp(ff)$x, rownames = "sampleName")
  xx <- inner_join(wide$annotation, xx)


  sh <- config$table$fkeysDepth()[2]
  point <- (if (!is.na(sh)) {
    geom_point(aes(shape = !!sym(sh)))
  }else{
    geom_point()
  })

  text <- geom_text(aes(label = !!sym(config$table$sampleName)),check_overlap = TRUE,
                    nudge_x = 0.25,
                    nudge_y = 0.25 )

  x <- ggplot(xx, aes(x = PC1, y = PC2,
                      color = !!sym(config$table$fkeysDepth()[1]),
                      text = !!sym(config$table$sampleName))) +
    point + if (add_txt) {text}
  return(x)
}

