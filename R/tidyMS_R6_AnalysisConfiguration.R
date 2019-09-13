library(R6)

# AnalysisParameters ----
#' Analysis parameters
#' @export
AnalysisParameters <- R6::R6Class("AnalysisParameters",
                                  public = list(
                                    qVal_individual_threshold  = 0.05,
                                    qVal_experiment_threshold = 0.01,
                                    qVal_minNumber_below_experiment_threshold = 3,
                                    is_intensity_transformed = FALSE, # important for some plotting functions
                                    min_nr_of_notNA = 1, # how many values per transition total
                                    min_nr_of_notNA_condition = 0, # how many not missing in condition
                                    min_peptides_protein = 2
                                  )
)

# AnalysisTableAnnotation ----
#' Create Annotation
#' @export
AnalysisTableAnnotation <- R6Class("AnalysisTableAnnotation",
                                   public = list(

                                     fileName = NULL,
                                     factors = list(), # ordering is important - first is considered the main
                                     factorLevel=1,

                                     sampleName = "sampleName",
                                     # measurement levels
                                     hierarchy = list(),
                                     hierarchyLevel = 1,
                                     retentionTime = NULL,
                                     isotopeLabel = character(),
                                     # do you want to model charge sequence etc?

                                     ident_qValue = character(), # rename to score (smaller better)
                                     ident_Score = character(), # larger better
                                     workIntensity = NULL, # could be list with names and functions

                                     opt_rt  = NULL,
                                     opt_mz = NULL,



                                     initialize = function(){
                                     },
                                     getFactorLevel = function(){
                                       if(length(self$factorLevel) == 0){
                                         return(length(self$factors))
                                       }else{
                                         return(self$factorLevel)
                                       }
                                     },
                                     setWorkIntensity = function(colName){
                                       self$workIntensity <- c(self$workIntensity, colName)
                                     },
                                     getWorkIntensity = function(){
                                       return(tail(self$workIntensity, n=1))
                                     },
                                     popWorkIntensity=function(){
                                       res <- self$workIntensity[length(self$workIntensity)]
                                       self$workIntensity <- self$workIntensity[-length(self$workIntensity)]
                                       return(res)
                                     },
                                     idRequired = function(){
                                       "Id Columns which must be in the input data frame"
                                       idVars <- c(
                                         self$fileName,
                                         unlist(self$factors),
                                         unlist(self$hierarchy),
                                         self$isotopeLabel
                                       )
                                       return(idVars)
                                     },
                                     hierarchyKeys = function(rev = FALSE){
                                       if(rev){
                                         return(rev(names(self$hierarchy)))
                                       }else{
                                         return(names(self$hierarchy))
                                       }
                                     },
                                     hkeysLevel = function(names = TRUE){
                                       res <- head( self$hierarchy,n=self$hierarchyLevel)
                                       res <- if(names){
                                         names(res)
                                       }else{
                                         res
                                       }
                                       return(res)
                                     },
                                     factorKeys = function(){
                                       return(names(self$factors))
                                     },
                                     fkeysLevel = function(){
                                       res <- head(self$factors, n= self$factorLevel)
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
                                     }
                                   )
)
# AnalysisConfiguration ----
#' Analysis Configuration
#' @export
AnalysisConfiguration <- R6Class("AnalysisConfiguration",

                                 public = list(
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
#' @examples
#' make_reduced_hierarchy_config(skylineconfig, "testintensity", skylineconfig$table$hierarchy[1:2])
#'
make_reduced_hierarchy_config <- function(config, workIntensity , hierarchy ){
  newConfig <-config$clone(deep=TRUE)
  newConfig$table$hierarchy = hierarchy
  newConfig$table$workIntensity = workIntensity
  return(newConfig)
}

#' create interaction column from factors
#' @export
#' @examples
#' xx <- data.frame(A = c("a","a","a"), B = c("d","d","e"))
#' make_interaction_column(xx, c("B","A"))
#' make_interaction_column(xx, c("A"))
make_interaction_column <- function(data, columns, sep="."){
  intr <- dplyr::select(data, columns)
  intr <- purrr::map2_dfc(columns, intr, paste0)
  colnames(intr) <- paste0("interaction_",columns)
  colname <- "interaction"
  data <- data %>% dplyr::mutate(!!colname := interaction(intr, sep=sep))
  return(data)
}


#' create interaction column from factors
#' @export
#' @examples
#'
#' skylineconfig$table$factorKeys()
#' skylineconfig$table$factorLevel <- 1
#' make_interaction_column_config(LFQService::sample_analysis,skylineconfig)
make_interaction_column_config <- function(data, config, sep="."){
  columns <- config$table$fkeysLevel()
  data <- make_interaction_column(data, columns, sep=sep)
  return(data)
}

# Functions - Configuration ----
#' Helper function to extract all value slots in an R6 object
#' @param r6class r6 class
#' @export
R6extractValues <- function(r6class){
  tmp <- sapply(r6class, class)
  slots <- tmp[! tmp %in% c("environment", "function")]
  res <- list()
  for(i in names(slots)){
    if("R6" %in% class(r6class[[i]])){
      res[[i]]  <- R6extractValues(r6class[[i]])
    }else{
      res[[i]] <- r6class[[i]]
    }
  }
  return(res)
}

#' Extracts columns relevant for a configuration from a data frame
#' @export
#' @examples
#'
#' skylineconfig <- createSkylineConfiguration(isotopeLabel="Isotope.Label.Type",
#'  ident_qValue="Detection.Q.Value")
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#' data(skylinePRMSampleData)
#'
#' sample_analysis <- setup_analysis(skylinePRMSampleData, skylineconfig)
#'
setup_analysis <- function(data, configuration ){
  table <- configuration$table
  for(i in 1:length(table$hierarchy))
  {
    data <- tidyr::unite(data, UQ(sym(table$hierarchyKeys()[i])), table$hierarchy[[i]],remove = FALSE, sep=configuration$sep)
  }
  data <- dplyr::select(data , -dplyr::one_of(dplyr::setdiff(unlist(table$hierarchy), table$hierarchyKeys() )))

  for(i in 1:length(table$factors))
  {
    if( length(table$factors[[i]]) > 1){
      data <- tidyr::unite(data, UQ(sym(table$factorKeys()[i])), table$factors[[i]],remove = FALSE, sep=configuration$sep)
    }else{
      newname <-table$factorKeys()[i]
      data <- dplyr::mutate(data, !!newname := !!sym(table$factors[[i]]))
    }
  }

  sampleName <- table$sampleName

  if(!sampleName  %in% names(data)){
    message("creating sampleName")

    data <- data %>%  tidyr::unite( UQ(sym( sampleName)) , unique(unlist(table$factors)), remove = TRUE , sep=configuration$sep) %>%
      dplyr::select(sampleName, table$fileName) %>% distinct() %>%
      dplyr::mutate_at(sampleName, function(x){ x<- make.unique( x, sep=configuration$sep )}) %>%
      dplyr::inner_join(data, by=table$fileName)
  } else{
    warning(sampleName, " already exists")
  }

  data <- data %>% dplyr::select(-dplyr::one_of(dplyr::setdiff(unlist(table$factors), table$factorKeys())))

  # Make implicit NA's explicit
  data <- data %>% dplyr::select(c(configuration$table$idVars(),configuration$table$valueVars()))
  data <- completeCases( data , configuration)

  return( data )
}

#' separates hierarchies into starting columns
#' @export
separate_hierarchy <- function(data, config){
  for(hkey in config$table$hkeysLevel()){
    data <- data %>% tidyr::separate( hkey, config$table$hierarchy[[hkey]], sep=config$sep, remove=FALSE)
  }
  return(data)
}

#' separates hierarchies into starting columns
#' @export
separate_factors <- function(data, config){
  for(fkey in config$table$factorKeys()){
    data<- data %>% tidyr::separate( fkey, config$table$factors[[fkey]], sep=config$sep, remove=FALSE)
  }
  return(data)
}



#' Complete cases
#' @export
#' @examples
#'
#'  config <- skylineconfig
#'  config$table$isotopeLabel <- "Isotope.Label.Type"
#'  data <- LFQService::sample_analysis
#'  data
#'  xx <- completeCases(sample_analysis,skylineconfig)
completeCases <- function(data, config){
  data <- tidyr::complete(
    data ,
    tidyr::nesting(!!!syms(c(config$table$hierarchyKeys(), config$table$isotopeLabel))),
    tidyr::nesting(!!!syms(c(config$table$fileName , config$table$sampleName, config$table$factorKeys() )))
  )
  return(data)
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
                                          log_y=FALSE
){
  if(length(isotopeLabel)){
    if(separate){
      formula <- paste(paste( isotopeLabel, collapse="+"), "~", paste(factor , collapse = "+"))
      p <- ggplot(data, aes_string(x = sample,
                                   y = intensity,
                                   group=fragment,
                                   color= peptide
      ))
    }else{
      formula <- sprintf("~%s",paste(factor, collapse=" + "))
      data <- tidyr::unite(data, "fragment_label", fragment, isotopeLabel, remove = FALSE)
      p <- ggplot(data, aes_string(x = sample,
                                   y = intensity,
                                   group="fragment_label",
                                   color= peptide
      ))
    }
    p <- p +  geom_point(aes_string(shape= isotopeLabel)) + geom_line(aes_string(linetype = isotopeLabel))
  }else{
    formula <- sprintf("~%s", paste(factor, collapse=" + "))
    p <- ggplot(data, aes_string(x = sample, y = intensity, group=fragment,  color= peptide))
    p <- p +  geom_point() + geom_line()
  }

  #p <- ggplot(data, aes_string(x = sample, y = intensity, group=fragment,  color= peptide, linetype = isotopeLabel))
  p <- p + facet_grid(as.formula(formula), scales = "free_x"   )
  p <- p + ggtitle(proteinName) + theme_classic()
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top")
  if(log_y){
    p <- p + scale_y_continuous(trans='log10')
  }
  return(p)
}

#' extracts the relevant information from the configuration to make the plot.
#' @export
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' conf <- LFQService::skylineconfig$clone(deep=TRUE)
#' xnested <- LFQService::sample_analysis %>%
#'  group_by_at(conf$table$hkeysLevel()) %>% tidyr::nest()
#'
#' LFQService::plot_hierarchies_line(xnested$data[[1]], xnested$protein_Id[[1]],conf )
#'
plot_hierarchies_line <- function(res, proteinName,
                                  configuration,
                                  factor_level=1,
                                  separate=FALSE){

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
#'
#' @examples
#' resDataStart <- LFQService::testDataStart2954$resDataStart
#' config <-  LFQService::testDataStart2954$config
#' res <- plot_hierarchies_line_df(resDataStart, config)
#' res[[1]]
#' config <- config$clone(deep=TRUE)
#' #TODO make it work for other hiearachy levels.
#' #config$table$hierarchyLevel=2
#' #res <- plot_hierarchies_line_df(resDataStart, config)
#filteredPep <- resDataStart
plot_hierarchies_line_df <- function(filteredPep, config){
  factor_level <- config$table$factorLevel

  hierarchy_ID <- "hierarchy_ID"
  filteredPep <- filteredPep %>% tidyr::unite(hierarchy_ID , !!!syms(config$table$hkeysLevel()))

  xnested <- filteredPep %>% dplyr::group_by_at(hierarchy_ID) %>% tidyr::nest()

  figs <- xnested %>%
    dplyr::mutate(plot = map2(data, !!sym(hierarchy_ID) ,
                              plot_hierarchies_line,
                              factor_level = factor_level, config ))
  return(figs$plot)
}



#' add quantline to plot
#' @export
#' @examples
#'
plot_hierarchies_add_quantline <- function(p, data, aes_y,  configuration){
  table <- configuration$table
  p + geom_line(data=data,
                aes_string(x = table$sampleName , y = aes_y, group=1),
                size=1.3,
                color="black",
                linetype="dashed") +
    geom_point(data=data,
               aes_string(x = table$sampleName , y = aes_y, group=1), color="black", shape=10)
}

#' plot peptides by factors and its levels.
#'
#' @export
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' conf <- LFQService::skylineconfig$clone(deep=TRUE)
#' xnested <- LFQService::sample_analysis %>%
#'  group_by_at(conf$table$hkeysLevel()) %>% tidyr::nest()
#'
#' plot_hierarchies_boxplot(xnested$data[[3]], xnested$protein_Id[[3]],conf )
#' #ddd <- xnested$data[[3]]
#' #proteinName <- xnested$protein_Id[[3]]
#' #config <- conf
#' #boxplot=TRUE
#' #factor_level = 1
#' plot_hierarchies_boxplot(xnested$data[[3]], xnested$protein_Id[[3]],conf, boxplot=FALSE )
plot_hierarchies_boxplot <- function(ddd, proteinName, config , factor_level = 1, boxplot=TRUE){

  stopifnot(factor_level <= length(config$table$factorKeys()))
  isotopeLabel <- config$table$isotopeLabel
  ddd %>%  tidyr::gather("factor","level",config$table$factorKeys()[1:factor_level]) -> ddlong
  ddlong <- as.data.frame(unclass(ddlong))
  p <- ggplot(ddlong, aes_string(x = tail(config$table$hierarchyKeys(),1),
                                 y = config$table$getWorkIntensity(),
                                 fill="level"
  )) +
    facet_wrap(~factor, ncol=1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle =45, hjust = 1, vjust=1)) +
    ggtitle(proteinName)
  if(!config$parameter$is_intensity_transformed){
    p <- p + scale_y_continuous(trans="log10")
  }
  if(boxplot){
    p <- p + geom_boxplot()
  }else{
    p <- p + ggbeeswarm::geom_quasirandom(dodge.width=0.7)
  }
  return(p)
}
#' generates peptide level plots for all Proteins
#' @export
#'
#' @examples
#' resDataStart <- LFQService::testDataStart2954$resDataStart
#' config <-  LFQService::testDataStart2954$config
#' res <- plot_hierarchies_boxplot_df(resDataStart, config)
#' res[[1]]
#' config <- config$clone(deep=TRUE)
#' #TODO make it work for other hiearachy levels.
#' config$table$hierarchyLevel=2
#' res <- plot_hierarchies_boxplot_df(resDataStart, config)
plot_hierarchies_boxplot_df <- function(filteredPep, config){
  factor_level <- config$table$factorLevel

  hierarchy_ID <- "hierarchy_ID"
  filteredPep <- filteredPep %>% tidyr::unite(hierarchy_ID , !!!syms(config$table$hkeysLevel()))

  xnested <- filteredPep %>% dplyr::group_by_at(hierarchy_ID) %>% tidyr::nest()

  figs <- xnested %>%
    dplyr::mutate(boxplot = map2(data, !!sym(hierarchy_ID),
                                 plot_hierarchies_boxplot,
                                 config ,
                                 factor_level = factor_level))
  return(figs$boxplot)
}

# Functions - summary ----

# Functions - summarize factors ----

#' table factors
#' @export
#' @examples
#' library(tidyverse)
#' conf <- LFQService::skylineconfig$clone(deep=TRUE)
#' configuration <- conf
#' data <- LFQService::sample_analysis
#' xx <- table_factors(data,configuration )
#' xx %>% dplyr::group_by(!!sym(configuration$table$factorKeys())) %>% dplyr::summarize(n = n())
#'
table_factors <- function(data, configuration){
  factorsTab <- data %>% dplyr::select(c(configuration$table$fileName, configuration$table$sampleName, configuration$table$factorKeys())) %>%
    distinct() %>%
    arrange(!!sym(configuration$table$sampleName))
  return(factorsTab)
}


# Functions - summarize hierarchies

#' Count distinct elements for each level of hierarchy
#'
#' @export
#' @examples
#' library(tidyverse)
#' library(LFQService)
#' skylineconfig <- createSkylineConfiguration(isotopeLabel="Isotope.Label.Type", ident_qValue="Detection.Q.Value")
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#'
#' sample_analysis <- setup_analysis(skylinePRMSampleData, skylineconfig)
#' hierarchyCounts(sample_analysis, skylineconfig)
hierarchyCounts <- function(x, configuration){
  hierarchy <- names( configuration$table$hierarchy )
  res <- x %>% dplyr::group_by_at(configuration$table$isotopeLabel) %>%
    dplyr::summarise_at( hierarchy, n_distinct )
  return(res)
}

#' Count distinct elements for each level of hierarchy per sample
#' @export
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' skylineconfig <- createSkylineConfiguration(isotopeLabel="Isotope.Label.Type", ident_qValue="Detection.Q.Value")
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#' data(skylinePRMSampleData)
#' sample_analysis <- setup_analysis(skylinePRMSampleData, skylineconfig)
#' res <- hierarchy_counts_sample(sample_analysis, skylineconfig)
#' res()
#' res("long")
#' res("plot")
hierarchy_counts_sample <- function(data,
                                    configuration)
{
  hierarchy <- configuration$table$hierarchyKeys()
  summary <-data %>% dplyr::filter(! is.na(!!sym(configuration$table$getWorkIntensity() ))) %>%
    dplyr::group_by_at(c(configuration$table$isotopeLabel, configuration$table$sampleName)) %>%
    dplyr::summarise_at( hierarchy, n_distinct )

  res <- function(value = c("wide", "long", "plot")){
    value <- match.arg(value)
    if(value == "wide"){
      return(summary)
    }else{
      long <- summary %>% tidyr::gather("key",
                                        "nr",-dplyr::one_of(configuration$table$isotopeLabel,
                                                            configuration$table$sampleName))
      if(value == "long"){
        return(long)
      }else if(value == "plot"){
        nudgeval <- max(long$nr) * 0.05
        ggplot(long, aes(x = sampleName, y = nr)) +
          geom_bar(stat = "identity", position = "dodge", colour="black", fill="white") +
          facet_wrap( ~ key, scales = "free_y", ncol = 1) +
          geom_text(aes(label = nr), nudge_y = nudgeval) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
      }
    }
  }
  return(res)
}




#' Light only version.
#' Summarize Protein counts
#' @export
#' @importFrom dplyr group_by_at
#' @examples
#' skylineconfig <- createSkylineConfiguration(isotopeLabel="Isotope.Label.Type",
#'   ident_qValue="Detection.Q.Value")
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#' data(skylinePRMSampleData)
#' sample_analysis <- setup_analysis(skylinePRMSampleData, skylineconfig)
#' summarizeProteins(sample_analysis, skylineconfig)
#'
summarizeProteins <- function( x, configuration ){
  rev_hierarchy <- configuration$table$hierarchyKeys(TRUE)

  precursorSum <- x %>% dplyr::select(rev_hierarchy) %>% distinct() %>%
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
#' Summarize peptide Counts
#' @export
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' skylineconfig <- createSkylineConfiguration(isotopeLabel="Isotope.Label.Type",
#'  ident_qValue="Detection.Q.Value")
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#' data(skylinePRMSampleData)
#' configuration <- skylineconfig
#' sample_analysis <- setup_analysis(skylinePRMSampleData, skylineconfig)
#'
#' #x <- sample_analysis
#' #configuration <- skylineconfig
#' #hierarchy = configuration$table$hkeysLevel()
#'
#' summarizeHierarchy(sample_analysis, skylineconfig)
#' summarizeHierarchy(sample_analysis, skylineconfig, factors=character())
#'
#' summarizeHierarchy(sample_analysis, skylineconfig,
#'  hierarchy = skylineconfig$table$hkeysLevel() )
#' summarizeHierarchy(sample_analysis, skylineconfig,
#'  hierarchy = NULL, factors=skylineconfig$table$fkeysLevel() )
#' skylineconfig$table$hierarchyLevel=1
#' summarizeHierarchy(sample_analysis, skylineconfig,
#'  factors = skylineconfig$table$fkeysLevel())
#' skylineconfig$table$hierarchyLevel=2
#' summarizeHierarchy(sample_analysis, skylineconfig)
#' skylineconfig$table$hierarchyLevel=3
#' summarizeHierarchy(sample_analysis, skylineconfig )
#' skylineconfig$table$hierarchyLevel=4
#' summarizeHierarchy(sample_analysis, skylineconfig )
#' summarizeHierarchy(testDataStart2954$resDataStart, testDataStart2954$config)
summarizeHierarchy <- function(x,
                               configuration,
                               hierarchy = configuration$table$hkeysLevel(),
                               factors=character())
{
  all_hierarchy <- c(configuration$table$isotopeLabel, configuration$table$hierarchyKeys() )
  #factors <- configuration$table$factorKeys()[ifelse(factor_level < 1, 0, 1): factor_level]

  precursor <- x %>% dplyr::select(factors,all_hierarchy) %>% dplyr::distinct()
  x3 <- precursor %>% dplyr::group_by_at(c(factors,hierarchy)) %>%
    dplyr::summarize_at( setdiff(all_hierarchy,hierarchy),
                         list( n = dplyr::n_distinct))
  return(x3)
}
# Functions - Missigness ----
#' compute missing statistics
#' @export
#' @examples
#' library(tidyverse)
#' library(LFQService)
#'
#'
#' skylineconfig$parameter$qVal_individual_threshold <- 0.01
#' xx <- LFQService::removeLarge_Q_Values(sample_analysis, skylineconfig)
#' xx <- completeCases(xx, skylineconfig)
#' interaction_missing_stats(xx, skylineconfig)
interaction_missing_stats <- function(x,
                                      configuration,
                                      factors = configuration$table$fkeysLevel(),
                                      workIntensity = configuration$table$getWorkIntensity())
{

  x <- completeCases(x, configuration)
  table <- configuration$table
  missingPrec <- x %>% group_by_at(c(factors,
                                     table$hierarchyKeys(),
                                     table$isotopeLabel
  ))
  missingPrec <- missingPrec %>%
    dplyr::summarize(nrReplicates = n(),
                     nrNAs = sum(is.na(!!sym(workIntensity))),
                     meanArea = mean(!!sym(workIntensity), na.rm=TRUE)) %>%
    arrange(desc(nrNAs)) %>% dplyr::ungroup()
  missingPrec
}
#' @export
#'
getMissingStats <- function(x,
                            configuration,
                            factors = configuration$table$fkeysLevel(),
                            workIntensity = configuration$table$getWorkIntensity()){
  message("DEPRECATED, use missing_stats")
  interaction_missing_stats(x, configuration, factors = factors, workIntensity = workIntensity)
}


#' Compute interaction averages and impute data using mean of lowest 10%
#'
#' used in Acetylation project p2916
#' @export
#' @return function
#'
missigness_impute_interactions <- function(mdataTrans,
                                           pepConfig,
                                           factors = pepConfig$table$fkeysLevel(),
                                           probs = 0.1){
  xx <- interaction_missing_stats(mdataTrans, pepConfig, factors=factors)
  xx <- make_interaction_column(xx, factors, sep=":")

  #xx %>% mutate(perc_missing= nrNAs/nrReplicates*100) -> xx
  xx %>% mutate(nrMeasured = nrReplicates - nrNAs) -> xx

  lowerMean <- function(meanArea, probs = probs){
    meanAreaNotNA <- na.omit(meanArea)
    small10 <- meanAreaNotNA[meanAreaNotNA < quantile(meanAreaNotNA, probs= probs)]
    meanArea[is.na(meanArea)] <- mean(small10)
    return(meanArea)
  }

  xx <- xx %>% group_by(interaction) %>% mutate(imputed = lowerMean(meanArea,probs=0.2))

  res <- function(value = c("long", "nrReplicates", "nrMeasured", "meanArea", "imputed", "allWide", "all" ), add.prefix = TRUE){
    value <- match.arg(value)
    if(value == "long"){
      return(xx)
    }else{
      xx %>% dplyr::select(-one_of(c("nrNAs", factors))) -> xx


      nrReplicates <- xx %>% dplyr::select( -meanArea, -nrMeasured, -imputed) %>%
        tidyr::spread(interaction, nrReplicates, sep=".nrReplicates.") %>%
        arrange(protein_Id) %>% ungroup()

      nrMeasured <- xx%>% dplyr::select( -meanArea, -nrReplicates, -imputed) %>%
        tidyr::spread(interaction, nrMeasured, sep=".nrMeasured.") %>%
        arrange(protein_Id) %>% ungroup()

      meanArea <- xx%>% dplyr::select(-nrReplicates, -nrMeasured, -imputed) %>%
        tidyr::spread(interaction, meanArea, sep=".meanArea.") %>%
        arrange(protein_Id) %>% ungroup()

      meanAreaImputed <- xx%>% dplyr::select(-nrReplicates, -nrMeasured, -meanArea) %>%
        tidyr::spread(interaction, imputed, sep=".imputed.") %>%
        arrange(protein_Id) %>% ungroup()


      allTables <- list(meanArea= meanArea,
                        nrMeasured = nrMeasured,
                        nrReplicates = nrReplicates,
                        meanAreaImputed = meanAreaImputed)
      if(value == "all"){
        allTables[["long"]] <- xx
        return(allTables)
      }else if(value == "allWide"){
        return(purrr::reduce(allTables,inner_join))
      }else if(value == "nrReplicates"){
        srepl <- if(add.prefix){"nrRep."}else{""}
        colnames(nrReplicates) <- gsub("interaction.nrReplicates.", srepl ,colnames(nrReplicates))
        nrReplicates <- tibble::add_column( nrReplicates, "value" = value, .before = 1)

        return(nrReplicates)
      }else if(value == "nrMeasured"){
        srepl <- if(add.prefix){"nrMeas."}else{""}
        colnames(nrMeasured) <- gsub("interaction.nrMeasured.", srepl ,colnames(nrMeasured))
        nrMeasured <- tibble::add_column( nrMeasured, "value" = value, .before = 1)
        return(nrMeasured)
      }else if(value == "meanArea"){
        srepl <- if(add.prefix){"mean."}else{""}
        colnames(meanArea) <- gsub("interaction.meanArea.", srepl ,colnames(meanArea))
        meanArea <- tibble::add_column( meanArea, "value" = value, .before = 1)
        return(meanArea)
      }else if(value == "imputed"){
        srepl <- if(add.prefix){"mean.imp."}else{""}
        colnames(meanAreaImputed) <- gsub("interaction.imputed.", srepl ,colnames(meanAreaImputed))
        meanAreaImputed<- tibble::add_column( meanAreaImputed, "value" = value, .before = 1)
        return(meanAreaImputed)
      }
    }
  }

  #  nrMeasured %>% dplyr::select(starts_with("interaction")) -> nrMeasuredM
  #  nrReplicates %>% dplyr::select(starts_with("interaction")) -> nrReplicatesM
  return(res)
}


#' compute per group averages and impute values
#' should generalize at some stage
#' @export
missigness_impute_factors_interactions <- function(data,
                                                   config,
                                                   value = c("nrReplicates", "nrMeasured", "meanArea", "imputed"),
                                                   probs=0.1,
                                                   add.prefix=FALSE){
  value <- match.arg(value)
  fac_fun <- list()
  fac_fun[["interaction"]] <- missigness_impute_interactions(data, config,probs = probs)
  if(config$table$factorLevel > 1 ){ # if 1 only then done
    for(factor in config$table$fkeysLevel()){
      fac_fun[[factor]] <- missigness_impute_interactions(data, config,factors = factor,probs = probs)
    }
  }
  fac_res <- list()
  for(fun_name in names(fac_fun)){
    fac_res[[fun_name]] <- fac_fun[[fun_name]](value, add.prefix=add.prefix)
  }
  intfact <- purrr::reduce(fac_res,dplyr::inner_join, by = c(config$table$hkeysLevel(), config$table$isotopeLabel, "value"))
  return(intfact)
}

#' Compute fold changes given Contrasts
#'
#' @export
missigness_impute_contrasts <- function(data,
                                        config,
                                        Contrasts,
                                        agg_fun=function(x){median(x, na.rm = TRUE)})
{
  for(i in 1:length(Contrasts)){
    message(names(Contrasts)[i], "=", Contrasts[i],"\n")
    data <- dplyr::mutate(data, !!names(Contrasts)[i] := !!rlang::parse_expr(Contrasts[i]))
  }

  if(!is.null(agg_fun)){
    data <- data %>% group_by_at(c("value" , config$table$hkeysLevel())) %>%
      summarise_if(is.numeric, agg_fun)
  }
  return(data)
}

#' Compute fold changes given Contrasts 2
#'
#' @export
workflow_missigness_impute_contrasts <- function(data,
                                                 config,
                                                 Contrasts){

  xx <- missigness_impute_factors_interactions(data, config, "imputed" )
  imputed <- missigness_impute_contrasts(xx, config, Contrasts)
  xx <- missigness_impute_factors_interactions(data,config,"meanArea" )
  mean <- missigness_impute_contrasts(xx, config, Contrasts)
  dd <- bind_rows(imputed, mean)
  dd_long <- dd %>% gather("contrast","int_val",colnames(dd)[sapply(dd, is.numeric)])


  res <- function(value = c("long", "wide","raw"), what = c("contrasts", "factors", "all")){
    value <- match.arg( value )
    what  <- match.arg( what  )

    if(what == "contrasts"){
      dd_long <- dplyr::filter(dd_long, contrast %in% names(Contrasts))
    }else if(what == "factors"){
      dd_long <- dplyr::filter(dd_long, ! contrast %in% names(Contrasts))
    }else if(what == "all"){

    }

    if(value == "long"){
      long_xxxx <- dd_long %>% spread(value, int_val)
      return(long_xxxx)
    }else if(value == "wide"){
      dd <- dd_long %>% unite(contrast.v , value, contrast, sep="~") %>% spread(contrast.v, int_val)
      xxx_imputed <- inner_join(LFQService::summarizeHierarchy(data,config),dd)
      return(xxx_imputed)
    }else if(value == "raw"){
      return(dd_long)
    }
  }
}
#' Histogram summarizing missigness
#' @export
#' @examples
#' library(tidyverse)
#' library(LFQService)
#' xx <- completeCases(sample_analysis,skylineconfig)
#' skylineconfig$parameter$qVal_individual_threshold <- 0.01
#' xx <- LFQService::removeLarge_Q_Values(sample_analysis, skylineconfig)
#' xx <- completeCases(xx, skylineconfig)
#' missignessHistogram(xx, skylineconfig)
#'
#' missingPrec <- interaction_missing_stats(xx, skylineconfig)
#'
#' setNa <- function(x){ifelse(x < 100, NA, x)}
#' sample_analysis %>% dplyr::mutate(Area = setNa(Area)) -> sample_analysis
#' missignessHistogram(sample_analysis, skylineconfig)
#'
missignessHistogram <- function(x, configuration, showempty = TRUE, factors = configuration$table$fkeysLevel()){
  table <- configuration$table
  missingPrec <- interaction_missing_stats(x, configuration , factors)
  missingPrec <- missingPrec %>%  dplyr::ungroup() %>% dplyr::mutate(nrNAs = as.factor(nrNAs))

  if(showempty){
    if(configuration$parameter$is_intensity_transformed)
    {
      missingPrec <- missingPrec %>% dplyr::mutate(meanArea = ifelse(is.na(meanArea),1,meanArea))
    }else{
      missingPrec <- missingPrec %>% dplyr::mutate(meanArea = ifelse(is.na(meanArea),-20,meanArea))
    }

  }

  factors <- table$fkeysLevel()
  formula <- paste(table$isotopeLabel, "~", paste(factors, collapse = "+"))
  message(formula)

  p <- ggplot(missingPrec, aes(x = meanArea, fill = nrNAs, colour = nrNAs)) +
    geom_histogram(alpha = 0.2, position = "identity") +
    facet_grid(as.formula(formula)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  if(!configuration$parameter$is_intensity_transformed)
  {
    p <- p + scale_x_log10()
  }
  p
}

#' cumulative sums of missing
#' @export
#' @examples
#'
#' setNa <- function(x){ifelse(x < 100, NA, x)}
#' sample_analysis %>% dplyr::mutate(Area = setNa(Area)) -> sample_analysis
#' res <- missingPerConditionCumsum(sample_analysis,skylineconfig)
#' names(res)
#' print(res$figure)
missingPerConditionCumsum <- function(x,configuration,factors = configuration$table$fkeysLevel()){
  table <- configuration$table
  missingPrec <- interaction_missing_stats(x, configuration,factors)

  xx <-missingPrec %>% group_by_at(c(table$isotopeLabel, factors,"nrNAs","nrReplicates")) %>%
    dplyr::summarize(nrTransitions =n())

  xxcs <-xx %>% group_by_at( c(table$isotopeLabel,factors)) %>% arrange(nrNAs) %>%
    dplyr::mutate(cumulative_sum = cumsum(nrTransitions))
  res <- xxcs  %>% dplyr::select(-nrTransitions)

  formula <- paste(table$isotopeLabel, "~", paste(factors, collapse = "+"))
  message(formula)

  nudgeval = max(res$cumulative_sum) * 0.05
  p <- ggplot(res, aes(x= nrNAs, y = cumulative_sum)) +
    geom_bar(stat="identity", color="black", fill="white") +
    geom_text(aes(label = cumulative_sum), nudge_y = nudgeval) +
    facet_grid(as.formula(formula))

  res <- res %>% tidyr::spread("nrNAs","cumulative_sum")
  return(list(data =res, figure=p))
}

#' Summarize missing in condtion as barplot
#' @export
#' @examples
#' setNa <- function(x){ifelse(x < 100, NA, x)}
#' sample_analysis %>% dplyr::mutate(Area = setNa(Area)) -> sample_analysis
#' res <- missingPerCondition(sample_analysis,skylineconfig)
#' names(res)
#' res$data
#' res$figure
#' print(res$figure)
#' configuration <- skylineconfig$clone(deep=TRUE)
#' x <- sample_analysis
#'
missingPerCondition <- function(x, configuration, factors = configuration$table$fkeysLevel()){
  table <- configuration$table
  missingPrec <- interaction_missing_stats(x, configuration, factors)
  hierarchyKey <- tail(configuration$table$hierarchyKeys(),1)
  hierarchyKey <- paste0("nr_",hierarchyKey)
  xx <-missingPrec %>% group_by_at(c(table$isotopeLabel,
                                     factors,"nrNAs","nrReplicates")) %>%
    dplyr::summarize( !!sym(hierarchyKey) := n())

  formula <- paste(table$isotopeLabel, "~", paste(factors, collapse = "+"))
  #message(formula)

  nudgeval = max(xx[[hierarchyKey]]) * 0.05

  p <- ggplot(xx, aes_string(x= "nrNAs", y = hierarchyKey)) +
    geom_bar(stat="identity", color="black", fill="white") +
    geom_text(aes(label = !!sym(hierarchyKey)), nudge_y = nudgeval) +
    facet_grid(as.formula(formula))
  xx <- xx %>% tidyr::spread("nrNAs",hierarchyKey)

  return(list(data = xx ,figure = p))
}

# Functions - Handling isotopes ----

#' spreads isotope label heavy light into two columns
#' @export
#' @examples
#' setNa <- function(x){ifelse(x < 100, NA, x)}
#' data(sample_analysis)
#' data(skylineconfig)
#' sample_analysis %>% dplyr::mutate(Area = setNa(Area)) -> sample_analysis
#' x<-spreadValueVarsIsotopeLabel(sample_analysis,skylineconfig)
#' head(x)
#'
#' x<-spreadValueVarsIsotopeLabel(sample_analysis_HL,skylineconfig_HL)
#' head(x[,5:ncol(x)])
spreadValueVarsIsotopeLabel <- function(resData, configuration){
  table <- configuration$table
  idVars <- table$idVars()
  resData2 <- resData %>% dplyr::select(c(table$idVars(), table$valueVars()) )
  resData2 <- resData2 %>% tidyr::gather(variable, value, - idVars  )
  resData2 <- resData2 %>%  tidyr::unite(temp, table$isotopeLabel, variable )
  HLData <- resData2 %>% tidyr::spread(temp,value)
  invisible(HLData)
}

# Computing protein Intensity summaries ---


#' compute tukeys median polish from peptide or precursor intensities
#' @family matrix manipulation
#' @param name if TRUE returns the name of the summary column
#' @export
#' @examples
#' medpolishPly(name=T)
medpolishPly <- function(x, name=FALSE){
  if(name){
    return("medpolish")
  }
  X <- medpolish(x,na.rm=TRUE, trace.iter=FALSE, maxiter = 10);
  res <- tibble("sampleName" = names(X$col) , medpolish = X$col + X$overall)
  res
}

#' realign data
#' @export
reestablishCondition <- function(data,
                                 medpolishRes,
                                 configuration
){
  table <- configuration$table
  xx <- data %>%  dplyr::select(c(table$sampleName,
                                  table$factorKeys(),
                                  table$fileName,
                                  table$isotopeLabel)) %>% distinct()
  res <- dplyr::inner_join(xx,medpolishRes, by=table$sampleName)
  res
}

#' Summarizes the intensities within hierarchy
#' @param func - a function working on a matrix of intensities for each protein.
#' @return retuns function object
#'
#' @export
#' @importFrom purrr map
#' @examples
#'
#' library(LFQService)
#' library(tidyverse)
#' config <- LFQService::skylineconfig$clone(deep=TRUE)
#' data <- LFQService::sample_analysis
#' x <- intensity_summary_by_hkeys(data, config, func = medpolishPly)
#'
#' res <- x("unnest")
#' x("unnest")$data %>% dplyr::select(config$table$hierarchyKeys()[1] , "medpolish") %>% tidyr::unnest()
#' config <- LFQService::skylineconfig$clone(deep=TRUE)
#' config$table$hierarchyLevel <- 1
#' x <- intensity_summary_by_hkeys(data, config, func = medpolishPly)
#'
#' x("unnest")$data
#' xnested<-x()
#' dd <- x("plot")
#' dd$medpolishPly[[1]]
#' dd$plot[[2]]
#'
intensity_summary_by_hkeys <- function( data, config, func)
{
  x <- as.list( match.call() )
  makeName <- make.names(as.character(x$func))
  config <- config$clone(deep=TRUE)

  xnested <- data %>% group_by_at(config$table$hkeysLevel()) %>% nest()

  xnested <- xnested %>%
    dplyr::mutate(spreadMatrix = map(data, extractIntensities, config))
  xnested <- xnested %>%
    dplyr::mutate(!!makeName := map(spreadMatrix, func))
  xnested <- xnested %>%
    dplyr::mutate(!!makeName := map2(data,!!sym(makeName),reestablishCondition, config ))


  res <- function(value = c("nested","unnest","wide","plot")){
    value <- match.arg(value)
    if(value == "nested"){
      return(xnested)
    }else if(value == "unnest" || value =="wide"){
      unnested <- xnested %>% dplyr::select(config$table$hkeysLevel(), makeName) %>% tidyr::unnest()
      newconfig <- make_reduced_hierarchy_config(config,
                                                 workIntensity = func(name=TRUE),
                                                 hierarchy = config$table$hkeysLevel(names=FALSE))

      if(value == "wide"){
        wide <- LFQService::toWideConfig(unnested, newconfig)
        wide$config <- newconfig
        return(wide)
      }

      return(list(data = unnested, config = newconfig))


    }else if(value == "plot"){
      hierarchy_ID <- "hierarchy_ID"
      xnested <- xnested %>% tidyr::unite(hierarchy_ID , !!!syms(config$table$hkeysLevel()))
      figs <- xnested %>%
        dplyr::mutate(plot = map2(data, !!sym(hierarchy_ID) ,
                                  plot_hierarchies_line,
                                  factor_level = config$table$factorLevel, config ))
      #plot_hierarchies_add_quantline(figs$plot[[1]], figs$medpolishPly[[1]], "medpolish", config)
      figs <- figs %>%
        dplyr::mutate(plot = map2(plot, !!sym(makeName) ,
                                  plot_hierarchies_add_quantline, func(name=TRUE), config ))
      return(figs)
    }

  }
  return(res)
}

#' applys func - a function working on matrix for each protein.
#' @export
applyToHierarchyBySample <- function(data, config, func = medpolishPly, unnest = FALSE){
  stop("DEPRECATED! use intensity_summary_by_hkeys")

}

#' median polish from normalized peptide intensities
#' @export
#' @examples
#' resultsV12954 <- LFQService::resultsV12954
#' res <- medpolish_protein_quants(resultsV12954$pepIntensityNormalized,
#' resultsV12954$config_pepIntensityNormalized )
#'
#' dim(res("unnest")$data)
#'
medpolish_protein_quants <- function(data, config){
  protintensity <- LFQService::intensity_summary_by_hkeys(data ,
                                                          config,
                                                          medpolishPly)
  return(protintensity)
}

#' write intensites into folder - for the moment protein
#' @export
#'
protein_quants_write <- function(protintensity,
                                 path_qc,
                                 suffix="",
                                 na_fraction = 0.3){

  suffix <- paste0("_",suffix)
  message("writing protein intensity data into: ", path_qc)

  unnest <- protintensity("unnest")

  lfq_write_table(separate_factors(separate_hierarchy(unnest$data, unnest$config), unnest$config),
                   path = file.path(path_qc,paste0("protein_intensities_long",suffix,".csv")))


  lfq_write_table(separate_hierarchy(protintensity("wide")$data, protintensity("wide")$config),
                   path = file.path(path_qc,paste0("protein_intensities_wide",suffix,".csv")))

  lfq_write_table(protintensity("wide")$annotation,
                   path = file.path(path_qc,paste0("protein_intensities_file_annotation",suffix,".csv")))


  pdf(file.path(path_qc,paste0("protein_intensities_heatmap_correlation",suffix,".pdf")), width = 10, height = 10)
  plot_heatmap_cor(protintensity("unnest")$data, protintensity("unnest")$config)
  dev.off()

  pdf(file.path(path_qc,paste0("protein_intensities_heatmap",suffix,".pdf")), width = 10, height = 10)
  LFQService::plot_heatmap(protintensity("unnest")$data, protintensity("unnest")$config, na_fraction = na_fraction)
  dev.off()

  res <- plot_pca(protintensity("unnest")$data,protintensity("unnest")$config)
  pdf(file.path(path_qc,paste0("protein_intensities_PCA",suffix,".pdf")), width=6, height=6)
  print(res)
  dev.off()
}



#' applys func - a funciton workin on matrix for each protein and returning a vector of the same length as the number of samples
#' @export
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' data <- sample_analysis
#' config <- skylineconfig$clone(deep=TRUE)
#'
#' res <- summarize_cv(data, config)
#' res$CV <- res$sd/res$mean
summarize_cv <- function(data, config, all = TRUE){
  intsym <- sym(config$table$getWorkIntensity())


  hierarchyFactor <- data %>%
    dplyr::group_by(!!!syms( c(config$table$hierarchyKeys(), config$table$fkeysLevel()) )) %>%
    dplyr::summarize(n = n(),
                     not_na = sum(!is.na(!!intsym)),
                     sd = sd(!!intsym, na.rm = T),
                     mean=mean(!!intsym, na.rm = T)) %>%  dplyr::ungroup()

  hierarchyFactor <- hierarchyFactor %>% dplyr::mutate_at(config$table$fkeysLevel(), funs(as.character) )

  if(all){
    hierarchy <- data %>%
      dplyr::group_by(!!!syms( config$table$hierarchyKeys() )) %>%
      dplyr::summarize(n = n(),
                       not_na = sum(!is.na(!!intsym)),
                       sd = sd(!!intsym,na.rm = T),
                       mean=mean(!!intsym,na.rm = T))

    hierarchy <- dplyr::mutate(hierarchy, !!config$table$factorKeys()[1] := "All")
    hierarchyFactor <- bind_rows(hierarchyFactor,hierarchy)
  }
  if(config$parameter$is_intensity_transformed == FALSE){
    hierarchyFactor %>% dplyr::mutate(CV = sd/mean * 100) -> hierarchyFactor
  }
  return(hierarchyFactor)
}

#' summarize stats output
#' @export
#' @examples
#' config <- skylineconfig$clone(deep=TRUE)
#' data <- sample_analysis
#' stats_res <- summarize_cv(data, config)
#' head(stats_res)
#' summarize_cv_quantiles(stats_res, config)
#' summarize_cv_quantiles(stats_res, config, stats="CV")
summarize_cv_quantiles <- function(stats_res ,config, stats = c("sd","CV"), probs = c(0.1, 0.25, 0.5, 0.75, 0.9)){
  stats <- match.arg(stats)
  toQuantiles <- function(x, probs_i = probs) {
    tibble(probs = probs, quantiles = quantile(x, probs_i , na.rm = T))
  }
  q_column <- paste0(stats,"_quantiles")

  xx2 <- stats_res %>%
    dplyr::group_by(!!!syms(config$table$fkeysLevel())) %>%
    tidyr::nest()


  sd_quantile_res2 <- xx2 %>%
    dplyr::mutate( !!q_column := purrr::map(data, ~toQuantiles(.[[stats]]) ))  %>%
    dplyr::select(!!!syms(c(config$table$fkeysLevel(),q_column))) %>%
    tidyr::unnest()

  xx <- sd_quantile_res2 %>% tidyr::unite("interaction",config$table$fkeysLevel())
  wide <- xx %>%  spread("interaction", quantiles)
  return(list(long=sd_quantile_res2, wide= wide))
}


#' Compute theoretical sample sizes from factor level standard deviations
#' @export
#' @examples
#'
#' data <- sample_analysis
#' config <- skylineconfig$clone(deep=TRUE)
#'
#' data2 <- transform_work_intensity(data, config, transformation = log2)
#'
#' res <- lfq_power_t_test(data2, config)
#' res
#' stats_res <- summarize_cv(data2, config, all=FALSE)
#' res <- lfq_power_t_test(data2, config, delta=2)
#' res
#' res <- lfq_power_t_test(data2, config, delta=c(0.5,1,2))
#' res
#'
lfq_power_t_test <- function(data,
                             config,
                             delta = 1,
                             power = 0.8,
                             sig.level = 0.05,
                             probs = seq(0.5,0.9, by = 0.1)){

  if(!config$parameter$is_intensity_transformed){
    warning("Intensities are not transformed yet.")
  }

  stats_res <- summarize_cv(data, config, all=FALSE)
  sd <- na.omit(stats_res$sd)
  if(length(sd) > 0){
    quantilesSD <- quantile(sd,probs)

    sampleSizes <- expand.grid(probs = probs, delta = delta)
    quantilesSD <- quantile( sd, sampleSizes$probs )
    sampleSizes <- add_column( sampleSizes, sd = quantilesSD, .before=2 )
    sampleSizes <- add_column( sampleSizes, quantile = names(quantilesSD), .before=1 )

    getSampleSize <- function(sd, delta){power.t.test(delta = delta, sd=sd, power=power, sig.level=sig.level)$n}

    sampleSizes <- sampleSizes %>% mutate( N_exact = purrr::map2_dbl(sd, delta, getSampleSize))
    sampleSizes <- sampleSizes %>% mutate( N = ceiling(N_exact))
    sampleSizes <- sampleSizes %>% mutate( FC = round(2^delta, digits=2))

    summary <- sampleSizes %>% dplyr::select( -N_exact, -delta) %>% spread(FC, N, sep="=")
    return(list(long = sampleSizes, summary = summary))
  }else{
    message("no standard deviation is available, check if model is saturated (factor level variable).")
    return(NULL)
  }
}


#' applys func - a funciton workin on matrix for each protein and returning a vector of the same length as the number of samples
#' @export
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' data <- sample_analysis
#' config <- skylineconfig$clone(deep=TRUE)
#' res <- summarize_cv(data, config)
#' plot_stat_density(res, config, stat="mean")
#' plot_stat_density(res, config, stat="sd")
#' plot_stat_density(res, config, stat="CV")
plot_stat_density <- function(data, config, stat = c("CV","mean","sd"), ggstat = c("density", "ecdf")){
  stat <- match.arg(stat)
  ggstat <- match.arg(ggstat)
  p <- ggplot(data, aes_string(x = stat, colour = config$table$factorKeys()[1] )) +
    geom_line(stat = ggstat)
  return(p)
}
#'plot_stat_density_median
#'@export
plot_stat_density_median <- function(data, config, stat = c("CV","sd"), ggstat = c("density", "ecdf")){
  stat <- match.arg(stat)
  ggstat <- match.arg(ggstat)
  data <- data %>% dplyr::filter_at(stat, all_vars(!is.na(.)))
  res <- data %>% dplyr::mutate(top = ifelse(mean > median(mean, na.rm=TRUE),"top 50","bottom 50")) -> top50
  p <- ggplot(top50, aes_string(x = stat, colour = config$table$factorKeys()[1])) +
    geom_line(stat = ggstat) + facet_wrap("top")
  return(p)
}

#' applys func - a funciton workin on matrix for each protein and returning a vector of the same length as the number of samples
#' @export
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' data <- sample_analysis
#' config <- skylineconfig$clone(deep=TRUE)
#' res <- summarize_cv(data, config)
#' plot_stat_violin(res, config, stat="mean")
#' plot_stat_violin(res, config, stat="sd")
#' plot_stat_violin(res, config, stat="CV")
plot_stat_violin <- function(data, config, stat = c("CV","mean","sd")){
  stat <- match.arg(stat)
  p <- ggplot(data, aes_string(x = config$table$factorKeys()[1], y = stat  )) +
    geom_violin()
  return(p)
}
#' plot_stat_violin_median
#' @export
#'
plot_stat_violin_median <- function(data, config , stat=c("CV","sd")){

  median.quartile <- function(x){
    out <- quantile(x, probs = c(0.25,0.5,0.75))
    names(out) <- c("ymin","y","ymax")
    return(out)
  }
  data <- data %>% dplyr::filter_at(stat, all_vars(!is.na(.)))

  res <- data %>%
    dplyr::mutate(top = ifelse(mean > median(mean, na.rm = TRUE),"top 50","bottom 50")) ->
    top50

  p <- ggplot(top50, aes_string(x = config$table$factorKeys()[1], y = stat)) +
    geom_violin() +
    stat_summary(fun.y=median.quartile,geom='point', shape=3) + stat_summary(fun.y=median,geom='point', shape=1) +
    facet_wrap("top")
  return(p)
}

#' stddev vs mean
#' @export
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' data <- sample_analysis
#' config <- skylineconfig$clone(deep=TRUE)
#' res <- summarize_cv(data, config)
#' plot_stdv_vs_mean(res, config)
#' datalog2 <- transform_work_intensity(data, config, transformation = log2)
#' statlog2 <- summarize_cv(datalog2, config)
#' plot_stdv_vs_mean(statlog2, config)
#' config$table$getWorkIntensity()
#' config$table$popWorkIntensity()
#' datasqrt <- transform_work_intensity(data, config, transformation = sqrt)
#' ressqrt <- summarize_cv(datasqrt, config)
#' plot_stdv_vs_mean(ressqrt, config)
plot_stdv_vs_mean <- function(data, config){
  p <- ggplot(data, aes(x = mean, y = abs(sd))) +
    geom_point() +
    geom_smooth(method="loess") +
    facet_wrap(config$table$factorKeys()[1], nrow=1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(p)
}


.string.to.colors <- function(string, colors=NULL){
  if (is.factor(string)){
    string = as.character(string)
  }
  if (!is.null(colors)){
    if (length(colors)!=length(unique(string))){
      stop("The number of colors must be equal to the number of unique elements.")
    }
    else {
      conv = cbind(unique(string), colors)
    }
  } else {
    conv = cbind(unique(string), rainbow(length(unique(string))))
  }
  unlist(lapply(string, FUN=function(x){conv[which(conv[,1]==x),2]}))
}

.number.to.colors = function(value, colors=c("red", "blue"), num=100){
  cols = colorRampPalette(colors)(num)
  cols = 	cols[findInterval(value, vec=seq(from=min(value), to=max(value), length.out=num))]
  cols
}


#' plot correlation heatmap with annotations
#'
#' @export
#' @importFrom heatmap3 heatmap3
#' @examples
#' library(tidyverse)
#' data <- sample_analysis
#' config <- skylineconfig$clone( deep=TRUE )
#' # LFQService::plot_heatmap_cor( data, config )
#' # plot_heatmap_cor( data, config, R2=TRUE )
#'
plot_heatmap_cor <- function(data,
                             config,
                             R2 = FALSE,
                             distfun = function(x) as.dist(1 - cor(t(x), use = "pa")),
                             ...){
  res <-  toWideConfig(data, config , as.matrix = TRUE)
  annot <- res$annotation
  res <- res$data


  cres <- cor(res, use = "pa")
  if(R2){
    cres <- cres^2
  }

  factors <- dplyr::select_at(annot, config$table$factorKeys())
  ColSideColors <- as.matrix(dplyr::mutate_all(factors, funs(.string.to.colors)))
  rownames(ColSideColors) <- annot$sampleName
  res <- heatmap3::heatmap3(cres,symm=TRUE,
                            scale="none",
                            ColSideColors = ColSideColors,
                            margin = c(8,3),
                            main = ifelse(R2, "R^2", "correlation"),
                            ...=...)
  invisible(res)
}

#' plot heatmap with annotations
#' @export
#' @importFrom heatmap3 heatmap3
#' @examples
#' library(tidyverse)
#' library(LFQService)
#' data <- sample_analysis
#' config <- skylineconfig$clone(deep=TRUE)
#' #plot_heatmap(data, config)
plot_heatmap <- function(data, config, na_fraction = 0.4){
  res <-  toWideConfig(data, config , as.matrix = TRUE)
  annot <- res$annotation
  res <- res$data

  factors <- dplyr::select_at(annot, config$table$factorKeys())
  ColSideColors <- as.matrix(dplyr::mutate_all(factors, funs(.string.to.colors)))
  rownames(ColSideColors) <- annot$sampleName
  res <- quantable::removeNArows(res, round(ncol(res)*na_fraction,digits = 0))
  res <- t(scale(t(res)))
  res <- heatmap3::heatmap3(res,
                            ColSideColors = ColSideColors,
                            labRow="",
                            showRowDendro =FALSE,
                            margin = c(8,3), scale="none")
  invisible(res)
}

#' plot PCA
#' @export
#' @import ggfortify
#' @examples
#'
#' library(tidyverse)
#' library(LFQService)
#' data <- sample_analysis
#' config <- skylineconfig$clone(deep=TRUE)
#' plot_pca(data, config)
plot_pca <- function(data , config){
  wide <- toWideConfig(data, config ,as.matrix = TRUE)
  ff <- na.omit(wide$data)
  ff <- t(ff)
  ids <- wide$annotation
  res <- autoplot(prcomp(ff), data=ids, colour=config$table$fkeysLevel()[1],
                  shape=if(!is.na(config$table$fkeysLevel()[2])){
                    config$table$fkeysLevel()[2]
                  })
  return(res)
}

#' plot heatmap of NA values
#' @export
#' @importFrom heatmap3 heatmap3 showLegend
#' @examples
#'
#' library(tidyverse)
#' library(LFQService)
#' data <- sample_analysis
#' config <- skylineconfig$clone(deep=TRUE)
#' plot_NA_heatmap(data, config, cexCol=1)
plot_NA_heatmap <- function(data, config, showRowDendro=FALSE, cexCol=1, limitrows=10000 ){
  res <-  toWideConfig(data, config , as.matrix = TRUE)
  annot <- res$annotation
  res <- res$data
  stopifnot(annot$sampleName == colnames(res))

  factors <- dplyr::select_at(annot, config$table$factorKeys())
  ColSideColors <- as.matrix(dplyr::mutate_all(factors, list(.string.to.colors)))
  rownames(ColSideColors) <- annot$sampleName

  res[!is.na(res)] <- 0
  res[is.na(res)] <- 1
  allrows <- nrow(res)
  res <- res[apply(res,1, sum) > 0,]

  message("rows with NA's: ", nrow(res), "; all rows :", allrows, "\n")
  if(nrow(res) > 0){
    if(nrow(res) > limitrows ){
      message("limiting nr of rows to:", limitrows,"\n")
      resPlot <- res[sample( 1:nrow(res),limitrows),]
    }else{
      resPlot <- res
    }
    res_plot <- heatmap3::heatmap3(resPlot,
                                   distfun = function(x){dist(x, method="binary")},
                                   scale="none",
                                   col=c("white","black"),
                                   labRow = "",
                                   ColSideColors = ColSideColors,
                                   showRowDendro = showRowDendro,
                                   cexCol = cexCol,
                                   margin = c(8,3),
                                   legendfun = function()
                                     showLegend(legend=c("NA"),
                                                col=c("black"),
                                                cex=1.5))

    invisible(list(res = res, res_plot=res_plot))

  }
}
