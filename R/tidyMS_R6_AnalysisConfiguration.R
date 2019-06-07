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
                                       res <- (self$hierarchy[1:self$hierarchyLevel])
                                       return(ifelse(names, names(res), res))
                                     },
                                     factorKeys = function(){
                                       return(names(self$factors))
                                     },
                                     fkeysLevel = function(){
                                       res <- (self$factors[1:self$factorLevel])
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
                                       c( self$getWorkIntensity(), self$ident_qValue, self$ident_Score)
                                     }
                                   )
)
# AnalysisConfiguration ----
#' Analysis Configuration
#' @export
AnalysisConfiguration <- R6Class("AnalysisConfiguration",
                                 public = list(
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
#' make_interaction_column(xx, c("A","B"))
#' make_interaction_column(xx, c("A"))
make_interaction_column <- function(data, columns, sep="."){
  intr <- dplyr::select(data, columns)
  intr <- purrr::map2_dfc(columns, intr, paste0)
  colnames(intr) <- paste0("interaction_",columns)
  data$interaction <- interaction(intr, sep=sep)
  return(data)
}


#' create interaction column from factors
#' @export
#' @examples
#'
#' skylineconfig$table$factorKeys()
#' skylineconfig$table$factorLevel <- 1
#' make_interaction_column_config(LFQService::sample_analysis,skylineconfig)
make_interaction_column_config <- function(data, config){
  columns <- config$table$factorKeys()[1:config$table$factorLevel]
  data <- make_interaction_column(data, columns)
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
#' skylineconfig <- createSkylineConfiguration(isotopeLabel="Isotope.Label.Type", ident_qValue="Detection.Q.Value")
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#' data(skylinePRMSampleData)
#'
#' sample_analysis <- setup_analysis(skylinePRMSampleData, skylineconfig)
#'
setup_analysis <- function(data, configuration ,sep="~"){
  table <- configuration$table
  for(i in 1:length(table$hierarchy))
  {
    data <- tidyr::unite(data, UQ(sym(table$hierarchyKeys()[i])), table$hierarchy[[i]],remove = FALSE)
  }
  data <- dplyr::select(data , -dplyr::one_of(dplyr::setdiff(unlist(table$hierarchy), table$hierarchyKeys() )))

  for(i in 1:length(table$factors))
  {
    if( length(table$factors[[i]]) > 1){
      data <- tidyr::unite(data, UQ(sym(table$factorKeys()[i])), table$factors[[i]],remove = FALSE)
    }else{
      newname <-table$factorKeys()[i]
      data <- dplyr::mutate(data, !!newname := !!sym(table$factors[[i]]))
    }
  }

  sampleName <- table$sampleName

  if(!sampleName  %in% names(data)){
    message("creating sampleName")

    data <- data %>%  tidyr::unite( UQ(sym( sampleName)) , unique(unlist(table$factors)), remove = TRUE ) %>%
      dplyr::select(sampleName, table$fileName) %>% distinct() %>%
      dplyr::mutate_at(sampleName, function(x){ x<- make.unique( x, sep=sep )}) %>%
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
plot_hierarchies_line <- function(res, proteinName, configuration, factor_level=1, separate=FALSE){
  rev_hnames <- configuration$table$hierarchyKeys(TRUE)
  fragment <- rev_hnames[1]
  peptide <- rev_hnames[1]

  if(length(rev_hnames) > 2){
    peptide <- rev_hnames[2]
  }
  res <- LFQService:::plot_hierarchies_line_default(res, proteinName = proteinName,
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


#' add quantline to plot
#' @export
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
#' LFQService::plot_hierarchies_boxplot(xnested$data[[3]], xnested$protein_Id[[3]],conf )
plot_hierarchies_boxplot <- function(ddd, proteinName, config , factor_level = 1, boxplot=TRUE){

  stopifnot(factor_level <= length(config$table$factorKeys()))
  isotopeLabel <- config$table$isotopeLabel
  ddd %>%  tidyr::gather("factor","level",config$table$factorKeys()[1:factor_level]) -> ddlong
  ddlong <- as.data.frame(unclass(ddlong))

  p <- ggplot(ddlong, aes_string(x = config$table$hkeysLevel(rev = TRUE)[1],
                                 y = config$table$getWorkIntensity(),
                                 fill="level"
  )) +
    facet_wrap(~factor, ncol=1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle =45, hjust = 0.5, vjust=0.5)) +
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
#' hierarchy_counts_sample(sample_analysis, skylineconfig)
#'
hierarchy_counts_sample <- function(data,
                                    configuration)
{
  hierarchy <- configuration$table$hierarchyKeys()
  data %>% dplyr::filter(! is.na(!!sym(configuration$table$getWorkIntensity() ))) -> xx
  res <- xx %>% dplyr::group_by_at(c(configuration$table$isotopeLabel, configuration$table$sampleName)) %>%
    dplyr::summarise_at( hierarchy, n_distinct )
  return(res)
}
#' Light only version.
#' Summarize Protein counts
#' @export
#' @importFrom dplyr group_by_at
#' @examples
#' skylineconfig <- createSkylineConfiguration(isotopeLabel="Isotope.Label.Type", ident_qValue="Detection.Q.Value")
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
#' skylineconfig <- createSkylineConfiguration(isotopeLabel="Isotope.Label.Type", ident_qValue="Detection.Q.Value")
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#' data(skylinePRMSampleData)
#' configuration <- skylineconfig
#' sample_analysis <- setup_analysis(skylinePRMSampleData, skylineconfig)
#'
#'
#' summarizeHierarchy(sample_analysis, skylineconfig)
#' summarizeHierarchy(sample_analysis, skylineconfig, factors=character())
#'
#' summarizeHierarchy(sample_analysis, skylineconfig, hierarchy = skylineconfig$table$hkeysLevel() )
#' summarizeHierarchy(sample_analysis, skylineconfig, hierarchy = NULL, factors=skylineconfig$table$fkeysLevel() )
#' skylineconfig$table$hierarchyLevel=1
#' summarizeHierarchy(sample_analysis, skylineconfig, factors = skylineconfig$table$fkeysLevel())
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
  all_hierarchy <- configuration$table$hierarchyKeys()
  #factors <- configuration$table$factorKeys()[ifelse(factor_level < 1, 0, 1): factor_level]

  precursor <- x %>% dplyr::select(factors,all_hierarchy) %>% dplyr::distinct()
  x3 <- precursor %>% group_by_at(c(factors,hierarchy)) %>%
    dplyr::summarize_at( all_hierarchy,
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
#' xx <- completeCases(sample_analysis,skylineconfig)
#'
#' skylineconfig$parameter$qVal_individual_threshold <- 0.01
#' xx <- LFQService::removeLarge_Q_Values(sample_analysis, skylineconfig)
#' xx <- completeCases(xx, skylineconfig)
#' interaction_missing_stats(xx, skylineconfig)
#'
#' configuration <- skylineconfig
#' configuration$table$fkeysLevel()
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
#' @export
#'
summarize_missigness_impute <- function(mdataTrans, pepConfig){
  xx <- interaction_missing_stats(mdataTrans, pepConfig)
  xx <-  make_interaction_column(xx, pepConfig$table$fkeysLevel(), sep=":")

  #xx %>% mutate(perc_missing= nrNAs/nrReplicates*100) -> xx
  xx %>% mutate(nrMeasured = nrReplicates - nrNAs) -> xx


  lowerMean <- function(meanArea, probs = 0.1){
    meanAreaNotNA <- na.omit(meanArea)
    small10 <- meanAreaNotNA[meanAreaNotNA < quantile(meanAreaNotNA, probs= probs)]
    meanArea[is.na(meanArea)] <- mean(small10)
    return(meanArea)
  }

  xx <- xx %>% group_by(interaction) %>% mutate(imputed = lowerMean(meanArea,probs=0.2))
  xx %>% dplyr::select(-one_of(c("nrNAs", pepConfig$table$fkeysLevel()))) -> xx


  nrReplicates <- xx%>% dplyr::select( -meanArea, -nrMeasured, -imputed) %>%
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

  allTables <- list(meanArea= meanArea, nrMeasured = nrMeasured,  nrReplicates = nrReplicates, meanAreaImputed = meanAreaImputed)
  res <- purrr::reduce(allTables,inner_join)

  #  nrMeasured %>% dplyr::select(starts_with("interaction")) -> nrMeasuredM
  #  nrReplicates %>% dplyr::select(starts_with("interaction")) -> nrReplicatesM
  return(list(long = xx, wide=res))
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
    dplyr::mutate(cs = cumsum(nrTransitions))
  res <- xxcs  %>% dplyr::select(-nrTransitions)

  formula <- paste(table$isotopeLabel, "~", paste(factors, collapse = "+"))
  message(formula)
  p <- ggplot(res, aes(x= nrNAs, y = cs)) + geom_bar(stat="identity") +
    facet_grid(as.formula(formula))

  res <- res %>% tidyr::spread("nrNAs","cs")
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
missingPerCondition <- function(x, configuration, factors =configuration$table$fkeysLevel()){
  table <- configuration$table
  missingPrec <- interaction_missing_stats(x, configuration, factors)
  hierarchyKey <- tail(configuration$table$hierarchyKeys(),1)
  hierarchyKey <- paste0("nr_",hierarchyKey)
  xx <-missingPrec %>% group_by_at(c(table$isotopeLabel,
                                     factors,"nrNAs","nrReplicates")) %>%
    dplyr::summarize( !!sym(hierarchyKey) := n())
  formula <- paste(table$isotopeLabel, "~", paste(factors, collapse = "+"))
  message(formula)

  p <- ggplot(xx, aes_string(x= "nrNAs", y = hierarchyKey)) + geom_bar(stat="identity")+
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

#' applys func - a function working on matrix for each protein.
#'
#' Returning a vector of the same length as the number of samples.
#' Can be used to compute column summaries per protein
#' @export
#' @importFrom purrr map
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' config <- LFQService::skylineconfig$clone(deep=TRUE)
#' data <- LFQService::sample_analysis
#' x <- applyToHierarchyBySample(data, config, medpolishPly)
#'
#' x %>% dplyr::select(skylineconfig$table$hierarchyKeys()[1] ,  medpolishPly) %>% tidyr::unnest()
#' config <- LFQService::skylineconfig$clone(deep=TRUE)
#' x <- applyToHierarchyBySample(data, config, medpolishPly, hierarchy_level = 2, unnest=TRUE)
#' config <- LFQService::skylineconfig$clone(deep=TRUE)
#' x <- applyToHierarchyBySample(data, config, medpolishPly, hierarchy_level = 2)
applyToHierarchyBySample <- function( data, config, func, hierarchy_level = 1, unnest = FALSE)
{
  config$table$hierarchyLevel <- hierarchy_level
  x <- as.list( match.call() )
  makeName <- make.names(as.character(x$func))
  xnested <- data %>% group_by_at(config$table$hkeysLevel()) %>% nest()
  xnested <- xnested %>% dplyr::mutate(spreadMatrix = map(data, extractIntensities, config))
  xnested <- xnested %>% dplyr::mutate(!!makeName := map(spreadMatrix, func))
  xnested <- xnested %>% dplyr::mutate(!!makeName := map2(data,!!sym(makeName),reestablishCondition, config ))
  if(unnest){
    unnested <- xnested %>% dplyr::select(config$table$hkeysLevel(), makeName) %>% tidyr::unnest()
    newconfig <- make_reduced_hierarchy_config(config,
                                               workIntensity = func(name=TRUE),
                                               hierarchy = config$table$hkeysLevel(names=FALSE))
    return(list(unnested = unnested, newconfig = newconfig))
  }
  return(xnested)
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
    dplyr::group_by(!!!syms( c(config$table$hierarchyKeys(), config$table$factorKeys()[1]) )) %>%
    dplyr::summarize(n = n(), sd = sd(!!intsym, na.rm = T), mean=mean(!!intsym, na.rm = T)) %>%  dplyr::ungroup()

  hierarchyFactor <- hierarchyFactor %>% dplyr::mutate_at(config$table$factorKeys()[1], funs(as.character) )

  if(all){
    hierarchy <- data %>%
      dplyr::group_by(!!!syms( config$table$hierarchyKeys() )) %>%
      dplyr::summarize(n = n(), sd = sd(!!intsym,na.rm = T), mean=mean(!!intsym,na.rm = T))
    hierarchy <- dplyr::mutate(hierarchy, !!config$table$factorKeys()[1] := "All")
    hierarchyFactor <- bind_rows(hierarchyFactor,hierarchy)
  }
  hierarchyFactor %>% dplyr::mutate(CV = sd/mean*100) -> res
  return(res)
}

#' Compute theoretical sample sizes from factor level standard deviations
#' @export
#' @examples
#'
#' data <- sample_analysis
#' config <- skylineconfig$clone(deep=TRUE)
#'
#' res <- lfq_power_t_test(data, config)
lfq_power_t_test <- function(dataTransformed, config, delta = 1, power=0.8,sig.level=0.05){
  if(!config$parameter$is_intensity_transformed){
    warning("Intensities are not transformed yet.")
  }
  stats_res <- summarize_cv(dataTransformed, config, all=FALSE)
  sd <- na.omit(stats_res$sd)
  if(length(sd) > 0){
    quantilesSD <- quantile(sd,seq(0.2,1,length=9))
    getSampleSize <- function(sd){power.t.test(delta = delta, sd=sd, power=power, sig.level=sig.level)$n}
    N <- sapply(quantilesSD, getSampleSize)
    sampleSizes <- data.frame( quantile = names(N) ,
                               sd = round(quantilesSD, digits=3),
                               N_exact = round(N, digits=3),
                               N = ceiling(N))
    rownames(sampleSizes) <- NULL
  }else{
    sampleSizes <- data.frame(quantile = NA ,
                              sd =NA,
                              N_exact = NA,
                              N = NA)
  }
  return(sampleSizes)
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
plot_stat_density <- function(data, config, stat = c("CV","mean","sd")){
  stat <- match.arg(stat)
  p <- ggplot(data, aes_string(x = stat, colour = config$table$factorKeys()[1] )) +
    geom_line(stat="density")
  return(p)
}
#'plot_stat_density_median
#'@export
plot_stat_density_median <- function(data, config, stat = c("CV","sd")){
  stat <- match.arg(stat)
  data <- data %>% dplyr::filter_at(stat, all_vars(!is.na(.)))
  res <- data %>% dplyr::mutate(top = ifelse(mean > median(mean, na.rm=TRUE),"top 50","bottom 50")) -> top50
  p <- ggplot(top50, aes_string(x = stat, colour = config$table$factorKeys()[1])) +
    geom_line(stat = "density") + facet_wrap("top")
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


#' plot correlation heatmap with annations
#' @export
#' @importFrom heatmap3 heatmap3
#' @examples
#' library(tidyverse)
#' data <- sample_analysis
#' config <- skylineconfig$clone( deep=TRUE )
#' # LFQService::plot_heatmap_cor( data, config )
#' # plot_heatmap_cor( data, config, R2=TRUE )
#'
plot_heatmap_cor <- function(data, config, R2 = FALSE, distfun = function(x) as.dist(1 - cor(t(x), use = "pa"))){
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
                            main = ifelse(R2, "R^2", "correlation"))
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
plot_heatmap <- function(data, config){
  res <-  toWideConfig(data, config , as.matrix = TRUE)
  annot <- res$annotation
  res <- res$data

  factors <- dplyr::select_at(annot, config$table$factorKeys())
  ColSideColors <- as.matrix(dplyr::mutate_all(factors, funs(.string.to.colors)))
  rownames(ColSideColors) <- annot$sampleName
  res <- quantable::removeNArows(res, round(ncol(res)*0.4,digits = 0))
  res <- t(scale(t(res)))
  res <- heatmap3::heatmap3(res,
                            ColSideColors = ColSideColors,
                            labRow="",
                            showRowDendro =FALSE,
                            margin = c(8,3), scale="none")
  invisible(res)
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
plot_NA_heatmap <- function(data, config, showRowDendro=FALSE, cexCol=1 ){
  res <-  toWideConfig(data, config , as.matrix = TRUE)
  annot <- res$annotation
  res <- res$data
  stopifnot(annot$sampleName == colnames(res))

  factors <- dplyr::select_at(annot, config$table$factorKeys())
  ColSideColors <- as.matrix(dplyr::mutate_all(factors, list(.string.to.colors)))
  rownames(ColSideColors) <- annot$sampleName

  res[!is.na(res)] <- 0
  res[is.na(res)] <- 1
  res <- res[apply(res,1, sum) > 0,]
  if(nrow(res) > 0){
    res_plot <- heatmap3::heatmap3(res,
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
