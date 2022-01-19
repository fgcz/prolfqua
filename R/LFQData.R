#LFQData ----
#'
#' LFQData R6 class
#' @export
#' @family LFQData
#' @examples
#' istar <- prolfqua_data('data_ionstar')$filtered()
#'
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(data, istar$config)
#' tmp <- lfqdata$to_wide()
#' stopifnot("data.frame" %in% class(tmp$data))
#' tmp <- lfqdata$to_wide(as.matrix = TRUE)
#' stopifnot("matrix" %in% class(tmp$data))
#' lfqdata$factors()
#' stopifnot(lfqdata$is_transformed()==FALSE)
#' lfqdata$summarize_hierarchy()
#' lfqdata$omit_NA()
#'
#' stopifnot("LFQData" %in% class(lfqdata$get_copy()))
#' stopifnot("LFQDataTransformer" %in% class(lfqdata$get_Transformer()))
#' stopifnot("LFQDataStats" %in% class(lfqdata$get_Stats()))
#' stopifnot("LFQDataSummariser" %in% class(lfqdata$get_Summariser()))
#' stopifnot("LFQDataPlotter" %in% class(lfqdata$get_Plotter()))
#' stopifnot("LFQDataWriter" %in% class(lfqdata$get_Writer()))
#' stopifnot("LFQDataAggregator" %in% class(lfqdata$get_Aggregator()))
#'
#'
LFQData <- R6::R6Class(
  "LFQData",

  public = list(
    #' @field config AnalysisConfiguration
    config = NULL,
    #' @field data data.frame or tibble matching AnalysisConfiguration.
    data = NULL,
    #' @field is_pep todo
    is_pep = FALSE,
    #' @field prefix e.g. "peptide_", "protein_", "compound_"
    prefix = "",
    #' @description
    #' initialize
    #' @param data data.frame
    #' @param config configuration
    #' @param is_pep todo
    #' @param prefix will be use as output prefix
    initialize = function(data, config, is_pep=TRUE, prefix = "ms_") {
      self$data <- data
      self$config <- config$clone(deep = TRUE)
      self$is_pep <- is_pep
      self$prefix <- prefix
    },
    #' @description
    #' get deep copy
    get_copy = function(){
      return(self$clone(deep = TRUE))
    },
    #' @description
    #' samples subset of data
    #' @param size size of subset default 100
    #' @param seed set seed
    get_sample = function(size = 100, seed = NULL){
      if(!is.null(seed)){ set.seed( seed )}
      subset <- prolfqua::sample_subset(size = size, self$data, self$config)
      return(LFQData$new(subset, self$config$clone(deep = TRUE)))
    },
    #' @description
    #' get subset of data
    #' @param x data frame with columns containing subjectId
    get_subset = function(x){
      x <- select(x, all_of(self$subjectId())) |> distinct()
      subset <- inner_join(x, self$data)
      return(LFQData$new(subset, self$config$clone(deep = TRUE)))
    },
    #' @description
    #' get subject ID columns
    subjectId = function(){
      return(self$config$table$hkeysDepth())
    },
    #' @description
    #' is data trasfromed
    #' @param is_transformed logical
    #' @return logical
    is_transformed = function(is_transformed){
      if (missing(is_transformed)) {
        return(self$config$table$is_intensity_transformed)
      }else{
        self$config$table$is_intensity_transformed = is_transformed
      }
    },
    #' @description
    #' return name of intensity column
    #' @return name of intensity column
    intensity_column = function(){
      return(self$config$table$getWorkIntensity())
    },
    #' @description
    #' some software is reporting NA's as 0, you must remove it from your data
    #' @param threshold default 4.
    #' @return self
    remove_small_intensities = function(threshold = 4){
      self$data <- prolfqua::remove_small_intensities( self$data, self$config, threshold = threshold )
      invisible(self)
    },
    #' @description
    #' remove proteins with less than X peptides
    #' @return self
    filter_proteins_by_peptide_count = function(){
      message("removing proteins with less than",
              self$config$parameter$min_peptides_protein,
              " peptpides")
      self$data <- prolfqua::filter_proteins_by_peptide_count(self$data, self$config)$data
      invisible(self)
    },
    #' @description
    #' Omit NA from intensities per hierarchy (e.g. protein or peptide), idea is to use it for normalization
    #' For instance if a peptide has a missing value in more then nrNA of the samples within a condition
    #' it will be removed
    #' @param nrNA number of NA values
    #' @param factorDepth you control for nrNA per condition or experiment etc. e.g. factorDepth = 0  then per experiment
    #' @return LFQData with NA omitted.
    omit_NA = function(nrNA = 0, factorDepth = FALSE){
      if (!factorDepth)
      {
        missing <- interaction_missing_stats(self$data, self$config)
      } else {
        cfg <- self$config$clone(deep = TRUE)
        cfg$table$factorDepth <- factorDepth
        missing <- interaction_missing_stats(self$data, cfg)
      }
      notNA <- missing$data |> dplyr::filter(nrNAs <= nrNA)
      sumN <- notNA |> group_by_at(self$config$table$hierarchyKeys()) |>
        summarise(n = n())
      notNA <- sumN |> dplyr::filter(n == max(n))

      notNA <- notNA |> dplyr::select(self$config$table$hierarchyKeys())
      notNAdata <- dplyr::inner_join( notNA, self$data) |> ungroup()
      return(LFQData$new(notNAdata, self$config$clone(deep = TRUE)))
    },
    #'
    #' @description
    #' some software is reporting NA's as 0, you must remove it from your data
    #' @param threshold default 4.
    #' @return self
    complete_cases = function(){
      self$data <- prolfqua::complete_cases(self$data, self$config)
      invisible(self)
    },
    #' @description
    #' converts the data to wide
    #' @param as.matrix return as data.frame or matrix
    #' @return list with data, annotation, and configuration
    to_wide = function(as.matrix = FALSE){
      wide <- prolfqua::toWideConfig(self$data, self$config, as.matrix = as.matrix)
      wide$config <- self$config$clone(deep = TRUE)
      return(wide)
    },
    #' @description
    #' Annotation table.
    #' @return data.frame
    factors = function(){
      prolfqua::table_factors(self$data, self$config)
    },
    #' @description
    #' name of response variable
    #' @return data.frame
    response = function(){
      self$config$table$getWorkIntensity()
    },
    #' @description
    #' number of elements at each level
    hierarchy_counts = function(){
      prolfqua::hierarchy_counts(self$data, self$config)
    },
    #' @description
    #' e.g. number of peptides per protein etc
    #' @return data.frame
    summarize_hierarchy = function(){
      prolfqua::summarize_hierarchy(self$data, self$config)
    },
    #' @description
    #' get Plotter
    #' @return LFQDataPlotter
    get_Plotter = function(){
      return(LFQDataPlotter$new(self, self$prefix))
    },
    #' @description
    #' get Writer
    #' @param format array of formats to write to supported are xlsx, csv and html
    #' @return LFQDataPlotter
    get_Writer = function(format = "xlsx"){
      return(LFQDataWriter$new(self, self$prefix, format = format))
    },
    #' @description
    #' get Summariser
    #' @return LFQDataSummarizer
    get_Summariser = function(){
      return(LFQDataSummariser$new(self))
    },
    #' @description
    #' Get \code{\link{LFQDataStats}}. For more details see \code{\link{LFQDataStats}}.
    #' @param stats default interaction, computes statistics within interaction.
    #' @return LFQDataStats
    get_Stats = function(stats = c("everything","interaction", "all")){
      stats <- match.arg(stats)
      return(LFQDataStats$new(self, stats = stats))
    },
    #' @description
    #' get Stats
    #' @return LFQDataTransformer
    get_Transformer = function() {
      return(LFQDataTransformer$new(self))
    },
    #' @description
    #' get Aggregator
    #' @return LFQDataAggregator
    get_Aggregator = function() {
      return(LFQDataAggregator$new(self))
    },
    #' @description
    #' get difference of self with other if other is subset of self
    #'
    #' @details
    #' Use to compare filtering results obtained from self, e.g. which proteins and peptides were removed (other)
    #'
    #' @param other a filtered LFQData set
    #' @return LFQData
    #'
    filter_difference = function(other){
      diffdata <- prolfqua::filter_difference(self$data,other$data,self$config )
      res <- LFQData$new(diffdata , self$config$clone(deep = TRUE))
      return(res)
    }
  )
)

.annotProtein <- function(
  startdata,
  Accession,
  revpattern = "^REV_",
  contpattern = "zz"){

  GRP2 <- list()
  distinctprotid <- startdata |> select(pID = !!sym(Accession)) |> distinct()
  distinctprotid <- distinctprotid |> mutate(
    proteinAnnot = case_when(grepl(revpattern,pID) ~ "REV",
                             grepl(contpattern,pID) ~ "CON",
                             TRUE ~ "FW"))
  GRP2$percentOfContaminants <-  round(mean(distinctprotid$proteinAnnot == "CON") * 100, digits = 2)
  GRP2$percentOfFalsePositives <- round(mean(distinctprotid$proteinAnnot == "REV") * 100, digits = 2)
  GRP2$totalNrOfProteins <- sum( table( distinctprotid$proteinAnnot ) )
  GRP2$NrOfProteinsNoDecoys <- sum(distinctprotid$proteinAnnot == "FW")
  return(list(stats = GRP2, distinctprotid = distinctprotid))
}

# RowAnnotProtein ----
#' Decorates LFQData with a row annotation and some protein specific functions.
#'
#' @export
#' @family LFQData
#' @examples
#'
#' istar <- prolfqua_data('data_ionstar')$filtered()
#'
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(data, istar$config)
#' pannot <- RowAnnotProtein$new( lfqdata )
#'
#' pannot$annotateREV()
#' pannot$annotateCON()
#'
#' pannot$nr_clean()
#' dd <- pannot$clean()
#' tmp <- lfqdata$get_subset(dd)
#' tmp$complete_cases()
#'
RowAnnotProtein <-
  R6::R6Class("RowAnnotProtein",
          public = list(
            #' @field row_annot data.frame containing further information
            row_annot = NULL,
            #' @field pID column with protein ids
            pID = character(),
            #' @description initialize
            #' @param lfqdata data frame from \code{\link{setup_analysis}}
            #' @param row_annot data frame with row annotation. Must have columns matching \code{config$table$hkeysDepth()}
            initialize = function(lfqdata, row_annot){
              stopifnot(lfqdata$config$table$hierarchyDepth == 1)

              self$pID = lfqdata$config$table$hkeysDepth()
              if (!missing(row_annot)) {
                stopifnot(self$pID %in% colnames(row_annot))
                self$row_annot <- row_annot
              } else {
                self$row_annot <- distinct(select(lfqdata$data,self$pID))
              }
            },
            #' @description
            #' annotate rev sequences
            #' @param pattern default "REV_"
            annotateREV = function(pattern = "REV_") {
              self$row_annot <- self$row_annot |> mutate(
                REV = case_when(grepl(pattern, !!sym(self$pID)) ~ TRUE,
                                         TRUE ~ FALSE))

              return(sum(self$row_annot$REV))
            },
            #' @description
            #' annotate contaminants
            #' @param pattern default "^zz|^CON"
            annotateCON = function(pattern = "^zz|^CON") {
              self$row_annot <- self$row_annot |> mutate(
                CON = case_when(grepl(pattern, !!sym(self$pID)) ~ TRUE,
                                TRUE ~ FALSE))
              return(sum(self$row_annot$CON))
            },
            #' @description
            #' return number of cleans
            nr_clean = function(){
              if (!("REV" %in% colnames(self$row_annot)) ) { stop("annotate REV") }
              if (!("CON" %in% colnames(self$row_annot)) ) { stop("annotate CON") }
              return(sum(!self$row_annot$REV & !self$row_annot$CON))
            },
            #' @description
            #' remove REV and CON sequences
            clean = function(){
              if (!("REV" %in% colnames(self$row_annot)) ) { stop("annotate REV") }
              if (!("CON" %in% colnames(self$row_annot)) ) { stop("annotate CON") }
              return(filter(self$row_annot , !self$row_annot$REV & !self$row_annot$CON) )
            }

          )
  )



# LFQDataTransformer ----

#' Decorate LFQData with Methods for transforming Intensities
#'
#' @export
#'
#' @examples
#' istar <- prolfqua_data('data_ionstar')$filtered()
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(data, istar$config)
#'
#' lfqcopy <- lfqdata$get_copy()
#' lfqTrans <- lfqcopy$get_Transformer()
#'
#' x <- lfqTrans$intensity_array(log2)
#' x$lfq$config$table$is_intensity_transformed
#' x <- x$intensity_matrix(robust_scale)
#' plotter <- x$lfq$get_Plotter()
#' plotter$intensity_distribution_density()
#'
#' # transform by asinh root and scale
#' lfqcopy <- lfqdata$get_copy()
#' lfqTrans <- lfqcopy$get_Transformer()
#' x <- lfqTrans$intensity_array(asinh)
#' mads1 <- mean(x$get_scales()$mads)
#' x <- lfqTrans$intensity_matrix(robust_scale)
#' mads2 <- mean(x$get_scales()$mads)
#' stopifnot(abs(mads1 - mads2) < 1e-8)
#' stopifnot(abs(mean(x$get_scales()$medians)) < 1e-8)
#' lfqcopy <- lfqdata$get_copy()
#' lfqTrans <- lfqcopy$get_Transformer()
#' lfqTrans$log2()
#' before <- lfqTrans$get_scales()
#' lfqTrans$robscale()
#' after <- lfqTrans$get_scales()
#' stopifnot(abs(mean(before$medians) - mean(after$medians)) < 1e-8)
#' stopifnot(abs(mean(before$mads) - mean(after$mads)) < 1e-8)
#'
#' # normalize data using vsn
#' lfqcopy <- lfqdata$get_copy()
#' lfqTrans <- lfqcopy$get_Transformer()
#' lfqTransCheck <- lfqcopy$get_Transformer()
#'
#' lfqTransCheck$log2()
#' lfqTransCheck$get_scales()
#' lfqTransCheck$lfq$get_Plotter()$intensity_distribution_density()
#'
#' if(require("vsn")){
#'  res <- lfqTrans$intensity_matrix( .func = vsn::justvsn)
#'  res$lfq$get_Plotter()$intensity_distribution_density()
#'  res$get_scales()
#' }
#'
#'
LFQDataTransformer <- R6::R6Class(
  "LFQDataTransformer",
  public = list(
    #' @field lfq LFQData object
    lfq = NULL,
    #' @description
    #' initialize
    #' @param lfqdata
    #' LFQData object to transform
    initialize = function(lfqdata){
      self$lfq = lfqdata$clone(deep = TRUE)
    },
    #' @description
    #' log2 transform data
    #' @param force if FALSE, then data already log2 transformed will not be transformed a second time. TRUE force log transformation.
    #' @return LFQDataTransformer
    log2 = function(force = FALSE){
      if (self$lfq$is_transformed() == FALSE | force ) {
        self$lfq$data  <-  prolfqua::transform_work_intensity(self$lfq$data ,
                                                              self$lfq$config, log2)
        self$lfq$is_transformed(TRUE)
      } else {
        warning("data already transformed. If you still want to log2 tranform, set force = TRUE")
      }
      invisible(self)

    },
    #' @description
    #' get mean and variance and standard deviation in each sample
    #' @return list with means and mads
    get_scales = function()
    {
      get_robscales(self$lfq$data, self$lfq$config)
    },
    #' @description
    #' robust scale data
    #' @param colname new name of transformed column
    #' @param preserveMean should original mean value be preserved TRUE, if FALSE then center at zero
    #' @return LFQDataTransformer (self)
    robscale = function(preserveMean = TRUE, colname = "transformedIntensity"){
      res <- self$robscale_subset(self$lfq, preserveMean = preserveMean, colname = colname)
      invisible(res)
    },
    #' @description
    #' log2 transform and robust scale data based on subset
    #' @param lfqsubset LFQData subset to use for normalization
    #' @param preserveMean should original mean value be preserved TRUE, if FALSE then center at zero
    #' @param colname - how to name the transformed intensities, default transformedIntensity
    #' @return LFQDataTransformer (self)
    #'
    robscale_subset = function(lfqsubset,
                               preserveMean = TRUE,
                               colname = "transformedIntensity"){
      message("data is : ",self$lfq$is_transformed())
      if (self$lfq$is_transformed() != lfqsubset$is_transformed()) {
        warning("the subset must have the same config as self")
        invisible(NULL)
      }
      scales <- prolfqua::scale_with_subset(self$lfq$data,
                                            lfqsubset$data,
                                            self$lfq$config,
                                            preserveMean = preserveMean)
      self$lfq$data  <- scales$data
      if (!is.null(colname)) {
        self$lfq$data <- self$lfq$data |>
          dplyr::rename(!!colname := self$lfq$config$table$getWorkIntensity())
        self$lfq$config$table$popWorkIntensity()
        self$lfq$config$table$setWorkIntensity(colname)
      }
      invisible(self)

    },
    #' @description
    #' Transforms intensities
    #' @param .func transformation function working with arrays e.g. log2, log10, asinh etc.
    #' @return LFQDataTransformer (self)
    #'
    intensity_array = function(.func = log2) {
      .call <- as.list( match.call() )
      r <- prolfqua::transform_work_intensity(
        self$lfq$data,
        self$lfq$config,
        .func = .func,
        .funcname = deparse(.call$.func))
      self$lfq$data <- r
      invisible(self)
    },
    #' @description
    #' pass a function which works with matrices see robust_scale for examples.
    #' @param .func any function taking a matrix and returning a matrix (columns sample, rows feature e.g. base::scale) default robust_scale
    #' @return LFQDataTransformer (self)
    #'
    intensity_matrix = function(.func = robust_scale){
      .call <- as.list( match.call() )
      r <- prolfqua::applyToIntensityMatrix(
        self$lfq$data,
        self$lfq$config,
        .func = .func,
        .funcname = deparse(.call$.func))
      self$lfq$data <- r
      invisible(self)
    }
  )
)

# LFQDataStats-----
#'
#' Decorates LFQData with methods to compute statistics of interactions
#'
#' compute stdv, mean and CV per peptide or protein and condition.
#'
#' @export
#' @family LFQData
#' @examples
#'
#' # study variance of not normalized data
#' #source("c:/Users/wewol/prog/prolfqua/R/LFQData.R")
#' runallfuncs <- function(x){
#'
#'   stopifnot("data.frame" %in% class(x$stats()))
#'   stopifnot(c("long", "wide") %in% names(x$stats_quantiles()))
#'   stopifnot("ggplot" %in% class(x$density()))
#'   stopifnot("ggplot" %in% class(x$density_median()))
#'   stopifnot("ggplot" %in% class(x$density("ecdf")))
#'   stopifnot("ggplot" %in% class(x$density_median("ecdf")))
#'   stopifnot("ggplot" %in% class(x$violin()))
#'   stopifnot("ggplot" %in% class(x$violin_median()))
#'   stopifnot("ggplot" %in% class(x$stdv_vs_mean(size = 400)))
#'   if(!x$lfq$is_transformed()){
#'     stopifnot(is.null(x$power_t_test()))
#'     stopifnot(is.null(x$power_t_test_quantiles()))
#'   }
#' }
#' istar <- prolfqua_data('data_ionstar')$filtered()
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(data, istar$config)
#' lfqstats <- lfqdata$get_Stats()
#' lfqstats$violin()
#' runallfuncs(lfqstats)
#' x<-lfqstats
#'
#' #study variance of normalized data
#'
#' istar <- prolfqua_data('data_ionstar')$normalized()
#' istar$config$table$is_intensity_transformed
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(data, istar$config)
#' lfqdata$is_transformed(TRUE)
#' lfqstats <- lfqdata$get_Stats()
#' runallfuncs(lfqstats)
#'
#' #Slightly different dataset
#'
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' stopifnot(nrow(bb$data) == 25780)
#' config <- bb$config$clone(deep = TRUE)
#' analysis <- bb$data
#'
#' lfqdata <- LFQData$new(analysis, config)
#' # estimates statistics for all samples
#' lfqstats <- lfqdata$get_Stats(stats = "all")
#' runallfuncs(lfqstats)
#' lfqstats <- lfqdata$get_Stats(stats = "everything")
#' runallfuncs(lfqstats)
LFQDataStats <- R6::R6Class(
  "LFQDataStats",
  public = list(
    #' @field lfq LFQData
    lfq = NULL,
    #' @field stat either CV or sd (if is_transformed)
    stat = "CV",
    #' @field statsdf frame with statistics.
    statsdf = NULL,
    #' @description
    #' create analyse variances and CV
    #' @param lfqdata LFQData object
    #' @param stats if interaction - within group stats, if all then overall CV, if pooled - then pooled variance using grouping information (t.b.d.)
    initialize = function(lfqdata, stats = c("everything", "interaction", "all")){
      stats <- match.arg(stats)
      self$lfq = lfqdata
      self$stat <- if (!self$lfq$is_transformed()) {"CV"} else {"sd"}
      if (stats == "interaction" ) {
        self$statsdf <- prolfqua::summarize_stats(self$lfq$data, self$lfq$config)
      } else if (stats == "all" ) {
        self$statsdf <-
          prolfqua::summarize_stats_all(self$lfq$data, self$lfq$config)
      } else if (stats == "everything" ) {

        self$statsdf <- bind_rows(
          prolfqua::summarize_stats(self$lfq$data, self$lfq$config),
          prolfqua::summarize_stats_all(self$lfq$data, self$lfq$config)
          )
      }
    },
    #' @description
    #' access data.frame with statistics
    #' @return data.frame with computed statistics
    stats = function(){
      self$statsdf
    },
    #' @description
    #' Determine CV or sd for the quantiles
    #' @param probs for which quantile to determine CV or sd
    stats_quantiles = function(probs = c(0.1, 0.25, 0.5, 0.75, 0.9)){
      res <- prolfqua::summarize_stats_quantiles(
        self$stats(),
        self$lfq$config,
        stats = self$stat,
        probs = probs)
      return(res)
    },
    #' @description
    #' plots density or ecdf
    #' @param ggstat either density or ecdf
    #' @return ggplot
    density = function(ggstat = c("density", "ecdf")){
      prolfqua::plot_stat_density(
        self$stats(),
        self$lfq$config,
        stat = self$stat,
        ggstat = ggstat)
    },
    #' @description
    #' plot density or ecdf of CV or sd for the 50% of low intensity data and 50% of high intensity data
    #' @param ggstat either density of ecdf
    #' @return ggplot
    density_median = function(ggstat = c("density", "ecdf")){
      prolfqua::plot_stat_density_median(self$stats(), self$lfq$config, stat = self$stat)
    },
    #' @description
    #' plot violinplot of CV or sd
    #' @param ggstat either density of ecdf
    #' @return ggplot
    violin = function(){
      prolfqua::plot_stat_violin(self$stats(), self$lfq$config, stat = self$stat)
    },
    #' @description
    #' plot violinplot of CV or sd for the 50% of low intensity data and 50% of high intensity data
    #'
    #' @return ggplot
    #'
    violin_median = function(){
      prolfqua::plot_stat_violin_median(self$stats(), self$lfq$config, stat = self$stat)
    },
    #' @description
    #' plot sd vs mean
    #' @param size number of points to sample (default 200)
    #' @return ggplot
    #'
    stdv_vs_mean = function(size= 200){
      prolfqua::plot_stdv_vs_mean(self$stats(), self$lfq$config, size = size)
    },
    #' @description
    #' compute sample size for entire dataset
    #' @param probs quantiles of sd for which sample size should be computed
    #' @param delta effect size
    #' @param power power of test
    #' @param sig.level significance level.
    power_t_test_quantiles = function(
      probs = c(0.1, 0.25, 0.5, 0.75, 0.9),
      delta = c(0.59,1,2),
      power = 0.8,
      sig.level = 0.05)
    {
      if (!self$lfq$is_transformed()) {
        warning("data is not transformed - aborting")
        return()
      }
      res <- self$stats_quantiles(probs)
      res <- lfq_power_t_test_quantiles_V2(res$long,
                                           delta = delta,
                                           power = power,
                                           sig.level = sig.level )
      return(res)
    },
    #' @description
    #' compute sample for each protein
    #' @param delta effect size
    #' @param power power of test
    #' @param sig.level significance level.
    power_t_test = function(
      delta = c(0.59,1,2),
      power = 0.8,
      sig.level = 0.05
    ){
      if (!self$lfq$is_transformed()) {
        warning("data is not transformed - aborting")
        return()
      }

      res <- prolfqua::lfq_power_t_test_proteins(self$stats(),
                                                 delta = delta,
                                                 power = power,
                                                 sig.level = sig.level,
                                                 min.n = 1.5)
      return(res)
    }
  )
)

# LFQDataSummariser ----
#' Summarize LFQData
#'
#' @export
#' @family LFQData
#' @examples
#'
#' istar <- prolfqua_data('data_ionstar')
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(data, istar$config)
#' sum <- lfqdata$get_Summariser()
#' sum
#' sum$hierarchy_counts()
#' sum$hierarchy_counts_sample("wide")
#' sum$hierarchy_counts_sample("long")
#' sum$plot_hierarchy_counts_sample()
#' sum$interaction_missing_stats()
#' sum$missingness_per_condition()
#' sum$missingness_per_condition_cumsum()
#' sum$plot_missingness_per_condition()
#' sum$plot_missingness_per_condition_cumsum()
#'
LFQDataSummariser <- R6::R6Class(
  "LFQDataSummariser",
  public = list(
    #' @field lfq LFQData
    lfq = NULL,
    #' @description
    #' initialize
    #' @param lfqdata LFQData
    initialize = function(lfqdata ) {
      self$lfq <- lfqdata
    },
    #' @description
    #' summarize hierarchy
    hierarchy_counts = function(){
      self$lfq$hierarchy_counts()
    },
    #' @description
    #' number of elements at each level in every sample
    #' @param value wide - wide format, long - long format, plot - ggplot
    hierarchy_counts_sample = function(value=c("wide","long")){
      value <- match.arg(value)
      fun <- prolfqua::hierarchy_counts_sample(self$lfq$data, self$lfq$config)
      return(fun(value))
    },
    #' @description
    #' barplot showing number of elements at each level in every sample
    #' @param value wide - wide format, long - long format, plot - ggplot
    plot_hierarchy_counts_sample = function(){
      fun <- prolfqua::hierarchy_counts_sample(self$lfq$data, self$lfq$config)
      return(fun("plot"))
    },
    #' @description
    #' missing per condition and protein
    interaction_missing_stats = function(){
      prolfqua::interaction_missing_stats(self$lfq$data, self$lfq$config)
    },
    #' @description
    #' missing stats per condition
    missingness_per_condition = function(){
      prolfqua::missingness_per_condition(self$lfq$data, self$lfq$config)$data
    },
    #' @description
    #' missing stats per condition as cumulative sum
    missingness_per_condition_cumsum = function(){
      prolfqua::missingness_per_condition_cumsum(self$lfq$data, self$lfq$config)$data
    },
    #' @description
    #' barplot with number of features with 1,2, etc missing in condition
    #' @return ggplot
    plot_missingness_per_condition = function(){
      prolfqua::missingness_per_condition(self$lfq$data, self$lfq$config)$figure
    },
    #' @description
    #' barplot with cumulative sum of features with 1,2, etc missing in condition
    #' @return ggplot
    plot_missingness_per_condition_cumsum = function(){
      prolfqua::missingness_per_condition_cumsum(self$lfq$data, self$lfq$config)$figure
    }
  )
)

# LFQDataPlotter ----
#' Create various visualization of the LFQdata
#' @export
#'
#' @family LFQData
#' @import dplyr
#' @examples
#'
#' #library(prolfqua)
#' data('data_IonstarProtein_subsetNorm')
#' istar <- data_IonstarProtein_subsetNorm
#'
#' istar$data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#'
#' lfqdata <- LFQData$new(
#'  istar$data,
#'  istar$config)
#' #LFQDataPlotter$debug("pca_plotly")
#' lfqplotter <- lfqdata$get_Plotter()
#' graphics.off()
#' lfqplotter$heatmap()
#' graphics.off()
#' lfqplotter$heatmap_cor()
#' graphics.off()
#' lfqplotter$pca()
#' lfqplotter$pca_plotly()
#'
#' #tmp <- lfqplotter$boxplots()
#' #tmp$boxplot[[1]]
#' lfqplotter$missigness_histogram()
#'
#' lfqplotter$NA_heatmap()
#' lfqplotter$intensity_distribution_density()
#' lfqplotter$intensity_distribution_violin()
#' lfqplotter$pairs_smooth()
#' lfqplotter$sample_correlation()
#' lfqplotter$raster()
#' prolfqua::plot_sample_correlation(istar$data, istar$config)
#'
LFQDataPlotter <- R6::R6Class(
  "LFQDataPlotter",
  list(
    #' @field lfq LFQData object
    #' @field prefix prefix to figure names when writing, e.g. protein_
    #' @field file_paths_pdf with paths to figures
    #' @field file_paths_html with paths to figures
    lfq = NULL,
    prefix = "",
    file_paths_pdf = list(),
    file_paths_html = list(),
    #' @description
    #' create LFQDataPlotter
    #' @param lfqdata LFQData
    #' @param prefix will be prepended to outputs written
    initialize = function(lfqdata, prefix = "ms_"){
      self$lfq = lfqdata
      self$prefix = prefix
      self$lfq$data <- na.omit(self$lfq$data)
    },

    #' @description
    #' plot intensities in raster
    #' @param arrange arrange by either mean or var
    #' @param not_na TRUE arrange by number of NA's, FALSE by arrange by intensity
    #' @param rownames show rownames (default FALSE - do not show.)
    #' @return ggplot
    raster = function(arrange = c("mean", "var") ,
                      not_na = FALSE,
                      rownames = FALSE ){
      arrange <- match.arg(arrange)
      fig <- prolfqua::plot_raster(self$lfq$data,
                                   self$lfq$config,
                                   arrange = arrange,
                                   not_na = not_na,
                                   show_rownames = rownames)
      return(fig)
    },
    #' @description
    #' heatmap of intensities
    #' @param na_fraction max fraction of NA's per row
    #' @param rownames show rownames (default FALSE - do not show.)
    #' @return pheatmap
    heatmap = function(na_fraction = 0.3, rownames = FALSE){
      fig <- prolfqua::plot_heatmap(self$lfq$data,
                                    self$lfq$config,
                                    na_fraction = na_fraction,
                                    show_rownames = rownames)
      return(fig)
    },
    #' @description
    #' heatmap of sample correlations
    #' @return pheatmap
    heatmap_cor = function(){
      fig <- prolfqua::plot_heatmap_cor(self$lfq$data, self$lfq$config)
      return(fig)
    },
    #' @description
    #' PCA plot
    #' @param add_txt show sample names
    #' @return ggplot
    pca = function(add_txt = TRUE){
      fig <- prolfqua::plot_pca(self$lfq$data, self$lfq$config, add_txt = add_txt)
      return(fig)
    },
    #' @description
    #' pca plot
    #' @param add_txt show sample names
    #' @return plotly
    pca_plotly = function(add_txt = FALSE){
      fig <- plotly::ggplotly(self$pca(add_txt = add_txt), tooltip = self$lfq$config$table$sampleName)
      return(fig)
    },
    #' @description
    #' boxplots for all proteins
    #' @return tibble with column boxplots containing ggplot objects
    boxplots = function(){
      bb <- prolfqua::plot_hierarchies_boxplot_df(self$lfq$data, self$lfq$config)
      return(bb)
    },
    #' @description
    #' histogram of intensities given number of missing in conditions
    #' @return ggplot
    missigness_histogram = function(){
      prolfqua::missigness_histogram(self$lfq$data, self$lfq$config)
    },

    #' @description
    #' heatmap of features with missing values
    #' @return ggplot
    NA_heatmap = function(){
      prolfqua::plot_NA_heatmap(self$lfq$data, self$lfq$config)
    },
    #' @description
    #' density distribution of intensities
    #' @param legend show legend TRUE, FALSE do not show.
    #' @return ggplot
    intensity_distribution_density = function(legend = TRUE){
      prolfqua::plot_intensity_distribution_density(self$lfq$data, self$lfq$config, legend = legend)
    },
    #' @description
    #' Violinplot showing distribution of intensities in all samples
    #' @return ggplot
    intensity_distribution_violin = function(){
      prolfqua::plot_intensity_distribution_violin(self$lfq$data, self$lfq$config)
    },
    #' @description
    #' pairsplot of intensities
    #' @param max maximal number of samples to show
    #' @return NULL
    pairs_smooth = function(max=10){
      dataTransformed <- self$lfq$data
      config <- self$lfq$config
      samples <- dplyr::select(self$lfq$data, config$table$sampleName) |>
        distinct() |>
        pull()
      if (length(samples) > max) {
        limit <- samples |> sample(max)
        ldata <- dataTransformed |>
          dplyr::filter(!!sym(config$table$sampleName) %in% limit)
        prolfqua::pairs_smooth( prolfqua::toWideConfig(ldata, config, as.matrix = TRUE)$data )
      }else{
        prolfqua::pairs_smooth( prolfqua::toWideConfig(dataTransformed, config, as.matrix = TRUE)$data )
      }
      NULL
    },
    #' @description
    #' plot of sample correlations
    #' @return NULL
    sample_correlation = function(){
      prolfqua::plot_sample_correlation(self$lfq$data, self$lfq$config)
    },
    #' @description
    #' write boxplots to file
    #' @param path_qc path to write to
    #' @param width fig width
    #' @param height fig height
    #'
    write_boxplots = function(path_qc, width = 6, height = 6){
      fpath <- file.path(path_qc,paste0(self$prefix, "boxplot.pdf"))
      message("generating boxplots")
      bb <- self$boxplots()

      message("writing ", fpath)
      pb <- progress::progress_bar$new(total = length(bb$boxplot))

      pdf(fpath, width = width, height = height)
      lapply(bb$boxplot, function(x){pb$tick(); print(x)})
      dev.off()
    },
    #' @description
    #' write pltly figures to path_qc
    #' @keywords static
    #' @param fig pltly figure
    #' @param path_qc path to write to
    #' @param fig_name file name (without extension)
    #' @return path the file was written to.
    write_pltly = function(fig,
                           path_qc,
                           fig_name){
      fname <- paste0(self$prefix, fig_name,".html")
      html_path <- file.path(".", path_qc, fname)
      message("writing ", html_path)
      htmlwidgets::saveWidget(widget = fig,
                              file = fname )
      file.rename(fname, html_path)
      self$file_paths_html[[fig_name]] <- html_path
      invisible(html_path)
    },
    #' @description
    #' write figure to pdf
    #' @param fig ggplot or pheatmap
    #' @param path_qc path to write to
    #' @param fig_name name of figure (no extension)
    #' @param width figure width
    #' @param height figure height
    #' @return path the file was written to
    write_pdf = function(fig,
                         path_qc,
                         fig_name,
                         width=7 , height=7 ){
      fpath <- file.path(path_qc,paste0(self$prefix,fig_name,".pdf"))
      message("writing ", fpath)
      graphics.off()
      pdf(fpath,
          width = width,
          height = height)
      print(fig)
      graphics.off()
      self$file_paths_pdf[[fig_name]] <- fpath
      invisible(fpath)
    },
    #' @description
    #' write heatmaps and pca plots to files
    #' @param path_qc path to write to
    #'
    write = function(path_qc){

      self$write_pdf(self$heatmap_cor(),
                     path_qc,
                     "intensities_heatmap_correlation",
                     width = 10, height = 10)
      self$write_pdf(self$heatmap(),
                     path_qc,
                     "intensities_heatmap",
                     width = 10, height = 10)
      self$write_pdf(self$pca(),
                     path_qc,
                     "intensities_PCA")
      self$write_pltly(self$pca_plotly(),
                       path_qc,
                       "intensities_PCA")

    }
  ))

# LFQDataWriter -----
#' Write LFQData, or provide outputs for writing.
#' @family LFQData
#'
#' returns long and wide format for writing
#' @export
#'
LFQDataWriter <- R6::R6Class(
  "LFQDataWriter",list(
    #' @field lfq LFQData
    #' @field format format to write to
    #' @field prefix prefix of filename
    #' @field file_paths list with paths were data was written to.
    #'
    lfq = NULL,
    format = "",
    prefix = "",
    file_paths = list(),
    #' @description
    #' initialize class
    #' @param lfqdata LFQData
    #' @param prefix prefix files with
    #' @param format which format to write to ("xlsx", "csv", "html")
    initialize = function(lfqdata,  prefix = "ms_", format="xlsx"){
      self$lfq = lfqdata
      self$format = format
      self$prefix = prefix
    },
    #' @description
    #' get Data in long format for writing
    #' @return tibble
    get_long = function(){
      #' gets data formatted for writing
      separate_factors(
        separate_hierarchy(self$lfq$data, self$lfq$config),
        self$lfq$config)
    },
    #' @description
    #' get Data in Wide format for writing
    #' @return list with data and annotation
    get_wide = function(){
      #' gets data formatted for writing
      wide <- self$lfq$to_wide()
      res <- list(data = separate_hierarchy(wide$data, self$lfq$config), annotation = wide$annotation)
      return(res)
    },
    #' @description
    #' write data to file
    #' @param path_qc path to write to
    write_long = function(path_qc) {
      fname <- paste0(self$prefix,"intensities_long")
      self$file_paths[[fname]] <-
        lfq_write_table(self$get_long(),
                        path = path_qc,
                        name = fname,
                        format = self$format)
    },
    #' @description
    #' write data to file
    #' @param path_qc path to write to
    write_wide = function(path_qc) {

      wide <- self$get_wide()
      fname <- paste0(self$prefix,"intensities_wide")
      self$file_paths[[fname]] <-
        lfq_write_table(wide$data,
                        path = path_qc,
                        name = fname,
                        format = self$format)

      fname <- paste0(self$prefix,"intensities_file_annotation")
      self$file_paths[[fname]] <-
        lfq_write_table(wide$annotation,
                        path = path_qc,
                        name = fname,
                        format = self$format)
    }
  ))

# LFQDataAggregator ----
#'
#' Decorates LFQData with methods to aggregate protein intensities
#'
#' aggregates intensities
#'
#' @export
#' @family LFQData
#'
#' Aggregate LFQ data
#' @examples
#' istar <- prolfqua_data('data_ionstar')$filtered()
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(data, istar$config)
#'
#' lfqTrans <- lfqdata$clone()$get_Transformer()
#' lfqTrans$log2()
#' lfqTrans <- lfqTrans$robscale()$lfq
#' lfqAggregator <- LFQDataAggregator$new(lfqTrans, "protein")
#'
#' lfqAggregator$medpolish()
#' pmed <- lfqAggregator$plot()
#' pmed$plots[[2]]
#' lfqAggregator$lmrob()
#' prob <- lfqAggregator$plot()
#' prob$plots[[2]]
#'
#' lfqCopy <- lfqdata$clone()
#' lfqCopy$is_transformed()
#' lfqAggregator <- LFQDataAggregator$new(lfqCopy, "protein")
#' lfqAggregator$sum_topN()
#' pSum <- lfqAggregator$plot()
#' pSum$plots[[2]]
#'
#' lfqAggregator$mean_topN()
#' pMean <- lfqAggregator$plot()
#' pMean$plots[[2]]
#' # lfqAggregator$write_plots(".")
#' protPlotter <- lfqAggregator$lfq_agg$get_Plotter()
#' protPlotter$heatmap()
#'
#'
LFQDataAggregator <- R6::R6Class(
  "LFQDataAggregator",
  public = list(
    #' @field lfq LFQData
    lfq = NULL,
    #' @field lfq_agg aggregation result
    lfq_agg = NULL,
    #' @field prefix to use for aggregation results e.g. protein
    prefix = character(),
    ## @field filepath
    ## filepath = character(),
    #' @description
    #' initialize
    #' @param lfq LFQData
    #' @param prefix default protein
    initialize = function(lfq, prefix = "protein"){
      if ( length(lfq$config$table$hierarchyKeys()) == 1 ) {
        stop("no hierarchies to aggregate from: ",  lfq$config$table$hierarchyKeys())
      }
      if (length(lfq$config$table$hierarchyKeys()) == lfq$config$table$hierarchyDepth) {
        stop("no hierarchies to aggregate from: ",
             lfq$config$table$hierarchyKeys(),
             ", hierarchyDepth :",
             lfq$config$table$hierarchyDepth)
      }
      self$lfq = lfq$clone(deep = TRUE)
      self$prefix = prefix
    },
    #' @description
    #' aggregate using median polish
    #' @param N top N by intensity
    #' @return LFQData
    medpolish = function(){
      if (!self$lfq$is_transformed()) {
        warning("You did not transform the intensities.",
                "medpolish works best with already variance stabilized intensities.",
                "Use LFQData$get_Transformer to transform the data.",
                lfq$config$table$workIntensity,)
      }
      res <- aggregate_intensity(self$lfq$data, self$lfq$config, .func = medpolishPlydf_config)
      self$lfq_agg <- LFQData$new(res$data, res$config, prefix = self$prefix)
      invisible(self$lfq_agg)
    },
    #' @description
    #' aggregate using robust regression
    #' @param N top N by intensity
    #' @return LFQData
    lmrob = function(){
      if (!self$lfq$is_transformed()) {
        warning("You did not transform the intensities.",
                "Robust regression works best with already variance stabilized intensities.",
                "Use LFQData$get_Transformer to transform the data.",
                lfq$config$table$workIntensity,)
      }

      res <- aggregate_intensity(self$lfq$data, self$lfq$config, .func = summarizeRobust_config)
      self$lfq_agg <- LFQData$new(res$data, res$config, prefix = self$prefix)
      invisible(self$lfq_agg)
    },
    #' @description
    #' aggregate topN using mean
    #' @param N top N by intensity
    #' @return LFQData
    mean_topN = function(N = 3){
      mean_f <- function(x, name = FALSE){
        if (name) {return("mean")}
        return(mean(x, na.rm = TRUE))
      }
      private$.topN(N = N, .func = mean_f)
    },
    #' @description
    #' aggregate topN using sum
    #' @param N top N by intensity
    #' @return LFQData
    sum_topN = function(N = 3){
      sum_f <- function(x, name = FALSE){
        if (name) { return("sum") }
        sum(x, na.rm = TRUE)}

      private$.topN(N = N, .func = sum_f)
    }
    ,
    #' @description
    #' creates aggreation plots
    #' @param show.legend default FALSE
    #' @return data.frame
    #'
    plot = function(show.legend = FALSE){
      if (is.null(self$lfq_agg)) {
        stop("please aggregate the data first")
      }
      df <- prolfqua::plot_aggregation(
        self$lfq$data,
        self$lfq$config,
        self$lfq_agg$data,
        self$lfq_agg$config,
        show.legend = show.legend)
      invisible(df)
    },
    #' @description
    #' writes plots to folder
    #' @param qcpath qcpath
    #' @param show.legend legend
    #' @param width figure width
    #' @param height figure height
    #' @return file path
    write_plots = function(qcpath, show.legend = FALSE, width = 6, height = 6){
      pl <- self$plot()
      pb <- progress::progress_bar$new(total = nrow(pl))
      filepath <- file.path(qcpath, paste0(self$prefix, "_aggregation_plot.pdf"))
      pdf(filepath , width = width, height = height)
      for (i in seq_len(nrow(pl))) {
        print(pl$plots[[i]])
        pb$tick()
      }
      dev.off()
      self$filepath = filepath
      invisible(filepath)
    }

  ),
  private = list(
    .topN = function(.func, N = 3){
      if (self$lfq$is_transformed()) {
        warning("You did transform the intensities.",
                "top N works with raw data.",
                self$lfq$config$table$workIntensity)
      }

      ranked <- rankPrecursorsByIntensity(self$lfq$data , self$lfq$config)
      resTOPN <- aggregateTopNIntensities(ranked, self$lfq$config, .func = .func, N = N)
      self$lfq_agg <- LFQData$new(resTOPN$data, resTOPN$config, prefix = self$prefix)
      invisible(self$lfq_agg)
    }
  )
)

#' converts LFQData object to SummarizedExperiment
#'
#' For compatibility with Bioconductor
#' @param lfqdata LFQData object
#' @return SummarizedExperiment (bioconductor)
#' @family LFQData
#' @export
#'
LFQDataToSummarizedExperiment <- function(lfqdata){
  if (requireNamespace("SummarizedExperiment")) {
    wide <- lfqdata$to_wide(as.matrix = TRUE)
    ann <- data.frame(wide$annotation)
    rownames(ann) <- wide$annotation$sampleName
    se <- SummarizedExperiment::SummarizedExperiment(S4Vectors::SimpleList(LFQ = wide$data), colData = ann)
    return(se)
  }
}
