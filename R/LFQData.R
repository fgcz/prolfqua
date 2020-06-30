#LFQData ----
#'
#' LFQData R6 class
#' @export
#' @family LFQData
#' @examples
#' library(tidyverse)
#' library(LFQService)
#' #source("c:/Users/wewol/prog/LFQService/R/LFQData.R")
#' istar <- LFQServiceData::dataIonstarFilteredPep
#' istar$config$table$is_intensity_transformed
#' istar$data <- istar$data %>% dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' tmp <- lfqdata$to_wide()
#' lfqdata$factors()
#'
#' lfqdata$is_transformed()
#'
#'
LFQData <- R6::R6Class(
  "LFQData",

  public = list(
    #' @field config AnalysisConfiguration
    #' @field data data.frame or tibble matching AnalysisConfiguration.
    #' @field prefix e.g. "peptide_", "protein_", "compound_"
    #' @field is_pep todo
    config = NULL,
    data = NULL,
    is_pep = FALSE,
    prefix = "",
    #' @description
    #' initialize
    #' @param data data.frame
    #' @param config configuration
    #' @param is_pep todo
    #' @param prefix will be use as output prefix
    initialize = function(data, config, is_pep=TRUE, prefix = "ms_") {
      self$data <- data
      self$config <- config
      self$is_pep <- is_pep
      self$prefix <- prefix
    },
    #' @description
    #' get deep copy
    clone_d = function(){
      self$clone(deep = TRUE)
    },
    is_transformed = function(is_transformed){
      if (missing(is_transformed)) {
        return(self$config$table$is_intensity_transformed)
      }else{
        self$config$table$is_intensity_transformed = is_transformed
      }
    },
    #' @description
    #'
    #' Render a QC document
    #'
    #' @param qc_path path to render to
    #' @keywords todo
    render = function(qc_path) {
      LFQService::render_MQSummary_rmd(
        self$data,
        self$config$clone(deep = TRUE),
        pep = self$is_pep,
        workdir = ".",
        dest_path = qc_path,
        dest_file_name = if (self$is_pep) {"peptide_intensities_qc"} else {"protein_intensities_qc"},
        format = "html"
      )
    },
    #' @description
    #' converts the data to wide
    #' @return data and annotation
    to_wide = function(){
      wide <- LFQService::toWideConfig(self$data, self$config)
      wide$config <- self$config
      return(wide)
    },
    #' @description
    #' Annotation table.
    #' @return data.frame
    factors = function(){
      LFQService::table_factors(self$data, self$config)
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
    #' get Stats
    #' @return LFQDataStats
    get_Stats = function(){
      return(LFQDataStats$new(self))
    },
    get_Transformer = function(){
      return(LFQDataTransformer$new(self))
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
      diffdata <- LFQService::filter_difference(self$data,other$data,self$config )
      res <- LFQData$new(diffdata , self$config$clone(deep = TRUE))
      return(res)
    }
  )
)

# LFQDataTransformer ----
#' methods for transforming Intensities
#' @export
#'
#' @examples
#' library(tidyverse)
#' library(LFQService)
#'
#'
#' rm(list = ls())
#' library(LFQService)
#' #source("c:/Users/wewol/prog/LFQService/R/LFQData.R")
#' istar <- LFQServiceData::dataIonstarFilteredPep
#' istar$data <- istar$data %>% dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(istar$data, istar$config)
#'
#' lfqcopy <- lfqdata$clone_d()
#' lfqTrans <- lfqcopy$get_Transformer()
#'
#' x <- lfqTrans$intensity_array(log2)
#' x$lfq$config$table$is_intensity_transformed
#' x <- x$intensity_matrix(robust_scale)
#' plotter <- x$lfq$get_Plotter()
#' plotter$intensity_distribution_density()
#'
#' # transform by asinh root and scale
#' lfqcopy <- lfqdata$clone_d()
#' lfqTrans <- lfqcopy$get_Transformer()
#' x <- lfqTrans$intensity_array(asinh)
#' x$lfq$config$table$is_intensity_transformed
#' x <- lfqTrans$intensity_matrix(robust_scale)
#' plotter <- x$lfq$get_Plotter()
#' plotter$intensity_distribution_density()
#'
LFQDataTransformer <- R6::R6Class(
  "LFQDataTransformer",
  public = list(
    lfq = NULL,

    initialize = function(lfqdata){
      self$lfq = lfqdata
    },
    #' @description
    #' log2 transform and robust scale datas
    #' @return LFQDataTransformer (self)
    log2_robscale = function(){
      r <- LFQService::normalize_log2_robscale(self$lfq$data, self$lfq$config)
      self$lfq$data <- r$data
      self$lfq$config <- r$config
      return(self)
    },
    #' @description
    #' log2 transform and robust scale data based on subset
    #' @param LFQData
    #' @return LFQDataTransformer (self)
    #'
    log2_robscale_subset = function(lfqsubset){
      self$lfq$data  <-  LFQService::transform_work_intensity(self$lfq$data , self$lfq$config, log2)
      self$lfq$data  <-  LFQService::scale_with_subset(self$lfq$data, lfqsubset$data, self$lfq$config)
      self$lfq$is_transformed(TRUE)
      return(self)
    },
    #' @description
    #' Transforms intensities
    #' @param .func transformation function working with arrays e.g. log2, log10, asinh etc.
    #' @return LFQDataTransformer (self)
    #'
    intensity_array = function(.func = log2) {
      .call <- as.list( match.call() )
      r <- LFQService::transform_work_intensity(
        self$lfq$data,
        self$lfq$config,
        .func = .func,
        .funcname = deparse(.call$.func))
      self$lfq$data <- r
      return(self)
    },
    #' @description
    #' @param .func any function taking a matrix and returning a matrix (columns sample, rows feature e.g. base::scale) default robust_scale
    #' @return LFQDataTransformer (self)
    #'
    intensity_matrix = function(.func = robust_scale){
      .call <- as.list( match.call() )
      r <- LFQService::applyToIntensityMatrix(
        self$lfq$data,
        self$lfq$config,
        .func = .func,
        .funcname = deparse(.call$.func))
      self$lfq$data <- r
      return(self)
    }
  )
)

# LFQDataStats-----
#' compute stdv, mean and CV per peptide or protein and condition.
#' @export
#' @family LFQData
#' @examples
#'
#' # study variance of not normalized data
#' #source("c:/Users/wewol/prog/LFQService/R/LFQData.R")
#' runallfuncs <- function(x){
#'
#'   x$cv()
#'   x$cv_quantiles()
#'   x$density()
#'   x$density_median()
#'   x$density("ecdf")
#'   x$density_median("ecdf")
#'   x$violin()
#'   x$violin_median()
#'   x$stdv_vs_mean(size = 400)
#'   x$power_t_test()
#'   x$power_t_test_quantiles()
#' }
#'
#' istar <- LFQServiceData::dataIonstarFilteredPep
#' istar$data <- istar$data %>% dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqstats <- lfqdata$get_Stats()
#' runallfuncs(lfqstats)
#' x<-lfqstats
#'
#' #study variance of normalized data
#' istar <- LFQServiceData::dataIonstarNormalizedPep
#' names(istar)
#' istar$data <- istar$data %>% dplyr::filter(protein_Id %in% sample(protein_Id, 200))
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$is_transformed(TRUE)
#' lfqstats <- lfqdata$get_Stats()
#' runallfuncs(lfqstats)
#'
#' #Slightly different dataset
#' data <- LFQServiceData::sample_analysis
#' config <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#'
#' lfqdata <- LFQData$new(data, config)
#'
#' lfqstats <- lfqdata$get_Stats()
#' runallfuncs(lfqstats)
#'
LFQDataStats <- R6::R6Class(
  "LFQDataStats",
  public = list(
    #' @field lfq LFQData
    #' @field stat either CV or sd (if is_transformed)
    #' @field is_transformed if TRUE data was transformed for stable variance
    lfq = NULL,
    stat = "CV",
    #' @description
    #' create analyse variances and CV
    #' @param lfqdata LFQData object
    initialize = function(lfqdata){
      self$lfq = lfqdata
      self$stat <- if (!self$lfq$is_transformed()) {"CV"}else{"sd"}
    },
    #' @description
    #' compute CV sd and mean of data
    #' @return data.frame
    cv = function(){
      LFQService::summarize_cv(self$lfq$data, self$lfq$config)
    },
    #' @description
    #' Determine CV or sd for the quantiles
    #' @param probs for which quantile to determine CV or sd
    cv_quantiles = function(probs = c(0.1, 0.25, 0.5, 0.75, 0.9)){
      res <- LFQService::summarize_cv_quantiles(
        self$cv(),
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
      LFQService::plot_stat_density(
        self$cv(),
        self$lfq$config,
        stat = self$stat,
        ggstat = ggstat)
    },
    #' @description
    #' plot density or ecdf of CV or sd for the 50% of low intensity data and 50% of high intensity data
    #' @param ggstat either density of ecdf
    #' @return ggplot
    density_median = function(ggstat = c("density", "ecdf")){
      LFQService::plot_stat_density_median(self$cv(), self$lfq$config, stat = self$stat)
    },
    #' @description
    #' plot violinplot of CV or sd
    #' @param ggstat either density of ecdf
    #' @return ggplot
    violin = function(){
      LFQService::plot_stat_violin(self$cv(), self$lfq$config, stat = self$stat)
    },
    #' @description
    #' plot violinplot of CV or sd for the 50% of low intensity data and 50% of high intensity data
    #'
    #' @return ggplot
    #'
    violin_median = function(){
      LFQService::plot_stat_violin_median(self$cv(), self$lfq$config, stat = self$stat)
    },
    #' @description
    #' plot sd vs mean
    #' @param size number of points to sample (default 200)
    #' @return ggplot
    #'
    stdv_vs_mean = function(size= 200){
      LFQService::plot_stdv_vs_mean(self$cv(), self$lfq$config, size = size)
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
      res <- self$cv_quantiles(probs)
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

      res <- LFQService::lfq_power_t_test_proteins(self$cv(),
                                delta = delta,
                                power = power,
                                sig.level = sig.level,
                                min.n = 1.5)
      return(res)
    }
  )
)

# LFQDataSummariser ----
#' generate dataset summaries.
#' @export
#' @family LFQData
#' @examples
#' library(tidyverse)
#'
#' istar <- LFQServiceData::dataIonstarFilteredPep
#' istar$data <- istar$data %>% dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' sum <- lfqdata$get_Summariser()
#' sum$hierarchy_counts()
#' sum$hierarchy_counts_sample("wide")
#' sum$hierarchy_counts_sample("long")
#' sum$hierarchy_counts_sample("plot")
#' sum$summarize_hierarchy()
#' sum$interaction_missing_stats()
#' sum$missingness_per_condition()
#' sum$missingness_per_condition_cumsum()
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
    #' number of elements at each level
    hierarchy_counts = function(){
      LFQService::hierarchy_counts(self$lfq$data, self$lfq$config)
    },
    #' @description
    #' number of elements at each level in every sample
    #' @param value wide - wide format, long - long format, plot - ggplot
    hierarchy_counts_sample = function(value = c("wide", "long", "plot")){
      value <- match.arg(value)
      fun <- LFQService::hierarchy_counts_sample(self$lfq$data, self$lfq$config)
      return(fun(value))
    },
    #' @description
    #' e.g. number of peptides per protein etc
    summarize_hierarchy = function(){
      LFQService::summarize_hierarchy(self$lfq$data, self$lfq$config)
    },
    #' @description
    #' e.g. number of peptides per protein overall
    summarize_protein = function(){
      LFQService::summarize_protein(self$lfq$data, self$lfq$config)
    },
    #' @description
    #' missing per condition and protein
    interaction_missing_stats = function(){
      LFQService::interaction_missing_stats(self$lfq$data, self$lfq$config)
    },
    #' @description
    #' missing stats per condition
    missingness_per_condition = function(){
      LFQService::missingness_per_condition(self$lfq$data, self$lfq$config)$data
    },
    #' @description
    #' missing stats per condition as cumulative sum
    missingness_per_condition_cumsum = function(){
      LFQService::missingness_per_condition_cumsum(self$lfq$data, self$lfq$config)$data
    }
  )
)

# LFQDataPlotter ----
#' Create various visualization of the data
#' @export
#'
#' @family LFQData
#' @import dplyr
#' @examples
#'
#' library(LFQService)
#' istar <- LFQServiceData::dataIonstarProtein
#'
#' istar$data <- istar$data %>% dplyr::filter(protein_Id %in% sample(protein_Id, 100))
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
#' lfqplotter$missingness_per_condition()
#' lfqplotter$missingness_per_condition_cumsum()
#'
#' lfqplotter$NA_heatmap()
#' lfqplotter$intensity_distribution_density()
#' lfqplotter$intensity_distribution_violin()
#' lfqplotter$pairs_smooth()
#' lfqplotter$sample_correlation()
#' LFQService::plot_sample_correlation(istar$data, istar$config)
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
    },
    #' @description
    #' heatmap of intensities
    #' @param na_fraction max fraction of NA's per row
    #' @return pheatmap
    heatmap = function(na_fraction = 0.3){
      fig <- LFQService::plot_heatmap(self$lfq$data,
                                      self$lfq$config,
                                      na_fraction = na_fraction)
      return(fig)
    },
    #' @description
    #' heatmap of sample correlations
    #' @return pheatmap
    heatmap_cor = function(){
      fig <- LFQService::plot_heatmap_cor(self$lfq$data, self$lfq$config)
      return(fig)
    },
    #' @description
    #' pca plot
    #' @return ggplot
    pca = function(){
      fig <- LFQService::plot_pca(self$lfq$data, self$lfq$config, add_txt = FALSE)
      return(fig)
    },
    #' @description
    #' pca plot
    #' @return plotly
    pca_plotly = function(){
      fig <- plotly::ggplotly(self$pca(), tooltip = self$lfq$config$table$sampleName)
      return(fig)
    },
    #' @description
    #' boxplots for all proteins
    #' @return tibble with column boxplots containing ggplot objects
    boxplots = function(){
      bb <- LFQService::plot_hierarchies_boxplot_df(self$lfq$data, self$lfq$config)
      return(bb)
    },
    #' @description
    #' histogram of intensities given number of missing in conditions
    #' @return ggplot
    missigness_histogram = function(){
      LFQService::missigness_histogram(self$lfq$data, self$lfq$config)
    },
    #' @description
    #' barplot with number of features with 1,2, etc missing in condition
    #' @return ggplot
    missingness_per_condition = function(){
      LFQService::missingness_per_condition(self$lfq$data, self$lfq$config)$figure
    },
    #' @description
    #' barplot with cumulative sum of features with 1,2, etc missing in condition
    #' @return ggplot
    missingness_per_condition_cumsum = function(){
      LFQService::missingness_per_condition_cumsum(self$lfq$data, self$lfq$config)$figure
    },
    #' @description
    #' heatmap of features with missing values
    #' @return ggplot
    NA_heatmap = function(){
      LFQService::plot_NA_heatmap(self$lfq$data, self$lfq$config)
    },
    #' @description
    #' density distribution of intensities
    #' @return ggplot
    intensity_distribution_density = function(){
      LFQService::plot_intensity_distribution_density(self$lfq$data, self$lfq$config)
    },
    #' @description
    #' Violinplot showing distribution of intensities in all samples
    #' @return ggplot
    intensity_distribution_violin = function(){
      LFQService::plot_intensity_distribution_violin(self$lfq$data, self$lfq$config)
    },
    #' @description
    #' pairsplot of intensities
    #' @param max maximal number of samples to show
    #' @return NULL
    pairs_smooth = function(max=10){
      dataTransformed <- self$lfq$data
      config <- self$lfq$config
      samples <- dplyr::select(self$lfq$data, config$table$sampleName) %>%
        distinct() %>%
        pull()
      if (length(samples) > max) {
        limit <- samples %>% sample(max)
        ldata <- dataTransformed %>%
          dplyr::filter(!!sym(config$table$sampleName) %in% limit)
        LFQService::pairs_smooth( LFQService::toWideConfig(ldata, config, as.matrix = TRUE)$data )
      }else{
        LFQService::pairs_smooth( LFQService::toWideConfig(dataTransformed, config, as.matrix = TRUE)$data )
      }
      NULL
    },
    #' @description
    #' plot of sample correlations
    #' @return NULL
    sample_correlation = function(){
      LFQService::plot_sample_correlation(self$lfq$data, self$lfq$config)
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
      pdf(fpath,
          width = width,
          height = height)
      #if ('pheatmap' %in% class(fig)) {
      #  print(plot(fig$gtable))
      #}else{
        print(fig)
      #}
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
#' Write LFQ data, or provide outputs for writing.
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
                        lfq_write_format = self$format)
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
                        lfq_write_format = self$format)

      fname <- paste0(self$prefix,"intensities_file_annotation")
      self$file_paths[[fname]] <-
      lfq_write_table(wide$annotation,
                      path = path_qc,
                      name = fname,
                      lfq_write_format = self$format)
    }
  ))

# LFQDataAggregator ----
#' LFQAggregator
#'
#' aggregates intensities
#'
#' @export
#' @family LFQData
#'
#' Aggregate LFQ data
#' @examples
#' library(tidyverse)
#' library(LFQService)
#'
#'
#' rm(list = ls())
#' istar <- LFQServiceData::dataIonstarFilteredPep
#' istar$data <- istar$data %>% dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(istar$data, istar$config)
#'
#' lfqTrans <- lfqdata$clone_d()$get_Transformer()$log2_robscale()
#' lfqAggregator <- LFQDataAggregator$new(lfqTrans$lfq, "protein")
#'
#' lfqAggregator$medpolish()
#' pmed <- lfqAggregator$plot()
#' pmed$plots[[2]]
#' lfqAggregator$lmrob()
#' prob <- lfqAggregator$plot()
#' prob$plots[[2]]
#'
#' lfqCopy <- lfqdata$clone_d()
#' lfqCopy$is_transformed()
#' lfqAggregator <- LFQDataAggregator$new(lfqCopy, "protein")
#' lfqAggregator$sum_topN()
#' pSum <- lfqAggregator$plot()
#' pSum$plots[[2]]
#'
#' lfqAggregator$mean_topN()
#' pMean <- lfqAggregator$plot()
#' pMean$plots[[2]]
#' lfqAggregator$write_plots("inst")
#' protPlotter <- lfqAggregator$lfq_agg$get_Plotter()
#'
LFQDataAggregator <- R6::R6Class(
  "LFQDataAggregator",
  public = list(
    #' @field lfq LFQData
    #' @field lfq_agg aggregation result
    #' @field prefix to use for aggregation results e.g. protein
    lfq = NULL,
    lfq_agg = NULL,
    prefix = character(),
    filepath = character(),
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
      self$lfq = lfq
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
      df <- LFQService::plot_aggregation(
        self$lfq$data,
        self$lfq$config,
        self$lfq_agg$data,
        self$lfq_agg$config,
        show.legend = show.legend)
      invisible(df)
    },
    #' @description
    #' Writes plots to folder
    #' @param qcpath qcpath
    #' @param legend legend
    #' @param width figure width
    #' @param height figure height
    #' @return file path
    write_plots = function(qcpath, show.legend = FALSE, width = 6, height = 6){
      pl <- self$plot()
      pb <- progress::progress_bar$new(total = nrow(pl))
      filepath <- file.path(qcpath, paste0(self$prefix, "_aggregation_plot.pdf"))
      pdf(filepath , width = width, height = height)
      for (i in 1:nrow(pl)) {
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
