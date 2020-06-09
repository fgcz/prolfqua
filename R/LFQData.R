#LFQData ----
#' LFQData R6 class
#' @export
#'
#' @examples
#' #source("c:/Users/wewol/prog/LFQService/R/LFQData.R")
#' istar <- LFQServiceData::dataIonstarFilteredPep
#' istar$data <- istar$data %>% filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' tmp <- lfqdata$to_wide()
#' lfqdata$factors()
#'
#'
LFQData <- R6::R6Class(
  "LFQData",

  public = list(
    config = NULL,
    data = NULL,
    is_pep = FALSE,
    prefix = "",
    initialize = function(data, config, is_pep=TRUE, prefix = "ms_") {
      self$data <- data
      self$config <- config
      self$is_pep <- is_pep
      self$prefix <- prefix
    },

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
    #' getter
    #' @return LFQDataPlotter
    get_Plotter = function(){
      return(LFQDataPlotter$new(self, self$prefix))
    },
    #' @description
    #' getter
    #' @return LFQDataPlotter
    get_Writer = function(format = "xlsx"){
      return(LFQDataWriter$new(self, self$prefix, format = format))
    },
    #' @description
    #' getter
    #' @return LFQDataSummarizer
    get_Summariser = function(){
      return(LFQDataSummariser$new(self, self$prefix))
    },
    #' @description
    #' getter
    #' @return LFQDataStats
    get_Stats = function(){
      return(LFQDataStats$new(self))
    }
  )
)

# LFQDataStats-----
#' compute stdv, mean and CV per peptide or protein and condition.
#' @export
#' @examples
#'
#' # study variance of not normalized data
#' #source("c:/Users/wewol/prog/LFQService/R/LFQData.R")
#' istar <- LFQServiceData::dataIonstarFilteredPep
#' istar$data <- istar$data %>% filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' sum <- lfqdata$get_Stats()
#' sum$cv()
#' sum$cv_quantiles()
#' sum$density()
#' sum$density_median()
#' sum$violin()
#' sum$violin_median()
#' sum$stdv_vs_mean(size = 400)
#' sum$power_t_test()
#' sum$power_t_test_quantiles()
#'
#' #study variance of normalized data
#' istar <- LFQServiceData::dataIonstarNormalizedPep
#' istar$data <- istar$data %>% filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' sum <- lfqdata$get_Stats()
#' sum$cv()
#' sum$cv_quantiles()
#' sum$density()
#' sum$density_median()
#' sum$violin()
#' sum$violin_median()
#' sum$stdv_vs_mean(size = 400)
#' sum$power_t_test()
#' sum$power_t_test_quantiles
#'
#' #Slightly different dataset
#' data <- LFQServiceData::sample_analysis
#' config <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#'
#'
LFQDataStats <- R6::R6Class(
  "LFQDataStats",
  public = list(
    lfq = NULL,
    stat = "CV",
    is_transformed = FALSE,
    initialize = function(lfqdata){
      self$lfq = lfqdata
      self$is_transformed <- self$lfq$config$parameter$is_intensity_transformed
      self$stat <- if (!self$is_transformed) {"CV"}else{"sd"}
    },
    cv = function(){
      LFQService::summarize_cv(self$lfq$data, self$lfq$config)
    },
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
    #' @param ggstat
    #' @return ggplot
    density = function(ggstat = c("density", "ecdf")){
      LFQService::plot_stat_density(
        self$cv(),
        self$lfq$config,
        stat = self$stat,
        ggstat = ggstat)
    },
    density_median = function(ggstat = c("density", "ecdf")){
      LFQService::plot_stat_density_median(self$cv(), self$lfq$config, stat = self$stat)
    },
    violin = function(){
      LFQService::plot_stat_violin(self$cv(), self$lfq$config, stat = self$stat)
    },
    violin_median = function(){
      LFQService::plot_stat_violin_median(self$cv(), self$lfq$config, stat = self$stat)
    },
    stdv_vs_mean = function(size= 200){
      LFQService::plot_stdv_vs_mean(self$cv(), self$lfq$config, size = size)
    },
    power_t_test_quantiles = function(
      probs = c(0.1, 0.25, 0.5, 0.75, 0.9),
      delta = c(0.59,1,2),
      power = 0.8,
      sig.level = 0.05)
    {
      if (!self$is_transformed) {
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
    power_t_test = function(
      delta = c(0.59,1,2),
      power = 0.8,
      sig.level = 0.05
    ){
      if (!self$is_transformed) {
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
#' @examples
#' library(tidyverse)
#'
#' istar <- LFQServiceData::dataIonstarFilteredPep
#' istar$data <- istar$data %>% filter(protein_Id %in% sample(protein_Id, 100))
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
    lfq = NULL,
    prefix = character(),
    initialize = function(lfqdata , prefix = "ms_") {
      self$lfq <- lfqdata
      self$prefix <- prefix
    },
    hierarchy_counts = function(){
      hierarchy_counts(self$lfq$data, self$lfq$config)
    },
    hierarchy_counts_sample = function(value = c("wide", "long", "plot")){
      value <- match.arg(value)
      fun <- LFQService::hierarchy_counts_sample(self$lfq$data, self$lfq$config)
      return(fun(value))
    },
    summarize_hierarchy = function(){
      LFQService::summarize_hierarchy(self$lfq$data, self$lfq$config)
    },
    summarize_protein = function(){
      LFQService::summarize_protein(self$lfq$data, self$lfq$config)
    },
    interaction_missing_stats = function(){
      LFQService::interaction_missing_stats(self$lfq$data, self$lfq$config)
    },
    missingness_per_condition = function(){
      LFQService::missingness_per_condition(self$lfq$data, self$lfq$config)$data
    },
    missingness_per_condition_cumsum = function(){
      LFQService::missingness_per_condition_cumsum(self$lfq$data, self$lfq$config)$data
    }
  )
)

# LFQDataPlotter ----
#' write intensites into folder - for the moment protein
#' @export
#'
#' @examples
#'
#' source("c:/Users/wewol/prog/LFQService/R/LFQData.R")
#' library(LFQService)
#' istar <- LFQServiceData::dataIonstarProtein
#'
#' istar$data <- istar$data %>% filter(protein_Id %in% sample(protein_Id, 100))
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
#' tmp <- lfqplotter$boxplots()
#' tmp$boxplot[[1]]
#' lfqplotter$missigness_histogram()
#' lfqplotter$missingness_per_condition()
#' lfqplotter$missingness_per_condition_cumsum()
#' dev.off()
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
    lfq = NULL,
    #' @field prefix Prepend prefix when writing figures e.g. protein_
    prefix = "",
    #' @field list with paths to figures
    file_paths = list(),
    initialize = function(lfqdata, prefix = "ms_"){
      self$lfq = lfqdata
      self$prefix = prefix
    },
    heatmap = function(na_fraction = 0.3){
      fig <- LFQService::plot_heatmap(self$lfq$data,
                                      self$lfq$config,
                                      na_fraction = na_fraction)
      return(fig)
    },
    heatmap_cor = function(){
      fig <- LFQService::plot_heatmap_cor(self$lfq$data, self$lfq$config)
      return(fig)
    },
    pca = function(){
      fig <- LFQService::plot_pca(self$lfq$data, self$lfq$config, add_txt = FALSE)
      return(fig)
    },
    pca_plotly = function(){
      fig <- plotly::ggplotly(self$pca(), tooltip = self$lfq$config$table$sampleName)
      return(fig)
    },
    boxplots = function(){
      bb <- plot_hierarchies_boxplot_df(self$lfq$data, self$lfq$config)
      return(bb)
    },
    missigness_histogram = function(){
      LFQService::missigness_histogram(self$lfq$data, self$lfq$config)
    },
    missingness_per_condition = function(){
      LFQService::missingness_per_condition(self$lfq$data, self$lfq$config)$figure
    },
    missingness_per_condition_cumsum = function(){
      LFQService::missingness_per_condition_cumsum(self$lfq$data, self$lfq$config)$figure
    },
    NA_heatmap = function(){
      LFQService::plot_NA_heatmap(self$lfq$data, self$lfq$config)
    },
    intensity_distribution_density = function(){
      LFQService::plot_intensity_distribution_density(self$lfq$data, self$lfq$config)
    },
    intensity_distribution_violin = function(){
      LFQService::plot_intensity_distribution_violin(self$lfq$data, self$lfq$config)
    },
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
    },
    sample_correlation = function(){
      LFQService::plot_sample_correlation(self$lfq$data, self$lfq$config)
    },
    write_boxplots = function(path_qc, width = 6, height = 6){
      fpath <- file.path(path_qc,paste0(self$prefix, "boxplot.pdf"))
      message("writing ", fpath)
      bb <- self$boxplots()

      pdf(fpath, width = width, height = height)
      lapply(bb$boxplot, print)
      dev.off()
    },
    write_pltly = function(fig,
                           path_qc,
                           fig_name){
      fname <- paste0(self$prefix, fig_name,".html")
      html_path <- file.path(".", path_qc, fname)
      message("writing ", html_path)
      htmlwidgets::saveWidget(widget = fig,
                              file = fname )
      file.rename(fname, html_path)
      self$file_paths[[fig_name]] <- html_path
      invisible(html_path)
    },
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
      self$file_paths[[fig_name]] <- fpath
      invisible(fpath)
    },
    write = function(path_qc, prefix = "ms"){

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
#'
#' returns long and wide format for writing
#' @export
#'
LFQDataWriter <- R6::R6Class(
  "LFQDataWriter",list(
    lfq = NULL,
    format = "",
    prefix = "",
    file_paths = list(),
    initialize = function(lfq,  prefix = "ms_", format="xlsx"){
      self$lfq = lfq
      self$format = format
      self$prefix = prefix
    },
    get_long = function(){
      #' gets data formatted for writing
      separate_factors(
        separate_hierarchy(self$lfq$data, self$lfq$config),
        self$lfq$config)
    },
    get_wide = function(){
      #' gets data formatted for writing
      wide <- self$lfq$to_wide()
      res <- list(data = separate_hierarchy(wide$data, self$lfq$config), annotation = wide$annotation)
      return(res)
    },
    write_long = function(path_qc) {
      fname <- paste0(self$prefix,"intensities_long")
      self$file_paths[[fname]] <-
        lfq_write_table(self$get_long(),
                        path = path_qc,
                        name = fname,
                        lfq_write_format = self$format)
    },
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

