#' LFQData R6 class
#' @export
#'
#' @examples
#'
#' data <- LFQService::data_c
#' config <- config <- LFQService::config_c
#' lfqdata <- LFQData$new(data, config)
#' tmp <- lfqdata$to_wide()
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
    get_Plotter = function(){
      return(LFQDataPlotter$new(self, self$prefix))
    },
    get_Writer = function(format = "xlsx"){
      return(LFQDataWriter$new(self, self$prefix, format = format))
    }
  )

)


#' write intensites into folder - for the moment protein
#' @export
#'
#' @examples
#' library(LFQService)
#' data <- LFQService::dataIonstarProtein
#' LFQService::plot_pca(data$data, data$config, add_txt = FALSE)
#'
#' dim(data$data)
#' dataR <- data$data %>% filter(protein_Id %in% sample(data$data$protein_Id,100))
#'
#' lfqdata <- LFQData$new(dataR, data$config)
#' lfqplotter <- LFQDataPlotter$new(lfqdata)
#' graphics.off()
#' lfqplotter$plot_heatmap()
#' graphics.off()
#' lfqplotter$plot_heatmap_cor()
#' graphics.off()
#' lfqplotter$plot_pca()
#' tmp <- lfqplotter$plot_prot_boxplots()
#' tmp$boxplot[[1]]
LFQDataPlotter <- R6::R6Class(
  "LFQDataPlotter",
  list(
    lfq = NULL,
    prefix = "",
    initialize = function(lfqdata, prefix = "ms_"){
      self$lfq = lfqdata
      self$prefix = prefix
    },
    plot_heatmap = function(na_fraction = 0.3){
      fig <- LFQService::plot_heatmap(self$lfq$data,
                                      self$lfq$config,
                                      na_fraction = na_fraction)
      return(fig)
    },
    plot_heatmap_cor = function(){
      fig <- LFQService::plot_heatmap_cor(self$lfq$data, self$lfq$config)
      return(fig)
    },
    plot_pca = function(){
      fig <- LFQService::plot_pca(self$lfq$data, self$lfq$config, add_txt = FALSE)
      return(fig)
    },
    plot_pca_plotly = function(path_qc){
      fig <- plotly::ggplotly(self$plot_pca(), tooltip = self$config$table$sampleName)
      return(fig)
    },
    plot_prot_boxplots = function(){
      bb <- plot_hierarchies_boxplot_df(self$lfq$data, self$lfq$config)
      return(bb)
    },
    write_prot_boxplots = function(path_qc, width = 6, height = 6){
      fpath <- file.path(path_qc,paste0(self$prefix, "boxplot.pdf"))
      message("writing ", fpath)
      pdf(fpath, width = width, height = height)
      bb <- self$plot_prot_boxplots()
      lapply(bb$boxplot, print)
      dev.off()
    },
    write_pltly = function(fig,
                           path_qc,
                           fig_name,
                           prefix = "ms_"){
      fname <- paste0(self$prefix, fig_name,".html")
      html_path <- file.path(".", path_qc, fname)
      message("writing ", html_path)
      htmlwidgets::saveWidget(widget = fig,
                              file = fname )
      file.rename(fname, html_path)
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
    },
    write = function(path_qc, prefix = "ms"){
      self$write_pdf(self$plot_heatmap_cor(),
                     path_qc,
                     "intensities_heatmap_correlation",
                     width = 10, height = 10)
      self$write_pdf(self$plot_heatmap(),
                     path_qc,
                     "intensities_heatmap",
                     width = 10, height = 10)
      self$write_pdf(self$plot_pca(),
                     path_qc,
                     "intensities_PCA")
      self$write_pltly(self$plot_pca_plotly(),
                       path_qc,
                       "intensities_PCA")
    }
  ))


#' LFQData write data decorator
#'
#' returns long and wide format for writing
#' @export
LFQDataWriter <- R6::R6Class(
  "LFQDataWriter",list(
    lfq = NULL,
    format = "",
    prefix = "",

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
    write_long = function(path_qc){
      path <- file.path(path_qc,paste0(self$prefix,"intensities_long",".csv"))
      message("writing ",path)
      lfq_write_table(self$get_long(),
                      path = path, lfq_write_format = self$format)
    },
    write_wide = function(path_qc){

      wide <- self$get_wide()
      path <- file.path(path_qc,paste0(self$prefix,"intensities_wide",".csv"))
      message("writing ",path)
      lfq_write_table(wide$data,
                      path = path, lfq_write_format = self$format)
      path <- file.path(path_qc,paste0(self$prefix,"intensities_file_annotation",".csv"))
      message("writing ",path)
      lfq_write_table(wide$annotation,
                      path = path, lfq_write_format = self$format)
    }
  ))

