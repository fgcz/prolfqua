#' LFQData R6 class
#' @export
#' @examples
#'
#' data <- LFQService::data_c
#' config <- config <- LFQService::config_c
#' lfqdata <- LFQData$new(data, config)
#' lfqdata$to_wide()
LFQData <- R6::R6Class(
  "LFQData",

  public = list(
    config = NULL,
    data = NULL,
    is_pep = FALSE,
    initialize = function(data, config, is_pep=TRUE) {
      self$data <- data
      self$config <- config
      self$is_pep <- is_pep
    },

    render = function(qc_path) {
      LFQService::render_MQSummary_rmd(
        self$data,
        self$config$clone(deep = TRUE),
        pep = self$is_pep,
        workdir = ".",
        dest_path = qc_path,
        dest_file_name = if (is_pep) {"peptide_intensities_qc"} else {"protein_intensities_qc"},
        format = "html"
      )

    },
    to_wide = function(){
      #' to wide
      wide <- LFQService::toWideConfig(self$data, self$config)
      wide$config <- self$config
      return(wide)
    }


  )

)


#' write intensites into folder - for the moment protein
#' @export
#'
#' @examples
#'
#' data <- LFQService::data_c
#' config <- config <- LFQService::config_c
#' debug(LFQService::plot_pca)
#' LFQService::plot_pca(data, config, add_txt = FALSE)
#'
#' lfqdata <- LFQData$new(data, config)
#' lfqplotter <- LFQDataPlotter$new(lfqdata)
#' lfqplotter$plot_heatmap()
#' lfqplotter$plot_heatmap_cor()
#' lfqplotter$plot_pca()
LFQDataPlotter <- R6::R6Class(
  "LFQDataWrite",
  list(
    data = NULL,
    config = NULL,
    initialize = function(lfqdata){
      self$data = lfqdata$data
      self$config = lfqdata$config
    },
    plot_heatmap = function(na_fraction = 0.3){
      fig <- LFQService::plot_heatmap(self$data, self$config, na_fraction = na_fraction)
      return(fig)
    },
    plot_heatmap_cor = function(){
      fig <- LFQService::plot_heatmap_cor(self$data, self$config)
      return(fig)
    },
    plot_pca = function(){
      fig <- LFQService::plot_pca(self$data, self$config, add_txt = FALSE)
      return(fig)
    },
    plot_pca_plotly = function(path_qc){
      fig <- plotly::ggplotly(self$plot_pca(), tooltip = self$config$table$sampleName)
      return(fig)
    },
    write_pltly = function(fig, fig_name, prefix = "ms_"){
      fname <- paste0(prefix, fig_name,".html")
      html_path <- file.path(".",path_qc,fname)
      htmlwidgets::saveWidget(widget = fig,
                              file = paste0(self$prefix,fig_name,self$suffix,".html") )
      file.rename(fname, html_path)
    },
    write_pdf = function(fig, path_qc, fig_name, prefix = "ms_" , width=7 , height=7 ){
      pdf(file.path(path_qc,paste0(prefix,fig_name,".pdf")),
          width = width,
          height = height)
      print(fig)
      dev.off()
    },
    write = function(path_qc, prefix = "ms"){
      write_pdf(self$plot_heatmap(), "intensities_heatmap_", width = 10, height = 10)
      write_pdf(self$plot_heatmap(), "intensities_heatmap_correlation_", width = 10, height = 10)
      write_pdf(self$plot_pca(), "intensities_PCA_")
      write_pltly(self$plot_pca_plotly(),"intensities_PCA_")
    }
  ))


#' LFQData write data decorator
#'
#' returns long and wide format for writing
#' @export
LFQDataWriter <- R6::R6Class(
  "LFQDataWriter",list(
    lfqdata = NULL,
    path_qc = "",
    format = "",
    prefix = "",

    initialize = function(lfqdata, path_qc, prefix = "ms_", format="xlsx"){
      self$lfqdata = lfqdata$data
      self$path_qc = path_qc
      self$format = format
      self$prefix = prefix
    },
    get_long = function(){
      #' gets data formatted for writing
      separate_factors(
        separate_hierarchy(self$lfqdata$data, self$lfqdata$config),
        self$lfqdata$config)
    },
    get_wide = function(){
      #' gets data formatted for writing
      wide <- self$lfqdata$to_wide()
      res <- list(data = separate_hierarchy(wide$data, self$lfqdata$config), annotation = wide$annotation)
      return(res)
    },
    write_long = function(){
      path <- file.path(self$path_qc,paste0(self$prefix,"intensities_long",".csv"))
      lfq_write_table(self$get_long(),
                      path = path, lfq_write_format = self$format)
    },
    write_wide = function(){
      wide <- self$get_wide()
      path <- file.path(self$path_qc,paste0(prefix,"intensities_wide",".csv"))
      lfq_write_table(wide$data,
                      path = path, lfq_write_format = self$format)
      path <- file.path(path_qc,paste0(prefix,"intensities_file_annotation",".csv"))
      lfq_write_table(wide$annotation,
                      path = path, lfq_write_format = self$format)
    }
  ))

