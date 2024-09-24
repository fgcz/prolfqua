#' LFQDataPlotter ----
#' Create various visualization of the LFQdata
#' @export
#'
#' @family LFQData
#' @import dplyr
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#'
#' lfqdata <- LFQData$new(
#'  istar$data,
#'  istar$config)
#' #LFQDataPlotter$debug("boxplots")
#' LFQDataPlotter$debug("pairs_smooth")
#' lfqplotter <- lfqdata$get_Plotter()
#'
#' stopifnot(class(lfqplotter$heatmap()) == "pheatmap")
#' stopifnot(class(lfqplotter$heatmap_cor()) == "pheatmap")
#' stopifnot("ggplot" %in% class(lfqplotter$pca()))
#' stopifnot("plotly" %in%  class(lfqplotter$pca_plotly()))
#' tmp <- lfqplotter$boxplots()
#' stopifnot("ggplot" %in%  class(tmp$boxplot[[1]]))
#' stopifnot("ggplot" %in% class(lfqplotter$missigness_histogram()))
#'
#' stopifnot(class(lfqplotter$NA_heatmap()) == "pheatmap")
#' class(lfqplotter$intensity_distribution_density())
#' class(lfqplotter$intensity_distribution_violin())
#' \dontrun{
#'   stopifnot(is.null(lfqplotter$pairs_smooth()))
#' }
#' stopifnot(class(lfqplotter$sample_correlation()) == "list")
#' stopifnot(class(lfqplotter$raster()) == "pheatmap")
#' stopifnot("upset" == class(lfqplotter$upset_missing()))
#' stopifnot(class(prolfqua::plot_sample_correlation(istar$data, istar$config)) == "list")
#'
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
    #'
    #' heatmap of intensities - columns are samples, rows are proteins or peptides.
    #'
    #' The abundances of each protein (row) are z-scored.
    #' Afterward, the mean abundance for each protein is zero,
    #' and the standard variation is one.
    #' z-scoring allows to compare (cluster) the proteins according
    #' to the difference in the expression in the samples.
    #' Without the z-scoring, the proteins would group according
    #' to their abundance, e.g., high abundant proteins would be one cluster.
    #'
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
    #' heatmap of sample correlations.
    #'
    #' The Spearman correlation among all samples
    #' is computed. Then the euclidean distance is used to compute the distances.
    #'
    #' @seealso \code{\link{plot_heatmap_cor}}
    #'
    #' @return pheatmap
    #'
    heatmap_cor = function(){
      fig <- prolfqua::plot_heatmap_cor(self$lfq$data, self$lfq$config)
      return(fig)
    },
    #' @description
    #' PCA plot
    #'
    #' A PCA is applied and the first and second principal component are shown.
    #'
    #' @seealso \code{\link{plot_pca}}
    #'
    #' @param add_txt show sample names
    #' @param PC default c(1,2) - first and second principal component
    #' @return ggplot
    pca = function(PC = c(1,2), add_txt = TRUE){
      fig <- prolfqua::plot_pca(self$lfq$data, self$lfq$config, PC = PC,add_txt = add_txt)
      return(fig)
    },
    #' @description
    #' pca plot
    #' @param add_txt show sample names
    #' @param PC default c(1,2) - first and second principal component
    #' @return plotly
    pca_plotly = function(PC = c(1,2), add_txt = FALSE){
      fig <- plotly::ggplotly(self$pca(add_txt = add_txt), tooltip = self$lfq$config$table$sampleName)
      return(fig)
    },
    #' @description
    #' boxplots for all proteins
    #' @param facet enable facet wrap if hierarchy_depth less then hierarchy lenght.
    #' @return tibble with column boxplots containing ggplot objects
    boxplots = function(facet = TRUE){
      config <- self$lfq$config
      if (config$table$hierarchyDepth < length(config$table$hierarchy) && facet) {
        bb <- prolfqua::plot_hierarchies_boxplot_df(
          self$lfq$data, config,
          hierarchy = config$table$hierarchy_keys_depth(),
          facet_grid_on = config$table$hierarchy_keys()[config$table$hierarchyDepth + 1])

      } else {
        bb <- prolfqua::plot_hierarchies_boxplot_df(self$lfq$data, config,
                                                    hierarchy = config$table$hierarchy_keys_depth(),
                                                    facet_grid_on = NULL)

      }
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
        prolfqua::pairs_smooth( prolfqua::tidy_to_wide_config(ldata, config, as.matrix = TRUE)$data )
      }else{
        prolfqua::pairs_smooth( prolfqua::tidy_to_wide_config(dataTransformed, config, as.matrix = TRUE)$data )
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
    #' upset plot based on presence absence information
    #' @return plot
    upset_missing = function(){
      pups <- prolfqua::UpSet_missing_stats(self$lfq$data, self$lfq$config)
      res <- UpSetR::upset(pups$data , order.by = "freq", nsets = pups$nsets)
      return(res)
    },
    #' @description
    #' write boxplots to file
    #' @param path_qc path to write to
    #' @param filename file to write into
    #' @param width fig width
    #' @param height fig height
    #'
    write_boxplots = function(path_qc, filename = NULL, width = 6, height = 6){
      if (is.null(filename)) {
        filename = self$prefix
      }
      fpath <- file.path(path_qc,paste0(filename, "boxplot.pdf"))
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
