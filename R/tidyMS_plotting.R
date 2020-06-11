#' visualize intensity distributions
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#'
#' @export
#'
#' @keywords internal
#' @import ggplot2
#' @family plotting
#' @examples
#' config <- LFQServiceData::skylineconfig$clone(deep=TRUE)
#' plot_intensity_distribution_violin(LFQServiceData::sample_analysis, config)
#' analysis <- transform_work_intensity(LFQServiceData::sample_analysis, config, log2)
#' plot_intensity_distribution_violin(analysis, config)
plot_intensity_distribution_violin <- function(pdata, config){
  p <- ggplot(pdata, aes_string(x = config$table$sampleName, y = config$table$getWorkIntensity() )) +
    geom_violin() +
    theme(axis.text.x = element_text(angle = ))
  if (!config$parameter$is_intensity_transformed) {
    p <- p + scale_y_continuous(trans = 'log10')
  }
  return(p)
}

#' visualize intensity distributions
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#'
#' @export
#' @keywords internal
#' @import ggplot2
#' @family plotting
#' @rdname plot_intensity_distribution_violin
#' @examples
#'
#' config <- LFQServiceData::skylineconfig$clone(deep=TRUE)
#' plot_intensity_distribution_density(LFQServiceData::sample_analysis, config)
#' analysis <- transform_work_intensity(LFQServiceData::sample_analysis, config, log2)
#' plot_intensity_distribution_density(analysis, config)
plot_intensity_distribution_density <- function(pdata, config){
  p <- ggplot(pdata, aes_string(x = config$table$getWorkIntensity(), colour = config$table$sampleName )) +
    geom_line(stat = "density")
  if (!config$parameter$is_intensity_transformed) {
    p <- p + scale_x_continuous(trans = 'log10')
  }
  return(p)
}

#' visualize correlation among samples
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#'
#' @export
#' @keywords internal
#' @family plotting
#' @rdname plot_sample_correlation
#' @examples
#' config <- LFQServiceData::skylineconfig$clone(deep=TRUE)
#' analysis <- remove_small_intensities(LFQServiceData::sample_analysis, config)
#' analysis <- transform_work_intensity(analysis, config, log2)
#' mm <- toWideConfig(analysis, config, as.matrix = TRUE)
#' plot_sample_correlation(analysis, config)
plot_sample_correlation <- function(pdata, config){
  matrix <- toWideConfig(pdata, config, as.matrix = TRUE)$data
  M <- cor(matrix, use = "pairwise.complete.obs")
  if (nrow(M) > 12) {
    corrplot::corrplot.mixed(M,upper = "ellipse",
                             lower = "pie",
                             diag = "u",
                             tl.cex = .6,
                             tl.pos = "lt",
                             tl.col = "black",
                             mar = c(2,5,5,2))
  } else{
    corrplot::corrplot.mixed(M,upper = "ellipse",
                             lower = "number",
                             lower.col = "black",
                             tl.cex = .6,
                             number.cex = .7,
                             diag = "u",
                             tl.pos = "lt",
                             tl.col = "black",
                             mar = c(2,5,5,2))

  }
  invisible(M)
}

#' plot peptides by factors and it's levels.
#'
#' @param pdata data.frame
#' @param title name to show
#' @param config AnalysisConfiguration
#' @param facet_grid_on on which variable to run facet_grid
#' @param beeswarm use beeswarm default FALSE
#'
#' @export
#' @keywords internal
#' @examples
#'
#' library(LFQService)
#' library(tidyverse)
#'
#' conf <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#' conf$table$hierarchyDepth
#' conf$table$hkeysDepth()
#'
#' xnested <- LFQServiceData::sample_analysis %>%
#'  group_by_at(conf$table$hkeysDepth()) %>% tidyr::nest()
#'
#'
#' p <- plot_hierarchies_boxplot(xnested$data[[3]],
#'  xnested$protein_Id[[3]],
#'   conf,
#'   facet_grid_on =  tail(conf$table$hierarchyKeys(),1))
#' p
#'
#' p <- plot_hierarchies_boxplot(xnested$data[[3]],
#'  xnested$protein_Id[[3]],
#'   conf )
#' p
#'
#' plot_hierarchies_boxplot(
#'  xnested$data[[3]],
#'  xnested$protein_Id[[3]],
#'   conf,
#'    beeswarm = FALSE )
#'
#'
#' config <- LFQServiceData::skylineconfig_HL$clone(deep = TRUE)
#' data <- LFQServiceData::skylineSRM_HL_data
#' data <- setup_analysis(data,config)
#' data <- LFQService::transform_work_intensity(data, config, log2)
#' res <- plot_hierarchies_boxplot_df(data, config)
#' res$boxplot[[1]]
#'
#' hierarchy = config$table$hkeysDepth()
#' xnested <- data %>% dplyr::group_by_at(hierarchy) %>% tidyr::nest()
#' p <- plot_hierarchies_boxplot(xnested$data[[1]], xnested$protein_Id[[1]],config, beeswarm = FALSE)
#' p
#' p <- plot_hierarchies_boxplot(xnested$data[[1]], xnested$protein_Id[[1]],config, beeswarm = TRUE)
#' p
#' p <- plot_hierarchies_boxplot(xnested$data[[1]], xnested$protein_Id[[1]],config, beeswarm = TRUE, facet_grid_on = "precursor_Id")
#' p
#'
plot_hierarchies_boxplot <- function(pdata,
                                     title,
                                     config,
                                     facet_grid_on = NULL ,
                                     beeswarm = TRUE){

  isotopeLabel <- config$table$isotopeLabel
  lil <- length(unique(pdata[[isotopeLabel]]))

  pdata <- LFQService::make_interaction_column( pdata , c(config$table$fkeysDepth()))
  color <- if (lil > 1) {isotopeLabel} else {NULL}
  p <- ggplot(pdata, aes_string(x = "interaction",
                              y = config$table$getWorkIntensity(),
                              color = color
  )) + theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    ggtitle(title)

  if (!config$parameter$is_intensity_transformed) {
    p <- p + scale_y_continuous(trans = "log10")
  }
  p <- p + geom_boxplot()
  if ( beeswarm ) {
    #p <- p + geom_point()
    p <- p + ggbeeswarm::geom_quasirandom(aes_string( color = color) , dodge.width = 0.7 )
  }
  if (!is.null( facet_grid_on ) && (facet_grid_on %in% colnames(pdata))) {
    p <- p + facet_grid( formula(paste0("~", facet_grid_on ) ))
  }
  p
  return(p)
}


#' generates peptide level plots for all Proteins
#' @export
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param hiearchy e.g. protein_Id default hkeysDepth
#' @param facet_grid_on default NULL
#'
#' @keywords internal
#' @examples
#' resDataStart <- LFQServiceData::testDataStart2954$resDataStart
#' config <-  LFQServiceData::testDataStart2954$config
#' res <- plot_hierarchies_boxplot_df(resDataStart, config)
#' res$boxplot[[1]]
#'
#' res <- plot_hierarchies_boxplot_df(resDataStart,
#'  config,
#'  facet_grid_on = "peptide_Id")
#' res$boxplot[[1]]
#' config <- config$clone(deep = TRUE)
#' #TODO plot on peptide level.
#' config$table$hierarchyDepth = 2
#'
#' res <- plot_hierarchies_boxplot_df(resDataStart, config)
#'
#' res$boxplot[[1]]
#'
plot_hierarchies_boxplot_df <- function(pdata,
                                        config,
                                        hierarchy = config$table$hkeysDepth(),
                                        facet_grid_on = NULL){

  xnested <- pdata %>% dplyr::group_by_at(hierarchy) %>% tidyr::nest()
  newcol <- paste(hierarchy, collapse = "+")
  xnested <- xnested %>% tidyr::unite(!!sym(newcol), one_of(hierarchy))

  figs <- xnested %>%
    dplyr::mutate(boxplot = map2(data, !!sym(newcol),
                                 plot_hierarchies_boxplot,
                                 config ,
                                 facet_grid_on = facet_grid_on))
  return(dplyr::select_at(figs, c(newcol, "boxplot")))
}
