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
#'
#' bb <- LFQService::ionstar$filtered()
#' config <- bb$config$clone(deep=TRUE)
#' analysis <- bb$data
#'
#' plot_intensity_distribution_violin(analysis, config)
#' analysis <- transform_work_intensity(analysis, config, log2)
#' plot_intensity_distribution_violin(analysis, config)
#'
plot_intensity_distribution_violin <- function(pdata, config){
  p <- ggplot(pdata, aes_string(x = config$table$sampleName, y = config$table$getWorkIntensity() )) +
    geom_violin() +
    theme(axis.text.x = element_text(angle = ))
  if (!config$table$is_intensity_transformed) {
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
#' @family plotting
#' @rdname plot_intensity_distribution_violin
#' @examples
#'
#' bb <- LFQService::ionstar$filtered()
#' config <- bb$config$clone(deep=TRUE)
#' analysis <- bb$data
#' plot_intensity_distribution_density(analysis, config)
#' analysis <- transform_work_intensity(analysis, config, log2)
#' plot_intensity_distribution_density(analysis, config)
#'
plot_intensity_distribution_density <- function(pdata, config){
  p <- ggplot(pdata, aes_string(x = config$table$getWorkIntensity(), colour = config$table$sampleName )) +
    geom_line(stat = "density")
  if (!config$table$is_intensity_transformed) {
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
#'
#' bb <- LFQService::ionstar$filtered()
#' config <- bb$config$clone(deep=TRUE)
#' analysis <- bb$data
#'
#' analysis <- remove_small_intensities(analysis, config)
#' analysis <- transform_work_intensity(analysis, config, log2)
#' mm <- toWideConfig(analysis, config, as.matrix = TRUE)
#' class(plot_sample_correlation(analysis, config))
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
#' bb <- LFQService::ionstar$filtered()
#' conf <- bb$config$clone(deep=TRUE)
#' analysis <- bb$data
#' conf$table$hierarchyDepth
#' conf$table$hkeysDepth()
#'
#' xnested <- analysis %>%
#'  group_by_at(conf$table$hkeysDepth()) %>% tidyr::nest()
#'
#' #debug(plot_hierarchies_boxplot)
#' p <- plot_hierarchies_boxplot(xnested$data[[3]],
#'  xnested$protein_Id[[3]],
#'   conf,
#'   facet_grid_on =  tail(conf$table$hierarchyKeys(),1))
#' stopifnot("ggplot" %in% class(p))
#'
#' p <- plot_hierarchies_boxplot(xnested$data[[3]],
#'  xnested$protein_Id[[3]],
#'   conf )
#' stopifnot("ggplot" %in% class(p))
#'
#' p <- plot_hierarchies_boxplot(
#'  xnested$data[[3]],
#'  xnested$protein_Id[[3]],
#'   conf,
#'    beeswarm = FALSE )
#' stopifnot("ggplot" %in% class(p))
#'
#' bb <- LFQServiceData::skylineSRM_HL_A
#' config <- bb$config_f()
#' analysis <- bb$analysis(bb$data, config)
#' data <- LFQService::transform_work_intensity(analysis, config, log2)
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
                                     beeswarm = TRUE,
                                     pb){
  if (!missing(pb)) { pb$tick() }

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

  if (!config$table$is_intensity_transformed) {
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
#'
#'  iostar <- LFQService::ionstar$filtered()
#'  iostar$data <- iostar$data %>%
#'    dplyr::filter(protein_Id %in% sample(protein_Id, 2))
#'  unique(iostar$data$protein_Id)
#'
#'  res <- plot_hierarchies_boxplot_df(iostar$data,iostar$config)
#'  res$boxplot[[1]]
#'  res <- plot_hierarchies_boxplot_df(iostar$data,iostar$config,iostar$config$table$hierarchyKeys()[1])
#'  res$boxplot[[1]]
#'  res <- plot_hierarchies_boxplot_df(iostar$data,iostar$config,
#'                                     iostar$config$table$hierarchyKeys()[1],
#'                                     facet_grid_on = iostar$config$table$hierarchyKeys()[2])
#'  res$boxplot[[1]]
#'
#'  iostar <- LFQServiceData::dataIonstarProtein
#'  iostar$data <- iostar$data %>%
#'    dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#'  unique(iostar$data$protein_Id)
#'
#'  res <- plot_hierarchies_boxplot_df(iostar$data,iostar$config)
#'  res$boxplot[[1]]
#'  res <- plot_hierarchies_boxplot_df(iostar$data,iostar$config,
#'                                     iostar$config$table$hierarchyKeys()[1])
#'  res$boxplot[[1]]
#'  res <- plot_hierarchies_boxplot_df(iostar$data,iostar$config,
#'                                     iostar$config$table$hierarchyKeys()[1],
#'                                     facet_grid_on = iostar$config$table$hierarchyKeys()[2])
#'  res$boxplot[[1]]
plot_hierarchies_boxplot_df <- function(pdata,
                                        config,
                                        hierarchy = config$table$hkeysDepth(),
                                        facet_grid_on = NULL){

  xnested <- pdata %>% dplyr::group_by_at(hierarchy) %>% tidyr::nest()
  newcol <- paste(hierarchy, collapse = "+")
  xnested <- xnested %>% tidyr::unite(!!sym(newcol), one_of(hierarchy))

  pb <- progress::progress_bar$new(total = nrow(xnested))

  figs <- xnested %>%
    dplyr::mutate(boxplot = map2(data, !!sym(newcol),
                                 plot_hierarchies_boxplot,
                                 config ,
                                 facet_grid_on = facet_grid_on, pb = pb))
  return(dplyr::select_at(figs, c(newcol, "boxplot")))
}



#' plot correlation heatmap with annotations
#'
#' @export
#' @keywords internal
#' @examples
#' library(tidyverse)
#' bb <- LFQService::ionstar$filtered()
#' config <- bb$config$clone(deep=TRUE)
#' analysis <- bb$data
#'
#' pheat_map <- LFQService::plot_heatmap_cor( analysis, config )
#' stopifnot("pheatmap" %in% class(pheat_map))
#' pheat_map <- plot_heatmap_cor( analysis, config, R2 = TRUE )
#' stopifnot("pheatmap" %in% class(pheat_map))
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
#' @examples
#' library(tidyverse)
#' bb <- LFQService::ionstar$filtered()
#' config <- bb$config$clone(deep=TRUE)
#' analysis <- bb$data
#'
#'
#'
#' graphics.off()
#' .Device
#'  p  <- plot_heatmap(analysis, config)
#' .Device
#'  print(p)
#'  .Device
#'  plot(1)
#'
#'  print(p)
#'

plot_heatmap <- function(data, config, na_fraction = 0.4, ...){
  res <-  toWideConfig(data, config , as.matrix = TRUE)
  annot <- res$annotation
  resdata <- res$data

  factors <- dplyr::select_at(annot, config$table$factorKeys())
  factors <- as.data.frame(factors)
  rownames(factors) <- annot$sampleName

  resdata <- t(scale(t(resdata)))
  resdata <- LFQService::removeNArows(resdata,floor(ncol(resdata)*na_fraction))


  # not showing row dendrogram trick
  res <- pheatmap::pheatmap(resdata,
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
#' @examples
#'
#' library(tidyverse)
#' library(LFQService)
#' bb <- LFQService::ionstar$filtered()
#' config <- bb$config$clone(deep=TRUE)
#' analysis <- bb$data
#'
#' tmp <- plot_NA_heatmap(analysis, config)
#' print(tmp)
#' xx <- plot_NA_heatmap(analysis, config, distance = "euclidean")
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
#'
#' bb <- LFQService::ionstar$filtered()
#' config <- bb$config$clone(deep=TRUE)
#' analysis <- bb$data
#'
#' tmp <- plot_pca(analysis, config, add_txt= TRUE)
#'
#' print(tmp)
#' tmp <- plot_pca(analysis, config, add_txt= FALSE)
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
    point +
    if (add_txt) {text}

  if (!is.na(sh)) {
    x <- x +  ggplot2::scale_shape_manual(values = 1:length(unique(xx[[sh]])))
  }
  return(x)
}

