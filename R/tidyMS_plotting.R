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
#' istar <-sim_lfq_data_peptide_config()
#'
#' config <- istar$config
#' analysis <- istar$data
#'
#' plot_intensity_distribution_violin(analysis, config)
#' analysis <- transform_work_intensity(analysis, config, log2)
#' plot_intensity_distribution_violin(analysis, config)
#'
plot_intensity_distribution_violin <- function(pdata, config){
  p <- ggplot(pdata, aes_string(x = config$table$sampleName, y = config$table$get_response() )) +
    geom_violin() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    stat_summary(fun.y = median, geom = "point", size = 1, color = "black")
  if (!config$table$is_response_transformed) {
    p <- p + scale_y_continuous(trans = 'log10')
  }
  return(p)
}

#' visualize intensity distributions
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param legend do not show legend
#' @export
#' @keywords internal
#' @family plotting
#' @rdname plot_intensity_distribution_violin
#' @examples
#'
#' istar <-sim_lfq_data_peptide_config()
#'
#' config <- istar$config
#' analysis <- istar$data
#' plot_intensity_distribution_density(analysis, config)
#' analysis <- transform_work_intensity(analysis, config, log2)
#' plot_intensity_distribution_density(analysis, config)
#'
plot_intensity_distribution_density <- function(pdata, config, legend = TRUE){
  p <- ggplot(pdata, aes_string(x = config$table$get_response(),
                                colour = config$table$sampleName )) +
    geom_line(stat = "density")
  if (!config$table$is_response_transformed) {
    p <- p + scale_x_continuous(trans = 'log10')
  }
  if (!legend) {
    p <- p + ggplot2::guides(colour = FALSE)
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
#' istar <-sim_lfq_data_peptide_config()
#'
#' config <- istar$config
#' analysis <- istar$data
#'
#' analysis <- remove_small_intensities(analysis, config)
#' analysis <- transform_work_intensity(analysis, config, log2)
#' mm <- tidy_to_wide_config(analysis, config, as.matrix = TRUE)
#' class(plot_sample_correlation(analysis, config))
#' plot_sample_correlation(analysis, config)
plot_sample_correlation <- function(pdata, config){
  matrix <- tidy_to_wide_config(pdata, config, as.matrix = TRUE)$data
  M <- cor(matrix, use = "pairwise.complete.obs")
  if (nrow(M) > 12) {
    res <- corrplot::corrplot.mixed(M,upper = "ellipse",
                             lower = "pie",
                             diag = "u",
                             tl.cex = .6,
                             tl.pos = "lt",
                             tl.col = "black",
                             mar = c(2,5,5,2))
  } else{
    res <- corrplot::corrplot.mixed(M,upper = "ellipse",
                             lower = "number",
                             lower.col = "black",
                             tl.cex = .6,
                             number.cex = .7,
                             diag = "u",
                             tl.pos = "lt",
                             tl.col = "black",
                             mar = c(2,5,5,2))

  }
  invisible(res)
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
#'
#' bb <- prolfqua_data('data_skylineSRM_HL_A')
#' config <- bb$config_f()
#' analysis <- bb$analysis(bb$data, config)
#' data <- prolfqua::transform_work_intensity(analysis, config, log2)
#' res <- plot_hierarchies_boxplot_df(data, config)
#' res$boxplot[[1]]
#'
#' hierarchy = config$table$hierarchy_keys_depth()
#' xnested <- data |> dplyr::group_by_at(hierarchy) |> tidyr::nest()
#' p <- plot_hierarchies_boxplot(xnested$data[[1]], xnested$protein_Id[[1]],config, beeswarm = FALSE)
#' p <- plot_hierarchies_boxplot(xnested$data[[1]], xnested$protein_Id[[1]],config, beeswarm = TRUE)
#' p <- plot_hierarchies_boxplot(xnested$data[[1]], xnested$protein_Id[[1]],config, beeswarm = TRUE, facet_grid_on = "precursor_Id")
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

  pdata <- prolfqua::make_interaction_column( pdata , c(config$table$factor_keys_depth()))
  color <- if (lil > 1) {isotopeLabel} else {NULL}
  p <- ggplot(pdata, aes_string(x = "interaction",
                                y = config$table$get_response(),
                                color = color
  )) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    ggtitle(title)

  if (!config$table$is_response_transformed) {
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
#' @param hiearchy e.g. protein_Id default hierarchy_keys_depth
#' @param facet_grid_on default NULL
#' @family plotting
#' @keywords internal
#' @examples
#'
#'
#'  istar <- sim_lfq_data_peptide_config()
#'  config <- istar$config
#'  analysis <- istar$data
#'  analysis <- analysis |>
#'    dplyr::filter(protein_Id %in% sample(protein_Id, 2))
#'
#'  res <- plot_hierarchies_boxplot_df(analysis,config)
#'  res$boxplot[[1]]
#'  res <- plot_hierarchies_boxplot_df(analysis,config,config$table$hierarchy_keys()[1])
#'  res$boxplot[[1]]
#'  res <- plot_hierarchies_boxplot_df(analysis,config,
#'                                     config$table$hierarchy_keys()[1],
#'                                     facet_grid_on = config$table$hierarchy_keys()[2])
#'  res$boxplot[[1]]
#'  res$boxplot[[2]]
#'
#'  iostar <- sim_lfq_data_protein_config()
#'  iostar$data <- iostar$data |>
#'    dplyr::filter(protein_Id %in% sample(protein_Id, 4))
#'  unique(iostar$data$protein_Id)
#'
#'  res <- plot_hierarchies_boxplot_df(iostar$data,iostar$config)
#'  res$boxplot[[1]]
#'  res <- plot_hierarchies_boxplot_df(iostar$data,iostar$config,
#'                                     iostar$config$table$hierarchy_keys()[1])
plot_hierarchies_boxplot_df <- function(pdata,
                                        config,
                                        hierarchy = config$table$hierarchy_keys_depth(),
                                        facet_grid_on = NULL){

  xnested <- pdata |> dplyr::group_by_at(hierarchy) |> tidyr::nest()
  newcol <- paste(hierarchy, collapse = "+")
  xnested <- xnested |> tidyr::unite(!!sym(newcol), one_of(hierarchy))

  pb <- progress::progress_bar$new(total = nrow(xnested))

  figs <- xnested |>
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
#' @family plotting
#' @examples
#'
#' istar <- sim_lfq_data_protein_config()
#' config <- istar$config
#' analysis <- istar$data
#'
#' pheat_map <- prolfqua::plot_heatmap_cor( analysis, config )
#' stopifnot("pheatmap" %in% class(pheat_map))
#' pheat_map <- plot_heatmap_cor( analysis, config, R2 = TRUE )
#' stopifnot("pheatmap" %in% class(pheat_map))
#'
plot_heatmap_cor <- function(data,
                             config,
                             R2 = FALSE,
                             color = colorRampPalette(c("white", "red"))(1024),
                             ...){

  res <-  tidy_to_wide_config(data, config , as.matrix = TRUE)
  annot <- res$annotation
  res <- res$data


  cres <- cor(res, use = "pa")
  if (R2) {
    cres <- cres^2
  }

  factors <- dplyr::select_at(annot, config$table$factor_keys())
  factors <- as.data.frame(factors)
  rownames(factors) <- annot[[config$table$sampleName]]

  #res <- pheatmap::pheatmap(cres,
  #                          scale = "none",
  #                          silent = TRUE)

  gg <- stats::hclust(stats::dist(cres))
  res <- pheatmap::pheatmap(cres[gg$order,],
                            scale = "none",
                            cluster_rows  = FALSE,
                            annotation_col = factors,
                            show_rownames = FALSE,
                            border_color = NA,
                            main = ifelse(R2, "R^2", "correlation"),
                            silent = TRUE,
                            color = color,
                            ... = ...)
  invisible(res)
}

.ehandler = function(e){
  warning("WARN :", e)
  # return string here
  as.character(e)
}

#' plot heatmap with annotations
#'
#' @export
#' @param na_fraction fraction of NA values per row
#' @param show_rownames if TRUE shows row names, default FALSE
#' @keywords internal
#' @family plotting
#' @examples
#'
#' istar <- sim_lfq_data_protein_config()
#' config <- istar$config
#' analysis <- istar$data
#'
#'
#' p  <- plot_heatmap(analysis, config)
#' stopifnot(class(p) == "pheatmap")
#' p2 <- plot_heatmap(analysis, config, show_rownames = TRUE)
#' stopifnot(class(p) == "pheatmap")
#'
plot_heatmap <- function(data,
                         config,
                         na_fraction = 0.4,
                         show_rownames = FALSE , ...){
  if (nrow(data) == 0 ) {
    warning("The dataset has :", nrow(data), "")
    return(NULL)
  }

  wide <-  tidy_to_wide_config(data, config , as.matrix = TRUE)
  annot <- wide$annotation

  factors <- dplyr::select_at(annot, config$table$factor_keys())
  factors <- as.data.frame(factors)
  rownames(factors) <- annot[[config$table$sampleName]]
  resdata <- t(scale(t(wide$data)))
  resdataf <- prolfqua::remove_NA_rows(resdata,floor(ncol(resdata)*na_fraction))

  if (nrow(resdataf) >= 3) {
    gg <- stats::hclust( stats::dist( resdataf ))
    res <- pheatmap::pheatmap(resdataf[gg$order,],
                              cluster_rows  = FALSE,
                              scale = "row",
                              annotation_col = factors,
                              show_rownames = show_rownames,
                              border_color = NA,
                              silent = TRUE,
                              ... = ...)


  } else {

    res <- tryCatch(pheatmap::pheatmap(resdata,
                                       cluster_rows  = FALSE,
                                       scale = "row",
                                       annotation_col = factors,
                                       show_rownames = show_rownames,
                                       border_color = NA,
                                       silent = TRUE,
                                       ... = ...), error = .ehandler)


  }
  invisible(res)
}

#' plot heatmap without any clustering (use to show NA's)
#' @param data dataframe
#' @param config dataframe configuration
#' @param arrange either mean or var
#' @param not_na if true than arrange by nr of NA's first and then by arrange
#' @param y.labels show y labels
#' @keywords internal
#'
#' @family plotting
#' @export
#' @examples
#'
#' istar <- sim_lfq_data_protein_config()
#' config <- istar$config
#' analysis <- istar$data
#' rs <- plot_raster(analysis, config, show_rownames=FALSE)
#' stopifnot(class(rs) == "pheatmap")
#' rs <- plot_raster(analysis[1,], config)
#' stopifnot(is.null(rs))
#' rs <- plot_raster(analysis, config, "var")
#' stopifnot(class(rs) == "pheatmap")
#' rs <- plot_raster(analysis, config, show_rownames = TRUE)
#' stopifnot(class(rs) == "pheatmap")
#'
plot_raster <- function(data,
                        config,
                        arrange = c("mean", "var"),
                        not_na = FALSE,
                        show_rownames = FALSE,
                        ...) {
  if (nrow(data) <= 1 ) {
    warning("The dataset has :", nrow(data), "")
    return(NULL)
  }
  arrange <- match.arg(arrange)
  res <-  tidy_to_wide_config(data, config , as.matrix = TRUE)
  annot <- res$annotation
  resdata <- res$data


  factors <- dplyr::select_at(annot, config$table$factor_keys())
  factors <- as.data.frame(factors)
  rownames(factors) <- annot[[config$table$sampleName]]


  if (arrange == "mean") {
    bb <- apply(resdata, 1, mean, na.rm = TRUE)
  } else if (arrange == "var") {
    bb <- apply(resdata, 1, stats::var, na.rm = TRUE)
  }
  if (not_na) {
    bNA <- apply(resdata, 1, function(x){sum(is.na(x))})
    resdata <- resdata[order(bNA, bb, decreasing = c(FALSE, TRUE)), , drop = FALSE]
  } else {
    resdata <- resdata[order(bb, decreasing =  TRUE), , drop = FALSE]
  }

  res <- pheatmap::pheatmap(resdata,
                            cluster_rows  = FALSE,
                            cluster_cols = FALSE,
                            annotation_col = factors,
                            show_rownames = show_rownames,
                            border_color = NA,
                            silent = TRUE,
                            ... = ...)

  invisible(res)
}





#' plot heatmap of NA values
#' @export
#' @keywords internal
#' @family plotting
#' @examples
#'
#'
#' istar <- sim_lfq_data_peptide_config()
#' config <- istar$config
#' analysis <- istar$data
#'
#'
#' tmp <- plot_NA_heatmap(analysis, config)
#' stopifnot(class(tmp) == "pheatmap")
#' tmp <- plot_NA_heatmap(analysis, config, distance = "euclidean")
#' stopifnot(class(tmp) == "pheatmap")
#'
#'
plot_NA_heatmap <- function(data,
                            config,
                            limitrows = 10000,
                            distance = "binary"){
  res <-  tidy_to_wide_config(data, config , as.matrix = TRUE)
  annot <- res$annotation
  res <- res$data
  stopifnot(annot[[config$table$sampleName]] %in% colnames(res))

  factors <- dplyr::select_at(annot, config$table$factor_keys())
  factors <- as.data.frame(factors)
  rownames(factors) <- annot[[config$table$sampleName]]

  res[!is.na(res)] <- 0
  res[is.na(res)] <- 1
  allrows <- nrow(res)
  res <- res[apply(res,1, sum) > 0, , drop = FALSE]

  message("rows with NA's: ", nrow(res), "; all rows :", allrows, "\n")

  if (nrow(res) > 1) {
    res <- if (nrow(res) > limitrows ) {
      message("limiting nr of rows to:", limitrows,"\n")
      res[sample( seq_len(nrow(res)),limitrows),]
    }else{
      res
    }

    gg <- stats::hclust( stats::dist( res, method = distance ))
    resclust <- pheatmap::pheatmap(res[gg$order,],
                                   cluster_rows  = FALSE,
                                   clustering_distance_cols = distance,
                                   scale = "none",
                                   annotation_col = factors,
                                   color = c("white", "black"),
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
#' @family plotting
#' @examples
#'
#'
#'
#'
#' istar <- sim_lfq_data_protein_config(with_missing = FALSE)
#' config <- istar$config
#' analysis <- istar$data
#' tmp <- plot_pca(analysis, config, add_txt= TRUE)
#' stopifnot("ggplot" %in% class(tmp) )
#' tmp <- plot_pca(analysis, config, add_txt= FALSE)
#' stopifnot("ggplot" %in% class(tmp) )
#' tmp
#' tmp <- plot_pca(analysis, config, PC = c(1,2))
#' stopifnot("ggplot" %in% class(tmp) )
#' tmp <- plot_pca(analysis, config, PC = c(2,40))
#' stopifnot(is.null(tmp))
#'
plot_pca <- function(data , config, PC = c(1,2), add_txt = FALSE, plotly = FALSE){
  stopifnot(length(PC) == 2)

  wide <- tidy_to_wide_config(data, config ,as.matrix = TRUE)
  ff <- na.omit(wide$data)
  ff <- t(ff)
  pca_result <- prcomp(ff)
  xx <- as_tibble(pca_result$x, rownames = config$table$sampleName)
  if (max(PC) > ncol(xx)) {
    warning("nr of PCs: ", ncol(xx), "\n")
    return(NULL)
    }
  variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
  xx <- inner_join(wide$annotation, xx)


  sh <- config$table$factor_keys()[2]
  point <- (if (!is.na(sh)) {
    geom_point(aes(shape = !!sym(sh)))
  }else{
    geom_point()
  })

  text <- geom_text(aes(label = !!sym(config$table$sampleName)),check_overlap = TRUE,
                    nudge_x = 0.25,
                    nudge_y = 0.25 )

    PCx <- paste0("PC", PC[1])
    PCy <- paste0("PC", PC[2])
    x <- ggplot(xx, aes(x = !!sym(PCx), y = !!sym(PCy),
                        color = !!sym(config$table$factor_keys()[1]),
                        text = !!sym(config$table$sampleName))) +
      labs(x = paste0(PCx," (", round(variance_explained[PC[1]]), "% variance)"),
           y = paste0(PCy," (", round(variance_explained[PC[2]]), "% variance)")) +
      point +
      if (add_txt) {text}
  if (!is.na(sh)) {
    x <- x +  ggplot2::scale_shape_manual(values = seq_along(unique(xx[[sh]])))
  }
  return(x)
}

#' plot screeplot
#' @export
#' @keywords internal
#' @family plotting
#' @examples
#' istar <- sim_lfq_data_protein_config(with_missing = FALSE)
#' config <- istar$config
#' analysis <- istar$data
#' tmp <- plot_screeplot(analysis, config, threshold_pc= NULL)
#' print(tmp)
#' tmp <- plot_screeplot(analysis, config, threshold_pc= 1)
#' print(tmp)
#' tmp <- plot_screeplot(analysis, config, nr_PC = 4, threshold_pc = NULL)
#' print(tmp)
plot_screeplot <- function(data , config, threshold_pc = 1, nr_PC = NULL) {
  wide <- tidy_to_wide_config(data, config ,as.matrix = TRUE)
  ff <- na.omit(wide$data)
  ff <- t(ff)
  pca_result <- prcomp(ff)
  variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
  xx <- data.frame(PC = paste("PC", 1:length(variance_explained), sep = "_"), percent_variance_explained = variance_explained)
  xx$PC <- factor(xx$PC, levels = xx$PC)

  if (!is.null(threshold_pc)) {
    xx <- xx |> dplyr::filter(percent_variance_explained > threshold_pc)
  }
  if (!is.null(nr_PC)) {
    minRow <- min(nr_PC, nrow(xx))
    xx <- xx[1:minRow,]
  }
  nudgeval = -2
  pl <- ggplot2::ggplot(xx, ggplot2::aes(x = PC, y = percent_variance_explained )) +
    ggplot2::geom_bar(stat = "identity", position = "dodge", colour = "black", fill = "white") +
    ggplot2::geom_text(aes(label = round(.data$percent_variance_explained)), nudge_y = nudgeval, angle = 65) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

  invisible(pl)
}
