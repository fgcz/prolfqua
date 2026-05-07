#' Print method for pheatmap objects
#'
#' Fixes pheatmap's print method which omits \code{grid.newpage()},
#' causing the heatmap to draw on top of the previous plot when
#' multiple plots are produced in the same knitr chunk or graphics device.
#'
#' @param x a pheatmap object
#' @param ... ignored
#' @return The computed result.
#' @export
#' @method print pheatmap
print.pheatmap <- function(x, ...) {
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
}

.show_sample_legend <- function(pdata, sample_name, legend, max_legend_samples) {
  if (is.logical(legend) && length(legend) == 1 && !is.na(legend)) {
    return(legend)
  }
  n_samples <- length(unique(stats::na.omit(pdata[[sample_name]])))
  n_samples <= max_legend_samples
}

.truncate_plot_labels <- function(labels, max_chars = 60) {
  if (is.null(max_chars) || !is.finite(max_chars)) {
    return(as.character(labels))
  }
  max_chars <- as.integer(max_chars)
  if (max_chars < 4) {
    stop("max_chars must be at least 4.", call. = FALSE)
  }

  labels <- as.character(labels)
  label_width <- nchar(labels, type = "chars", allowNA = FALSE)
  too_long <- label_width > max_chars
  if (!any(too_long)) {
    return(labels)
  }

  keep_chars <- max_chars - 3
  left_chars <- ceiling(keep_chars / 2)
  right_chars <- floor(keep_chars / 2)
  labels[too_long] <- paste0(
    substr(labels[too_long], 1, left_chars),
    "...",
    substr(
      labels[too_long],
      label_width[too_long] - right_chars + 1,
      label_width[too_long]
    )
  )
  labels
}

#' visualize intensity distributions
#' @param pdata data.frame
#' @param sample_name character — sample column name
#' @param response character — intensity column name
#' @param is_transformed logical — is intensity log-transformed?
#'
#' @export
#'
#' @keywords internal
#' @import ggplot2
#' @family plotting
#' @examples
#' istar <- sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' plot_intensity_distribution_violin(
#'   lfq$data_long(), lfq$sample_name(), lfq$response(), lfq$is_transformed())
#'
plot_intensity_distribution_violin <- function(pdata, sample_name, response, is_transformed = FALSE) {
  p <- ggplot(pdata, aes(x = .data[[sample_name]], y = .data[[response]])) +
    geom_violin() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    stat_summary(fun = median, geom = "point", size = 1, color = "black")
  if (!is_transformed) {
    p <- p + scale_y_continuous(trans = "log10")
  }
  return(p)
}

#' visualize intensity distributions
#' @param pdata data.frame
#' @param sample_name character — sample column name
#' @param response character — intensity column name
#' @param is_transformed logical — is intensity log-transformed?
#' @param legend show legend. If `NA`, hide automatically when the number of
#'   samples is larger than `max_legend_samples`.
#' @param max_legend_samples maximum number of samples for automatic legend
#'   display.
#' @export
#' @keywords internal
#' @family plotting
#' @rdname plot_intensity_distribution_violin
#' @examples
#' istar <- sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' plot_intensity_distribution_density(
#'   lfq$data_long(), lfq$sample_name(), lfq$response(), lfq$is_transformed())
#'
plot_intensity_distribution_density <- function(
  pdata,
  sample_name,
  response,
  is_transformed = FALSE,
  legend = NA,
  max_legend_samples = 16
) {
  p <- ggplot(pdata, aes(x = .data[[response]], colour = .data[[sample_name]])) +
    geom_line(stat = "density")
  if (!is_transformed) {
    p <- p + scale_x_continuous(trans = "log10")
  }
  if (!.show_sample_legend(pdata, sample_name, legend, max_legend_samples)) {
    p <- p + ggplot2::guides(colour = "none")
  }
  return(p)
}

#' visualize correlation among samples
#' @param matrix numeric matrix — wide-format intensity data (samples as columns)
#'
#' @export
#' @keywords internal
#' @family plotting
#' @rdname plot_sample_correlation
#' @examples
#' istar <- sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' lfq <- lfq$get_Transformer()$log2()$lfq
#' plot_sample_correlation(lfq$data_wide(as.matrix = TRUE)$data)
plot_sample_correlation <- function(matrix) {
  M <- cor(matrix, use = "pairwise.complete.obs")
  if (nrow(M) > 12) {
    res <- corrplot::corrplot.mixed(
      M,
      upper = "ellipse",
      lower = "pie",
      diag = "u",
      tl.cex = .6,
      tl.pos = "lt",
      tl.col = "black",
      mar = c(2, 5, 5, 2)
    )
  } else {
    res <- corrplot::corrplot.mixed(
      M,
      upper = "ellipse",
      lower = "number",
      lower.col = "black",
      tl.cex = .6,
      number.cex = .7,
      diag = "u",
      tl.pos = "lt",
      tl.col = "black",
      mar = c(2, 5, 5, 2)
    )
  }
  invisible(res)
}

#' plot peptides by factors and it's levels.
#'
#' @param pdata data.frame
#' @param title name to show
#' @param lfqdata LFQData object
#' @param facet_grid_on on which variable to run facet_grid
#' @param beeswarm use beeswarm default FALSE
#' @param show_mean show mean values
#' @param pb progress bar
#'
#' @export
#' @keywords internal
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' res <- plot_hierarchies_boxplot_df(lfq$data_long(), lfq)
#' res$boxplot[[1]]
#'
#' xnested <- lfq$data_long() |>
#'   dplyr::group_by(across(all_of(lfq$relevant_hierarchy_keys()))) |> tidyr::nest()
#' p <- plot_hierarchies_boxplot(xnested$data[[1]], xnested$protein_Id[[1]],
#'   lfq, beeswarm = FALSE, show_mean = TRUE)
#' p <- plot_hierarchies_boxplot(xnested$data[[1]], xnested$protein_Id[[1]],
#'   lfq, beeswarm = TRUE)
#' p
plot_hierarchies_boxplot <- function(
  pdata,
  title,
  lfqdata,
  facet_grid_on = NULL,
  beeswarm = TRUE,
  show_mean = TRUE,
  pb
) {
  if (!missing(pb)) {
    pb$tick()
  }

  isotope_col <- lfqdata$isotope_label()
  nr_children_col <- lfqdata$nr_children_col()
  response <- lfqdata$response()

  lil <- length(unique(pdata[[isotope_col]]))

  pdata <- prolfqua::make_interaction_column(pdata, c(lfqdata$relevant_factor_keys()))
  pdata$size <- ifelse(pdata[[nr_children_col]] == 0, 2, pdata[[nr_children_col]])
  pdata[[nr_children_col]] <- as.factor(pdata[[nr_children_col]])
  color <- if (lil > 1) {
    isotope_col
  } else {
    NULL
  }
  p <- ggplot(pdata, aes(x = .data[["interaction"]], y = .data[[response]])) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    ggtitle(title)
  if (!is.null(color)) {
    p <- p + aes(colour = .data[[color]])
  }

  if (!lfqdata$is_transformed()) {
    p <- p + scale_y_continuous(trans = "log10")
  }
  p <- p + geom_boxplot()
  if (beeswarm) {
    if (length(levels(pdata[[nr_children_col]])) > 1) {
      shape_values <- c(4, rep(16, length(levels(pdata[[nr_children_col]])) - 1))
    } else {
      shape_values <- 16
    }
    p <- p +
      ggbeeswarm::geom_quasirandom(
        aes(
          size = .data[["size"]],
          shape = .data[[nr_children_col]]
        ),
        dodge.width = 0.7
      ) +
      scale_shape_manual(values = shape_values) +
      scale_size_continuous(range = range(pdata$size, na.rm = TRUE))
  }
  if (show_mean) {
    p <- p + stat_summary(fun = mean, geom = "point", position = position_dodge(0.7), size = 3, shape = 3)
    p <- p +
      stat_summary(
        fun = mean,
        geom = "text",
        aes(label = round(after_stat(y), 2)),
        position = position_dodge(0.7),
        vjust = -1,
        size = 3
      )
  }
  if (!is.null(facet_grid_on) && (facet_grid_on %in% colnames(pdata))) {
    p <- p + facet_grid(formula(paste0("~", facet_grid_on)))
  }
  return(p)
}


#' generates peptide level plots for all Proteins
#' @export
#' @param pdata data.frame
#' @param lfqdata LFQData object
#' @param hierarchy e.g. protein_Id default relevant_hierarchy_keys
#' @param facet_grid_on default NULL
#' @family plotting
#' @keywords internal
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' res <- plot_hierarchies_boxplot_df(lfq$data_long(), lfq)
#' res$boxplot[[1]]
#'
#' lfq2 <- LFQData$new(
#'   istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 2)),
#'   istar$config)
#' res <- plot_hierarchies_boxplot_df(lfq2$data_long(), lfq2)
#' res$boxplot[[1]]
plot_hierarchies_boxplot_df <- function(
  pdata,
  lfqdata,
  hierarchy = lfqdata$relevant_hierarchy_keys(),
  facet_grid_on = NULL
) {
  xnested <- pdata |> dplyr::group_by(across(all_of(hierarchy))) |> tidyr::nest()
  newcol <- paste(hierarchy, collapse = "+")
  xnested <- xnested |> tidyr::unite(!!sym(newcol), all_of(hierarchy))

  pb <- progress::progress_bar$new(total = nrow(xnested))

  figs <- xnested |>
    dplyr::mutate(
      boxplot = map2(
        data,
        !!sym(newcol),
        plot_hierarchies_boxplot,
        lfqdata = lfqdata,
        facet_grid_on = facet_grid_on,
        pb = pb
      )
    )
  return(dplyr::select(figs, all_of(c(newcol, "boxplot"))))
}


#' plot correlation heatmap with annotations
#'
#' @param matrix numeric matrix — wide-format intensity data
#' @param annotation data.frame — sample annotation
#' @param factor_keys character vector — factor column names for annotation
#' @param sample_name character — sample name column
#' @param R2 logical — plot R-squared instead of correlation
#' @param color color palette
#' @param ... passed to pheatmap
#' @export
#' @keywords internal
#' @family plotting
#' @examples
#'
#' istar <- sim_lfq_data_protein_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' wide <- lfq$data_wide(as.matrix = TRUE)
#' pheat_map <- plot_heatmap_cor(wide$data, wide$annotation,
#'   lfq$factor_keys(), lfq$sample_name())
#' stopifnot("pheatmap" %in% class(pheat_map))
#' pheat_map <- plot_heatmap_cor(wide$data, wide$annotation,
#'   lfq$factor_keys(), lfq$sample_name(), R2 = TRUE)
#' stopifnot("pheatmap" %in% class(pheat_map))
#'
plot_heatmap_cor <- function(
  matrix,
  annotation,
  factor_keys,
  sample_name,
  R2 = FALSE,
  color = colorRampPalette(c("white", "red"))(1024),
  ...
) {
  cres <- cor(matrix, use = "pa")
  if (R2) {
    cres <- cres^2
  }

  factors <- dplyr::select(annotation, all_of(factor_keys))
  factors <- as.data.frame(factors)
  rownames(factors) <- annotation[[sample_name]]

  gg <- stats::hclust(stats::dist(cres))
  res <- pheatmap::pheatmap(
    cres[gg$order, ],
    scale = "none",
    cluster_rows = FALSE,
    annotation_col = factors,
    show_rownames = FALSE,
    border_color = NA,
    main = ifelse(R2, "R^2", "correlation"),
    silent = TRUE,
    color = color,
    ... = ...
  )
  invisible(res)
}


#' plot heatmap with annotations
#'
#' @param matrix numeric matrix — wide-format intensity data
#' @param annotation data.frame — sample annotation
#' @param factor_keys character vector — factor column names for annotation
#' @param sample_name character — sample name column
#' @param na_fraction fraction of NA values per row
#' @param show_rownames if TRUE shows row names, default FALSE
#' @param max_rownames_chars maximum displayed row label length
#' @param ... passed to pheatmap
#' @export
#' @keywords internal
#' @family plotting
#' @examples
#'
#' istar <- sim_lfq_data_protein_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' wide <- lfq$data_wide(as.matrix = TRUE)
#' p <- plot_heatmap(wide$data, wide$annotation, lfq$factor_keys(), lfq$sample_name())
#' stopifnot(class(p) == "pheatmap")
#'
plot_heatmap <- function(
  matrix,
  annotation,
  factor_keys,
  sample_name,
  na_fraction = 0.4,
  show_rownames = FALSE,
  max_rownames_chars = 60,
  ...
) {
  if (nrow(matrix) == 0) {
    warning("The dataset has :", nrow(matrix), "")
    return(NULL)
  }

  factors <- dplyr::select(annotation, all_of(factor_keys))
  factors <- as.data.frame(factors)
  rownames(factors) <- annotation[[sample_name]]
  resdata <- t(scale(t(matrix)))
  resdataf <- prolfqua::remove_na_rows(resdata, floor(ncol(resdata) * na_fraction))

  if (nrow(resdataf) >= 3) {
    gg <- stats::hclust(stats::dist(resdataf))
    plot_data <- resdataf[gg$order, ]
    res <- pheatmap::pheatmap(
      plot_data,
      cluster_rows = FALSE,
      scale = "row",
      annotation_col = factors,
      show_rownames = show_rownames,
      labels_row = .truncate_plot_labels(rownames(plot_data), max_rownames_chars),
      border_color = NA,
      silent = TRUE,
      ... = ...
    )
  } else {
    res <- tryCatch(
      pheatmap::pheatmap(
        resdata,
        cluster_rows = FALSE,
        scale = "row",
        annotation_col = factors,
        show_rownames = show_rownames,
        labels_row = .truncate_plot_labels(rownames(resdata), max_rownames_chars),
        border_color = NA,
        silent = TRUE,
        ... = ...
      ),
      error = .error_handler # nolint object_usage_linter. defined in utilities.R
    )
  }
  invisible(res)
}


#' plot heatmap without any clustering (use to show NA's)
#' @param matrix numeric matrix — wide-format intensity data
#' @param annotation data.frame — sample annotation
#' @param factor_keys character vector — factor column names for annotation
#' @param sample_name character — sample name column
#' @param arrange either mean or var
#' @param not_na if true than arrange by nr of NA's first and then by arrange
#' @param show_rownames logical, show row names in heatmap
#' @param max_rownames_chars maximum displayed row label length
#' @param ... additional arguments passed to pheatmap
#' @keywords internal
#'
#' @family plotting
#' @export
#' @examples
#'
#' istar <- sim_lfq_data_protein_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' wide <- lfq$data_wide(as.matrix = TRUE)
#' rs <- plot_raster(wide$data, wide$annotation, lfq$factor_keys(), lfq$sample_name())
#' stopifnot(class(rs) == "pheatmap")
#' rs <- plot_raster(wide$data, wide$annotation, lfq$factor_keys(), lfq$sample_name(), "var")
#' stopifnot(class(rs) == "pheatmap")
#'
plot_raster <- function(
  matrix,
  annotation,
  factor_keys,
  sample_name,
  arrange = c("mean", "var"),
  not_na = FALSE,
  show_rownames = FALSE,
  max_rownames_chars = 60,
  ...
) {
  if (nrow(matrix) <= 1) {
    warning("The dataset has :", nrow(matrix), "")
    return(NULL)
  }
  arrange <- match.arg(arrange)

  factors <- dplyr::select(annotation, all_of(factor_keys))
  factors <- as.data.frame(factors)
  rownames(factors) <- annotation[[sample_name]]

  if (arrange == "mean") {
    bb <- apply(matrix, 1, mean, na.rm = TRUE)
  } else if (arrange == "var") {
    bb <- apply(matrix, 1, stats::var, na.rm = TRUE)
  }
  if (not_na) {
    na_counts <- apply(matrix, 1, function(x) {
      sum(is.na(x))
    })
    matrix <- matrix[order(na_counts, bb, decreasing = c(FALSE, TRUE)), , drop = FALSE]
  } else {
    matrix <- matrix[order(bb, decreasing = TRUE), , drop = FALSE]
  }

  res <- pheatmap::pheatmap(
    matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_col = factors,
    show_rownames = show_rownames,
    labels_row = .truncate_plot_labels(rownames(matrix), max_rownames_chars),
    border_color = NA,
    silent = TRUE,
    ... = ...
  )

  invisible(res)
}


#' plot heatmap of NA values
#' @param matrix numeric matrix — wide-format intensity data
#' @param annotation data.frame — sample annotation
#' @param factor_keys character vector — factor column names for annotation
#' @param sample_name character — sample name column
#' @param limitrows max rows to display
#' @param distance distance method for clustering
#' @export
#' @keywords internal
#' @family plotting
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' wide <- lfq$data_wide(as.matrix = TRUE)
#' tmp <- plot_na_heatmap(wide$data, wide$annotation, lfq$factor_keys(), lfq$sample_name())
#' stopifnot(class(tmp) == "pheatmap")
#'
plot_na_heatmap <- function(matrix, annotation, factor_keys, sample_name, limitrows = 10000, distance = "binary") {
  stopifnot(annotation[[sample_name]] %in% colnames(matrix))

  factors <- dplyr::select(annotation, all_of(factor_keys))
  factors <- as.data.frame(factors)
  rownames(factors) <- annotation[[sample_name]]

  matrix[!is.na(matrix)] <- 0
  matrix[is.na(matrix)] <- 1
  allrows <- nrow(matrix)
  matrix <- matrix[apply(matrix, 1, sum) > 0, , drop = FALSE]

  message("rows with NA's: ", nrow(matrix), "; all rows :", allrows, "\n")

  if (nrow(matrix) > 1) {
    matrix <- if (nrow(matrix) > limitrows) {
      message("limiting nr of rows to:", limitrows, "\n")
      matrix[sample(seq_len(nrow(matrix)), limitrows), ]
    } else {
      matrix
    }

    gg <- stats::hclust(stats::dist(matrix, method = distance))
    resclust <- pheatmap::pheatmap(
      matrix[gg$order, ],
      cluster_rows = FALSE,
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
#' @param matrix numeric matrix — wide-format intensity data
#' @param annotation data.frame — sample annotation
#' @param sample_name character — sample name column
#' @param factor_keys character vector — factor column names (first for color, second for shape)
#' @param PC which principal components to plot
#' @param add_txt show sample labels
#' @param nudge label nudge distance
#' @export
#' @keywords internal
#' @family plotting
#' @examples
#'
#' istar <- sim_lfq_data_protein_config(with_missing = TRUE, weight_missing = .8, Nprot = 3000)
#' lfq <- LFQData$new(istar$data, istar$config)
#' wide <- lfq$data_wide(as.matrix = TRUE)
#' tmp <- plot_pca(wide$data, wide$annotation, lfq$sample_name(), lfq$factor_keys(),
#'   add_txt = TRUE, nudge = 0.01)
#' stopifnot("ggplot" %in% class(tmp))
#' tmp <- plot_pca(wide$data, wide$annotation, lfq$sample_name(), lfq$factor_keys())
#' stopifnot("ggplot" %in% class(tmp))
#'
plot_pca <- function(matrix, annotation, sample_name, factor_keys, PC = c(1, 2), add_txt = FALSE, nudge = 0.1) {
  stopifnot(length(PC) == 2)

  ff <- matrix
  if (any(is.na(ff))) {
    n_before <- nrow(ff)
    ff <- na.omit(ff)
    n_after <- nrow(ff)
    message(
      "PCA: removed ",
      n_before - n_after,
      " of ",
      n_before,
      " features with missing values. ",
      "For PCA with all features, impute first using impute_with_zcomp()."
    )
  }
  ff <- t(ff)
  pca_result <- prcomp(ff)
  xx <- as_tibble(pca_result$x, rownames = sample_name)
  variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

  if (max(PC) > (ncol(xx) - 1)) {
    warning("nr of PCs: ", ncol(xx), "\n")
    return(NULL)
  }

  xx <- inner_join(annotation, xx)

  sh <- factor_keys[2]
  point <- (if (!is.na(sh)) {
    geom_point(aes(shape = !!sym(sh)))
  } else {
    geom_point()
  })

  text <- geom_text(aes(label = !!sym(sample_name)), check_overlap = TRUE, nudge_x = nudge, nudge_y = nudge)

  pc_x <- paste0("PC", PC[1])
  pc_y <- paste0("PC", PC[2])
  x <- ggplot(
    xx,
    aes(x = !!sym(pc_x), y = !!sym(pc_y), color = !!sym(factor_keys[1]), text = !!sym(sample_name))
  ) +
    labs(
      x = paste0(pc_x, " (", round(variance_explained[PC[1]]), "% variance)"),
      y = paste0(pc_y, " (", round(variance_explained[PC[2]]), "% variance)")
    ) +
    point +
    if (add_txt) {
      text
    }
  if (!is.na(sh)) {
    x <- x + ggplot2::scale_shape_manual(values = seq_along(unique(xx[[sh]])))
  }
  return(x)
}
