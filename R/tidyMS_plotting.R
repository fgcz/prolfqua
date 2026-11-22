#' Prism-style ggplot theme used by all prolfqua plots
#'
#' White panel, black axis lines and outward ticks, bold axis titles, bold facet labels on a light grey strip,
#' legend without box. Add it to your own ggplots to match the look of prolfqua figures.
#' @param grid add faint major gridlines, for plots read against reference values (scatter, density)
#' @param base_size base font size
#' @return a ggplot2 theme
#' @export
#' @family plotting
#' @examples
#' ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
#'   ggplot2::geom_point() +
#'   theme_prolfqua(grid = TRUE)
theme_prolfqua <- function(grid = FALSE, base_size = 12) {
  th <- ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.6),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.6),
      axis.ticks.length = ggplot2::unit(4, "pt"),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "plain")
    )
  if (grid) {
    th <- th + ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.3))
  }
  th
}

#' Decide whether sample counts allow a readable legend
#' @param pdata long-format measurement table
#' @param sample_name column identifying samples
#' @param legend whether to show the plot legend
#' @param max_legend_samples largest sample count shown in a legend
#' @noRd
.show_sample_legend <- function(pdata, sample_name, legend, max_legend_samples) {
  if (is.logical(legend) && length(legend) == 1 && !is.na(legend)) {
    return(legend)
  }
  n_samples <- length(unique(stats::na.omit(pdata[[sample_name]])))
  n_samples <= max_legend_samples
}

#' Label widths for .truncate_plot_labels() and .suffix_plot_labels(); NULL when max_chars is NULL or not finite.
#' @param labels character plot labels
#' @param max_chars maximum displayed label width
#' @param min_chars minimum supported label width
#' @noRd
.plot_label_widths <- function(labels, max_chars, min_chars) {
  if (is.null(max_chars) || !is.finite(max_chars)) {
    return(NULL)
  }
  max_chars <- as.integer(max_chars)
  if (max_chars < min_chars) {
    abort_bad_argument("max_chars", paste0("be at least ", min_chars))
  }
  width <- nchar(labels, type = "chars", allowNA = FALSE)
  list(max_chars = max_chars, width = width, too_long = width > max_chars)
}

#' Shorten long plot labels with an ellipsis
#' @param labels character plot labels
#' @param max_chars maximum displayed label width
#' @noRd
.truncate_plot_labels <- function(labels, max_chars = 60) {
  labels <- as.character(labels)
  lw <- .plot_label_widths(labels, max_chars, 4)
  if (is.null(lw) || !any(lw$too_long)) {
    return(labels)
  }
  too_long <- lw$too_long
  keep_chars <- lw$max_chars - 3
  right_chars <- floor(keep_chars / 2)
  labels[too_long] <- paste0(
    substr(labels[too_long], 1, ceiling(keep_chars / 2)),
    "...",
    substr(labels[too_long], lw$width[too_long] - right_chars + 1, lw$width[too_long])
  )
  labels
}

#' Shorten plot labels while preserving distinguishing suffixes
#' @param labels character plot labels
#' @param max_chars maximum displayed label width
#' @noRd
.suffix_plot_labels <- function(labels, max_chars = 20) {
  labels <- as.character(labels)
  lw <- .plot_label_widths(labels, max_chars, 1)
  if (!is.null(lw)) {
    width <- lw$width[lw$too_long]
    labels[lw$too_long] <- substr(labels[lw$too_long], width - lw$max_chars + 1, width)
  }
  labels
}

# light gray for missing (NA) cells in ComplexHeatmap output
.HEATMAP_NA_COL <- "#D3D3D3"

# point shapes for a second PCA factor: filled first, then hollow, then line shapes
.PCA_SHAPES <- c(16, 17, 15, 18, 1, 2, 0, 5, 6, 3, 4, 8, 7, 9, 10, 12, 13, 11, 14)

#' Colour-blind safe discrete palette for factor levels
#'
#' Okabe-Ito without black for up to 8 levels, viridis beyond; both stay distinguishable under the common forms of
#' colour blindness.
#' @param n number of colours in the returned palette
#' @noRd
.factor_palette <- function(n) {
  if (n > 8L) {
    return(grDevices::hcl.colors(n, "viridis"))
  }
  unname(grDevices::palette.colors(palette = "Okabe-Ito"))[-1L][seq_len(n)]
}

#' Deterministic level -> colour map for each (non-numeric) factor column, shared
#'
#' by .heatmap_top_annotation() and plot_pca() so a given group always gets the
#' same colour across the heatmap, na-heatmap, and PCA legends. Without this,
#' ComplexHeatmap assigns random colours to discrete annotations on every call.
#' @param annotation one row of annotations per sample
#' @param factor_keys configured explanatory-variable columns
#' @noRd
.factor_annotation_colors <- function(annotation, factor_keys) {
  cols <- list()
  for (key in factor_keys) {
    values <- annotation[[key]]
    if (is.numeric(values)) {
      next
    }
    levs <- levels(factor(values))
    cols[[key]] <- stats::setNames(.factor_palette(length(levs)), levs)
  }
  cols
}

#' Build a ComplexHeatmap top annotation from the sample factor columns.
#'
#' ComplexHeatmap aligns the top annotation positionally to the matrix columns
#' (unlike pheatmap's name matching), so the annotation rows must be reordered to
#' match the column order of the matrix handed to Heatmap().
#' @param annotation one row of annotations per sample
#' @param factor_keys configured explanatory-variable columns
#' @param sample_name column identifying samples
#' @param matrix_colnames sample names in heatmap column order
#' @noRd
.heatmap_top_annotation <- function(annotation, factor_keys, sample_name, matrix_colnames) {
  factors <- dplyr::select(annotation, all_of(factor_keys))
  factors <- as.data.frame(factors)
  rownames(factors) <- annotation[[sample_name]]
  factors <- factors[matrix_colnames, , drop = FALSE]
  col <- .factor_annotation_colors(annotation, factor_keys)
  ComplexHeatmap::HeatmapAnnotation(df = factors, col = col, show_annotation_name = TRUE)
}

#' Heatmap() arguments shared by plot_heatmap_cor(), plot_heatmap() and plot_raster(); rows arrive already ordered.
#' @param matrix numeric feature-by-sample matrix
#' @param name heatmap legend title
#' @param col heatmap colour mapping
#' @param top_annotation sample annotation for the heatmap
#' @param max_sample_label_chars maximum sample label width
#' @param ... additional arguments passed to the underlying plotting function
#' @noRd
.annotated_heatmap <- function(matrix, name, col, top_annotation, max_sample_label_chars, ...) {
  ComplexHeatmap::Heatmap(
    matrix,
    name = name,
    col = col,
    na_col = .HEATMAP_NA_COL,
    cluster_rows = FALSE,
    top_annotation = top_annotation,
    show_column_names = TRUE,
    column_labels = .suffix_plot_labels(colnames(matrix), max_sample_label_chars),
    border = FALSE,
    heatmap_legend_param = list(title = name),
    ...
  )
}

#' Blue-white-red diverging colour function centred at zero (colour-blind safe), for row z-scored
#'
#' abundance heatmaps. The range is symmetric around zero based on the data.
#' @param matrix numeric feature-by-sample matrix
#' @noRd
.abundance_col_fun <- function(matrix) {
  m <- max(abs(matrix), na.rm = TRUE)
  if (!is.finite(m) || .is_zero(m)) {
    m <- 1
  }
  circlize::colorRamp2(c(-m, 0, m), c("#2166AC", "#F7F7F7", "#B2182B"))
}

#' Compute pairwise heatmap distances from finite values
#' @param matrix numeric feature-by-sample matrix
#' @noRd
.finite_heatmap_dist <- function(matrix) {
  tryCatch(
    {
      distance <- stats::dist(matrix)
      if (!all(is.finite(distance))) {
        return(NULL)
      }
      distance
    },
    error = function(e) NULL
  )
}

#' Check whether heatmap columns support distance-based clustering
#' @param matrix numeric feature-by-sample matrix
#' @noRd
.can_cluster_heatmap_columns <- function(matrix) {
  if (ncol(matrix) < 2) {
    return(FALSE)
  }
  !is.null(.finite_heatmap_dist(t(matrix)))
}

#' Cluster heatmap rows using distances from finite values
#' @param matrix numeric feature-by-sample matrix
#' @noRd
.cluster_heatmap_rows <- function(matrix) {
  if (nrow(matrix) < 3) {
    return(matrix)
  }
  distance <- .finite_heatmap_dist(matrix)
  if (is.null(distance)) {
    return(matrix)
  }
  matrix[stats::hclust(distance)$order, , drop = FALSE]
}

#' White-to-red color function for correlation heatmaps, ranging over the
#'
#' observed values (R-squared lives in \[0, 1\]).
#' @param cres sample correlation matrix
#' @param R2 show squared correlations
#' @noRd
.cor_col_fun <- function(cres, R2 = FALSE) {
  rng <- range(cres, na.rm = TRUE)
  lo <- if (R2) 0 else max(-1, rng[1])
  hi <- rng[2]
  if (!is.finite(lo) || !is.finite(hi) || lo == hi) {
    lo <- 0
    hi <- 1
  }
  circlize::colorRamp2(c(lo, hi), c("white", "red"))
}

#' Sequential blue-yellow-red color function for the raster of raw abundances.
#' @param matrix numeric feature-by-sample matrix
#' @noRd
.raster_col_fun <- function(matrix) {
  rng <- stats::quantile(matrix, c(0, 1), na.rm = TRUE)
  if (!all(is.finite(rng)) || rng[1] == rng[2]) {
    rng <- c(0, 1)
  }
  circlize::colorRamp2(
    seq(rng[1], rng[2], length.out = 3),
    c("#2c7bb6", "#ffffbf", "#d7191c")
  )
}

#' Visualize intensity distributions
#' @param pdata data.frame
#' @param sample_name character - sample column name
#' @param response character - intensity column name
#' @param is_transformed logical - is intensity log-transformed?
#' @return `plot_intensity_distribution_violin()`: a ggplot object with one violin of the intensities per sample
#'   and the median marked.
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
    theme_prolfqua() +
    geom_violin(linewidth = 0.3) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    stat_summary(fun = median, geom = "point", size = 1, color = "black")
  if (!is_transformed) {
    p <- p + scale_y_log10()
  }
  p
}

#' visualize intensity distributions
#' @param pdata data.frame
#' @param sample_name character - sample column name
#' @param response character - intensity column name
#' @param is_transformed logical - is intensity log-transformed?
#' @param legend show legend. If `NA`, hide automatically when the number of
#'   samples is larger than `max_legend_samples`.
#' @param max_legend_samples maximum number of samples for automatic legend
#'   display.
#' @return `plot_intensity_distribution_density()`: a ggplot object with one density line of the intensities per
#'   sample.
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
  # with plotly, ggplotly() of this plot links the sample curves; without it the plot is static
  density_data <- if (.has_package("plotly")) {
    plotly::highlight_key(
      pdata,
      key = as.character(pdata[[sample_name]]),
      group = paste("prolfqua_density", sample_name, sep = "_")
    )
  } else {
    pdata
  }
  p <- ggplot(
    density_data,
    aes(
      x = .data[[response]],
      colour = .data[[sample_name]],
      group = .data[[sample_name]]
    )
  ) +
    theme_prolfqua(grid = TRUE) +
    geom_line(stat = "density")
  if (!is_transformed) {
    p <- p + scale_x_log10()
  }
  if (!.show_sample_legend(pdata, sample_name, legend, max_legend_samples)) {
    p <- p + ggplot2::guides(colour = "none")
  }
  p
}

#' Visualize correlation among samples
#' @param matrix numeric matrix - wide-format intensity data (samples as columns)
#' @param method correlation method, as for [plot_heatmap_cor()]
#' @param ... passed to [corrplot::corrplot.mixed()], e.g. `tl.cex` and `number.cex` for text size (default 1)
#' @return invisibly, the list returned by [corrplot::corrplot.mixed()] with the correlation matrix `corr` and the
#'   cell positions `corrPos`; called for its side effect of drawing the plot
#'
#' @export
#' @keywords internal
#' @family plotting
#' @rdname plot_sample_correlation
#' @examplesIf requireNamespace("corrplot", quietly = TRUE)
#' istar <- sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' lfq <- lfq$get_Transformer()$log2()$lfq
#' plot_sample_correlation(lfq$data_wide(as.matrix = TRUE)$data)
plot_sample_correlation <- function(matrix, method = "pearson", ...) {
  rlang::check_installed("corrplot", reason = "for sample correlation plots.")
  M <- .sample_correlation(matrix, method)
  # up to 12 samples the lower panel shows the correlation values, otherwise pies (corrplot defaults otherwise)
  numbers <- nrow(M) <= 12
  res <- corrplot::corrplot.mixed(
    M,
    upper = "ellipse",
    lower = if (numbers) "number" else "pie",
    lower.col = if (numbers) "black",
    diag = "u",
    tl.pos = "lt",
    tl.col = "black",
    mar = c(2, 5, 5, 2),
    ...
  )
  invisible(res)
}

#' Boxplot of feature intensities by factor levels
#'
#' @param pdata data.frame
#' @param title name to show
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param facet_grid_on on which variable to run facet_grid
#' @param beeswarm use beeswarm default FALSE
#' @param show_mean show mean values
#' @param pb progress bar
#' @return a ggplot object with boxplots of the intensities in `pdata` per combination of factor levels, optionally
#'   with the individual points and the means
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
  p <- ggplot(pdata, aes(x = .data[["interaction"]], y = .data[[response]])) +
    theme_prolfqua() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    ggtitle(title)
  if (lil > 1) {
    p <- p + aes(colour = .data[[isotope_col]])
  }

  if (!lfqdata$is_transformed()) {
    p <- p + scale_y_log10()
  }
  p <- p + geom_boxplot()
  if (beeswarm) {
    if (nlevels(pdata[[nr_children_col]]) > 1) {
      shape_values <- c(4, rep(16, nlevels(pdata[[nr_children_col]]) - 1))
    } else {
      shape_values <- 16
    }
    point_aes <- aes(size = .data[["size"]], shape = .data[[nr_children_col]])
    points <- if (.has_package("ggbeeswarm")) {
      ggbeeswarm::geom_quasirandom(point_aes, dodge.width = 0.7)
    } else {
      geom_point(point_aes, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7))
    }
    p <- p +
      points +
      scale_shape_manual(values = shape_values) +
      scale_size_continuous(range = c(1, 4))
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
    p <- p + facet_grid(cols = vars(.data[[facet_grid_on]]))
  }
  p
}


#' Generate feature-level boxplots for all proteins
#' @export
#' @param pdata data.frame
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param hierarchy e.g. protein_Id default relevant_hierarchy_keys
#' @param facet_grid_on default NULL
#' @return a tibble with one row per element of `hierarchy`, a column (named after the `hierarchy` columns joined
#'   by `+`) identifying it, and a list column `boxplot` of ggplot objects made by [plot_hierarchies_boxplot()]
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
  xnested <- xnested |> tidyr::unite(!!newcol, all_of(hierarchy))

  pb <- .make_progress(nrow(xnested))

  figs <- xnested |>
    dplyr::mutate(
      boxplot = map2(
        data,
        .data[[newcol]],
        plot_hierarchies_boxplot,
        lfqdata = lfqdata,
        facet_grid_on = facet_grid_on,
        pb = pb
      )
    )
  dplyr::select(figs, all_of(c(newcol, "boxplot")))
}


#' Sample-by-sample correlation matrix
#' @param matrix numeric feature-by-sample matrix
#' @param method a stats::cor() method name or a function of the matrix returning a correlation matrix
#' @noRd
#' @srrstats {G3.1, G3.1a} Sample correlations accept any `stats::cor()` method or a user function (examples of
#'   `plot_heatmap_cor()`; test-tidyMS_plotting.R). Coefficient covariances come from the estimator of the chosen
#'   `build_contrast_analysis()` method (least squares, M-estimation, rank-based, mixed, moderated, Firth) and are
#'   returned by `vcov()`.
.sample_correlation <- function(matrix, method) {
  if (is.function(method)) {
    result <- method(matrix)
    if (!is.matrix(result) || !is.numeric(result) || !identical(dim(result), rep(ncol(matrix), 2L))) {
      abort_bad_argument("method", "return a numeric sample-by-sample matrix when it is a function")
    }
    return(result)
  }
  method <- .arg_match_prolfqua(method, c("pearson", "spearman", "kendall"), "method")
  stats::cor(matrix, use = "pairwise.complete.obs", method = method)
}

#' Plot a sample correlation heatmap with annotations
#'
#' @param matrix numeric matrix - wide-format intensity data
#' @param annotation data.frame - sample annotation
#' @param factor_keys character vector - factor column names for annotation
#' @param sample_name character - sample name column
#' @param R2 logical - plot R-squared instead of correlation
#' @param color color palette
#' @param max_sample_label_chars maximum displayed sample label length. Labels
#'   keep their suffix because sample prefixes are often shared.
#' @param method correlation method: `"pearson"` (default), `"spearman"` or `"kendall"` for [stats::cor()] on
#'   pairwise complete observations, or a function of the matrix returning a sample-by-sample correlation matrix
#' @param ... passed to [ComplexHeatmap::Heatmap()]
#' @return invisibly, a ComplexHeatmap `Heatmap` object of the clustered sample correlations (or R^2 values)
#' @export
#' @keywords internal
#' @family plotting
#' @examplesIf requireNamespace("ComplexHeatmap", quietly = TRUE) && requireNamespace("circlize", quietly = TRUE)
#'
#' istar <- sim_lfq_data_protein_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' wide <- lfq$data_wide(as.matrix = TRUE)
#' pheat_map <- plot_heatmap_cor(wide$data, wide$annotation,
#'   lfq$factor_keys(), lfq$sample_name())
#' stopifnot(methods::is(pheat_map, "Heatmap"))
#' pheat_map <- plot_heatmap_cor(wide$data, wide$annotation,
#'   lfq$factor_keys(), lfq$sample_name(), R2 = TRUE)
#' stopifnot(methods::is(pheat_map, "Heatmap"))
#' # rank correlation, or any function returning a correlation matrix
#' spearman <- plot_heatmap_cor(wide$data, wide$annotation,
#'   lfq$factor_keys(), lfq$sample_name(), method = "spearman")
#' complete <- plot_heatmap_cor(wide$data, wide$annotation,
#'   lfq$factor_keys(), lfq$sample_name(), method = function(x) stats::cor(x, use = "complete.obs"))
#'
plot_heatmap_cor <- function(
  matrix,
  annotation,
  factor_keys,
  sample_name,
  R2 = FALSE,
  color = colorRampPalette(c("white", "red"))(1024),
  max_sample_label_chars = 20,
  method = "pearson",
  ...
) {
  rlang::check_installed(c("ComplexHeatmap", "circlize"), reason = "for heatmaps.")
  cres <- .sample_correlation(matrix, method)
  if (R2) {
    cres <- cres^2
  }

  gg <- stats::hclust(stats::dist(cres))
  plot_data <- cres[gg$order, ]
  title <- if (R2) "R^2" else "correlation"
  res <- .annotated_heatmap(
    plot_data,
    title,
    .cor_col_fun(plot_data, R2),
    .heatmap_top_annotation(annotation, factor_keys, sample_name, colnames(plot_data)),
    max_sample_label_chars,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    column_title = title,
    ...
  )
  invisible(res)
}


#' Plot a row z-score heatmap with annotations
#'
#' @param matrix numeric matrix - wide-format intensity data
#' @param annotation data.frame - sample annotation
#' @param factor_keys character vector - factor column names for annotation
#' @param sample_name character - sample name column
#' @param na_fraction maximum fraction of missing values per row. Rows that
#'   cannot be meaningfully z-scored are always removed.
#' @param show_rownames if TRUE shows row names, default FALSE
#' @param max_rownames_chars maximum displayed row label length
#' @param max_sample_label_chars maximum displayed sample label length. Labels
#'   keep their suffix because sample prefixes are often shared.
#' @param ... passed to [ComplexHeatmap::Heatmap()]
#' @return invisibly, a ComplexHeatmap `Heatmap` object of the row z-scored data; `NULL` with a warning if
#'   `matrix` has no rows or no row can be z-scored
#' @export
#' @keywords internal
#' @family plotting
#' @examplesIf requireNamespace("ComplexHeatmap", quietly = TRUE) && requireNamespace("circlize", quietly = TRUE)
#'
#' istar <- sim_lfq_data_protein_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' wide <- lfq$data_wide(as.matrix = TRUE)
#' p <- plot_heatmap(wide$data, wide$annotation, lfq$factor_keys(), lfq$sample_name())
#' stopifnot(methods::is(p, "Heatmap"))
#'
plot_heatmap <- function(
  matrix,
  annotation,
  factor_keys,
  sample_name,
  na_fraction = 0.4,
  show_rownames = FALSE,
  max_rownames_chars = 60,
  max_sample_label_chars = 20,
  ...
) {
  rlang::check_installed(c("ComplexHeatmap", "circlize"), reason = "for heatmaps.")
  if (nrow(matrix) == 0) {
    warning("plot_heatmap: the dataset has no rows.")
    return(NULL)
  }

  resdata <- t(scale(t(matrix)))
  na_threshold <- floor(ncol(resdata) * na_fraction)
  can_z_score <- rowSums(is.finite(resdata)) >= 2
  keep_rows <- rowSums(is.na(resdata)) <= na_threshold & can_z_score
  resdataf <- resdata[keep_rows, , drop = FALSE]
  if (nrow(resdataf) == 0) {
    warning("No features have enough observed variation for a row z-scored heatmap.")
    return(NULL)
  }
  plot_data <- .cluster_heatmap_rows(resdataf)
  res <- .annotated_heatmap(
    plot_data,
    "row z-score",
    .abundance_col_fun(plot_data),
    .heatmap_top_annotation(annotation, factor_keys, sample_name, colnames(plot_data)),
    max_sample_label_chars,
    cluster_columns = .can_cluster_heatmap_columns(plot_data),
    show_row_names = show_rownames,
    row_labels = .truncate_plot_labels(rownames(plot_data), max_rownames_chars),
    use_raster = FALSE,
    ...
  )
  invisible(res)
}


#' Plot a heatmap without clustering (use to show NAs)
#' @param matrix numeric matrix - wide-format intensity data
#' @param annotation data.frame - sample annotation
#' @param factor_keys character vector - factor column names for annotation
#' @param sample_name character - sample name column
#' @param arrange either mean or var
#' @param not_na if true than arrange by nr of NA's first and then by arrange
#' @param show_rownames logical, show row names in heatmap
#' @param max_rownames_chars maximum displayed row label length
#' @param max_sample_label_chars maximum displayed sample label length. Labels
#'   keep their suffix because sample prefixes are often shared.
#' @param ... additional arguments passed to [ComplexHeatmap::Heatmap()]
#' @return invisibly, a ComplexHeatmap `Heatmap` object of the abundances with rows ordered by `arrange` and
#'   unclustered columns; `NULL` with a warning if `matrix` has fewer than two rows
#' @keywords internal
#'
#' @family plotting
#' @export
#' @examplesIf requireNamespace("ComplexHeatmap", quietly = TRUE) && requireNamespace("circlize", quietly = TRUE)
#'
#' istar <- sim_lfq_data_protein_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' wide <- lfq$data_wide(as.matrix = TRUE)
#' rs <- plot_raster(wide$data, wide$annotation, lfq$factor_keys(), lfq$sample_name())
#' stopifnot(methods::is(rs, "Heatmap"))
#' rs <- plot_raster(wide$data, wide$annotation, lfq$factor_keys(), lfq$sample_name(), "var")
#' stopifnot(methods::is(rs, "Heatmap"))
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
  max_sample_label_chars = 20,
  ...
) {
  rlang::check_installed(c("ComplexHeatmap", "circlize"), reason = "for heatmaps.")
  if (nrow(matrix) <= 1) {
    warning("plot_raster: the dataset needs at least two rows; got ", nrow(matrix), ".")
    return(NULL)
  }
  arrange <- match.arg(arrange)

  bb <- apply(matrix, 1, if (arrange == "mean") mean else stats::var, na.rm = TRUE)
  if (not_na) {
    matrix <- matrix[order(rowSums(is.na(matrix)), bb, decreasing = c(FALSE, TRUE)), , drop = FALSE]
  } else {
    matrix <- matrix[order(bb, decreasing = TRUE), , drop = FALSE]
  }

  res <- .annotated_heatmap(
    matrix,
    "abundance",
    .raster_col_fun(matrix),
    .heatmap_top_annotation(annotation, factor_keys, sample_name, colnames(matrix)),
    max_sample_label_chars,
    cluster_columns = FALSE,
    show_row_names = show_rownames,
    row_labels = .truncate_plot_labels(rownames(matrix), max_rownames_chars),
    use_raster = FALSE,
    ...
  )
  invisible(res)
}


#' Plot a heatmap of missing values
#' @param matrix numeric matrix - wide-format intensity data
#' @param annotation data.frame - sample annotation
#' @param factor_keys character vector - factor column names for annotation
#' @param sample_name character - sample name column
#' @param limitrows max rows to display
#' @param distance distance method for clustering
#' @return a ComplexHeatmap `Heatmap` object marking the missing values of the rows that have any; `NULL` if fewer
#'   than two such rows
#' @export
#' @keywords internal
#' @family plotting
#' @examplesIf requireNamespace("ComplexHeatmap", quietly = TRUE) && requireNamespace("circlize", quietly = TRUE)
#'
#' istar <- sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' wide <- lfq$data_wide(as.matrix = TRUE)
#' tmp <- plot_na_heatmap(wide$data, wide$annotation, lfq$factor_keys(), lfq$sample_name())
#' stopifnot(methods::is(tmp, "Heatmap"))
#'
plot_na_heatmap <- function(matrix, annotation, factor_keys, sample_name, limitrows = 10000, distance = "binary") {
  rlang::check_installed(c("ComplexHeatmap", "circlize"), reason = "for heatmaps.")
  if (!all(annotation[[sample_name]] %in% colnames(matrix))) {
    abort_bad_argument("annotation", "name only samples that are columns of `matrix`")
  }

  matrix[!is.na(matrix)] <- 0
  matrix[is.na(matrix)] <- 1
  matrix <- matrix[rowSums(matrix) > 0, , drop = FALSE]

  if (nrow(matrix) > 1) {
    matrix <- if (nrow(matrix) > limitrows) {
      inform_prolfqua("limiting nr of rows to: ", limitrows)
      matrix[sample(seq_len(nrow(matrix)), limitrows), ]
    } else {
      matrix
    }

    gg <- stats::hclust(stats::dist(matrix, method = distance))
    plot_data <- matrix[gg$order, ]
    resclust <- ComplexHeatmap::Heatmap(
      plot_data,
      name = "missing",
      col = circlize::colorRamp2(c(0, 1), c("white", "black")),
      cluster_rows = FALSE,
      cluster_columns = TRUE,
      clustering_distance_columns = distance,
      top_annotation = .heatmap_top_annotation(annotation, factor_keys, sample_name, colnames(plot_data)),
      show_row_names = FALSE,
      show_column_names = TRUE,
      border = FALSE,
      use_raster = FALSE,
      show_heatmap_legend = FALSE
    )
    resclust
  } else {
    NULL
  }
}


#' Plot a principal component analysis of the samples
#' @param matrix numeric matrix - wide-format intensity data
#' @param annotation data.frame - sample annotation
#' @param sample_name character - sample name column
#' @param factor_keys character vector - factor column names (first for color, second for shape)
#' @param PC which principal components to plot
#' @param add_txt show sample labels
#' @param nudge label nudge distance
#' @return a ggplot object of the sample scores on the principal components `PC`, coloured by the first and shaped
#'   by the second of `factor_keys`
#' @export
#' @keywords internal
#' @family plotting
#' @examples
#'
#' istar <- sim_lfq_data_protein_config(with_missing = TRUE, weight_missing = .8, n_proteins = 3000)
#' lfq <- LFQData$new(istar$data, istar$config)
#' wide <- lfq$data_wide(as.matrix = TRUE)
#' tmp <- plot_pca(wide$data, wide$annotation, lfq$sample_name(), lfq$factor_keys(),
#'   add_txt = TRUE, nudge = 0.01)
#' stopifnot("ggplot" %in% class(tmp))
#' tmp <- plot_pca(wide$data, wide$annotation, lfq$sample_name(), lfq$factor_keys())
#' stopifnot("ggplot" %in% class(tmp))
#'
plot_pca <- function(matrix, annotation, sample_name, factor_keys, PC = c(1, 2), add_txt = FALSE, nudge = 0.1) {
  if (length(PC) != 2) {
    abort_bad_argument("PC", "be two principal component numbers", not = toString(PC))
  }
  impute_hint <- "missing values first, e.g. AggregateLimpa$new(lfqdata, impute_only = TRUE)$aggregate()."

  # Duplicated sample names become duplicated row names after transpose and
  # cause a cartesian expansion in the join below - error early instead.
  dup_samples <- unique(colnames(matrix)[duplicated(colnames(matrix))])
  if (length(dup_samples) > 0) {
    abort_prolfqua(sprintf("PCA: duplicated sample names are not allowed: %s.", toString(dup_samples)))
  }

  if (anyNA(matrix)) {
    n_before <- nrow(matrix)
    matrix <- na.omit(matrix)
    inform_prolfqua(sprintf(
      "PCA: removed %s of %s features with missing values. To keep all features, impute %s",
      n_before - nrow(matrix),
      n_before,
      impute_hint
    ))
  }

  # prcomp() errors on a 0-row matrix, which happens when every feature had at
  # least one missing value. Fail with an actionable message instead.
  if (nrow(matrix) == 0) {
    abort_prolfqua("PCA: no features without missing values remain after NA filtering. Impute ", impute_hint)
  }

  # Centered data has rank at most (n_samples - 1), so component max(PC) only
  # exists when there are at least max(PC) + 1 samples. Guard before prcomp() so
  # the caller (e.g. pca_plotly()) gets a clear error rather than a NULL plot.
  if (ncol(matrix) <= max(PC)) {
    abort_prolfqua(sprintf(
      "PCA: need at least %s samples to plot principal components %s and %s; got %s.",
      max(PC) + 1,
      PC[1],
      PC[2],
      ncol(matrix)
    ))
  }

  pca_result <- prcomp(t(matrix), center = TRUE, scale. = FALSE)
  xx <- as_tibble(pca_result$x, rownames = sample_name)
  variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

  if (max(PC) > (ncol(xx) - 1)) {
    abort_prolfqua(sprintf(
      "PCA: requested principal component %s but only %s principal component(s) could be computed %s",
      max(PC),
      ncol(xx) - 1,
      "(too few features or samples)."
    ))
  }

  xx <- inner_join(annotation, xx, by = sample_name, relationship = "many-to-one")

  sh <- factor_keys[2]
  if (!is.na(sh) && length(unique(stats::na.omit(xx[[sh]]))) > length(.PCA_SHAPES)) {
    inform_prolfqua(
      "plot_pca: `",
      sh,
      "` has more than ",
      length(.PCA_SHAPES),
      " levels; points are not shaped by it."
    )
    sh <- NA
  }
  point <- (if (!is.na(sh)) {
    geom_point(aes(shape = .data[[sh]]))
  } else {
    geom_point()
  })

  text <- ggrepel::geom_text_repel(
    aes(label = .data[[sample_name]]),
    min.segment.length = 0,
    box.padding = 0.35,
    point.padding = 0.25,
    seed = 1,
    nudge_x = nudge,
    nudge_y = nudge,
    show.legend = FALSE
  )

  pc_x <- paste0("PC", PC[1])
  pc_y <- paste0("PC", PC[2])
  # Share the heatmap annotation palette so a group keeps one colour everywhere.
  color_palette <- .factor_annotation_colors(annotation, factor_keys[1])[[factor_keys[1]]]
  x <- ggplot(
    xx,
    aes(x = .data[[pc_x]], y = .data[[pc_y]], color = .data[[factor_keys[1]]], text = .data[[sample_name]])
  ) +
    theme_prolfqua(grid = TRUE) +
    labs(
      x = paste0(pc_x, " (", round(variance_explained[PC[1]]), "% variance)"),
      y = paste0(pc_y, " (", round(variance_explained[PC[2]]), "% variance)")
    ) +
    point +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.12)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.12)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 24, 8, 24)) +
    if (add_txt) {
      text
    }
  if (!is.null(color_palette)) {
    x <- x + ggplot2::scale_colour_manual(values = color_palette)
  }
  if (!is.na(sh)) {
    x <- x + ggplot2::scale_shape_manual(values = .PCA_SHAPES)
  }
  x
}
