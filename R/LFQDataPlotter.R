#' Keep the `n_most_variable` most variable features for the intensity heatmap.
#'
#' Ranks features by their prolfqua variability statistic (CV for untransformed
#' data, sd for transformed data; see LFQDataStats) and keeps the most variable.
#' The heatmap clusters rows with stats::hclust, which errors above 65536
#' features, so this keeps clustering feasible (and the plot legible) for
#' peptide-list / entrapment searches with tens of thousands of degenerate
#' protein groups. `n_most_variable = NULL` / non-finite keeps all features.
#' Returns a subset of wide$data in the original row order.
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param wide wide-format data, row identifiers and sample annotations
#' @param n_most_variable maximum number of variable features to retain
#' @noRd
.select_most_variable_features <- function(lfqdata, wide, n_most_variable) {
  mat <- wide$data
  if (is.null(n_most_variable) || !is.finite(n_most_variable) || nrow(mat) <= n_most_variable) {
    return(mat)
  }
  lfqstats <- lfqdata$get_Stats(stats = "all")
  metric <- lfqstats$stat # "CV" (untransformed) or "sd" (transformed)
  keys <- lfqdata$hierarchy_keys()
  # wide$rowdata is in the same row order as the matrix, so joining the
  # per-feature statistic onto it aligns the ranking with the matrix rows.
  ranked <- dplyr::left_join(
    wide$rowdata,
    dplyr::select(lfqstats$stats(), dplyr::all_of(c(keys, metric))),
    by = keys,
    relationship = "one-to-one"
  )
  keep <- utils::head(order(ranked[[metric]], decreasing = TRUE, na.last = TRUE), n_most_variable)
  mat[sort(keep), , drop = FALSE]
}

#' Hierarchy level below the modelled one to facet boxplots on, or NULL
#' @param lfq LFQData
#' @param facet facet when there is a lower level
#' @noRd
.boxplot_facet <- function(lfq, facet) {
  rhk <- lfq$relevant_hierarchy_keys()
  hk <- lfq$hierarchy_keys()
  if (length(rhk) < length(hk) && facet) hk[length(rhk) + 1] else NULL
}

#' LFQDataPlotter ----
#' Create various visualization of the LFQdata
#' @return An R6 class generator.
#' @export
#'
#' @family LFQData
#' @import dplyr
#' @examplesIf requireNamespace("ComplexHeatmap", quietly = TRUE) && requireNamespace("circlize", quietly = TRUE)
#'
#' istar <- sim_lfq_data_peptide_config()
#'
#' lfqdata <- LFQData$new(
#'  istar$data,
#'  istar$config)
#' lfqplotter <- lfqdata$get_Plotter()
#'
#' stopifnot(methods::is(lfqplotter$heatmap(), "Heatmap"))
#' stopifnot(methods::is(lfqplotter$heatmap_cor(), "Heatmap"))
#' stopifnot("ggplot" %in% class(lfqplotter$pca()))
#' stopifnot("plotly" %in%  class(lfqplotter$pca_plotly()))
#' tmp <- lfqplotter$boxplots()
#' stopifnot("ggplot" %in%  class(tmp$boxplot[[1]]))
#' stopifnot("ggplot" %in% class(lfqplotter$missingness_histogram()))
#'
#' stopifnot(methods::is(lfqplotter$na_heatmap(), "Heatmap"))
#' class(lfqplotter$intensity_distribution_density())
#' class(lfqplotter$intensity_distribution_violin())
#' stopifnot(is.null(lfqplotter$pairs_smooth()))
#' stopifnot(class(lfqplotter$sample_correlation()) == "list")
#' stopifnot(methods::is(lfqplotter$raster(), "Heatmap"))
#' stopifnot("upset" == class(lfqplotter$upset_missing()))
#'
#'
LFQDataPlotter <- R6::R6Class(
  "LFQDataPlotter",
  list(
    #' @field lfq LFQData object
    #' @field prefix prefix to figure names when writing, e.g. protein_
    lfq = NULL,
    prefix = "",
    #' @description
    #' create LFQDataPlotter
    #' @param lfqdata [LFQData] object with the data and its configuration
    #' @param prefix will be prepended to outputs written
    initialize = function(lfqdata, prefix = "ms_") {
      self$lfq <- lfqdata$clone(deep = TRUE)
      self$prefix <- prefix
    },

    #' @description
    #' plot intensities in raster
    #' @param arrange arrange by either mean or var
    #' @param not_na TRUE arrange by number of NA's, FALSE by arrange by intensity
    #' @param rownames show rownames (default FALSE - do not show.)
    #' @param max_rownames_chars maximum displayed row label length
    #' @param max_sample_label_chars maximum displayed sample label length.
    #'   Labels keep their suffix because sample prefixes are often shared.
    #' @return ComplexHeatmap::Heatmap
    raster = function(
      arrange = c("mean", "var"),
      not_na = FALSE,
      rownames = FALSE,
      max_rownames_chars = 60,
      max_sample_label_chars = 20
    ) {
      arrange <- match.arg(arrange)
      wide <- self$lfq$data_wide(as.matrix = TRUE)
      fig <- .plot_raster(
        wide$data,
        wide$annotation,
        self$lfq$factor_keys(),
        self$lfq$sample_name(),
        arrange = arrange,
        not_na = not_na,
        show_rownames = rownames,
        max_rownames_chars = max_rownames_chars,
        max_sample_label_chars = max_sample_label_chars
      )
      fig
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
    #' By default only the \code{n_most_variable} most variable features are shown and the rows are clustered.
    #' Row clustering uses \code{stats::hclust}, which errors above 65536 features, so the rows are ranked by
    #' their variability statistic (CV for untransformed data, sd for transformed data; see
    #' \code{\link{LFQDataStats}}) and the most variable are kept. This keeps the heatmap feasible and legible for
    #' peptide-list / entrapment searches with tens of thousands of features.
    #'
    #' With \code{row_order}, the heatmap shows exactly these features in this order, without clustering. Use it
    #' to compare two states of the same data row by row, e.g. before and after imputation: imputation changes the
    #' per-feature sd, so the most variable set and the clustering would differ between the two heatmaps.
    #'
    #' @param na_fraction maximum fraction of missing values per row. The default
    #'   keeps every row that can be meaningfully z-scored.
    #' @param rownames show rownames (default FALSE - do not show.)
    #' @param max_rownames_chars maximum displayed row label length
    #' @param max_sample_label_chars maximum displayed sample label length.
    #'   Labels keep their suffix because sample prefixes are often shared.
    #' @param n_most_variable keep the \code{n_most_variable} most variable features (default 1000);
    #'   \code{NULL} or \code{Inf} keeps all features. Not together with \code{row_order}.
    #' @param row_order feature IDs (rownames of \code{LFQData$data_wide()}) to show, in this order, without row
    #'   clustering; default \code{NULL} selects and clusters the most variable features.
    #' @param top_n deprecated, use \code{n_most_variable}
    #' @return ComplexHeatmap::Heatmap
    heatmap = function(
      na_fraction = 1,
      rownames = FALSE,
      max_rownames_chars = 60,
      max_sample_label_chars = 20,
      n_most_variable = 1000,
      row_order = NULL,
      top_n
    ) {
      if (!missing(top_n)) {
        n_most_variable <- .deprecate_arg("top_n", "n_most_variable", top_n)
      }
      wide <- self$lfq$data_wide(as.matrix = TRUE)
      if (is.null(row_order)) {
        data <- .select_most_variable_features(self$lfq, wide, n_most_variable)
      } else {
        if (!missing(n_most_variable) || !missing(top_n)) {
          abort_bad_argument("row_order", "be NULL when `n_most_variable` is given; use one of the two")
        }
        unknown <- setdiff(row_order, rownames(wide$data))
        if (length(unknown) > 0L) {
          abort_bad_argument(
            "row_order",
            "hold feature IDs of `data_wide()`",
            not = paste0("unknown IDs ", toString(unknown))
          )
        }
        data <- wide$data[row_order, , drop = FALSE]
      }
      fig <- .plot_heatmap(
        data,
        wide$annotation,
        self$lfq$factor_keys(),
        self$lfq$sample_name(),
        na_fraction = na_fraction,
        show_rownames = rownames,
        max_rownames_chars = max_rownames_chars,
        max_sample_label_chars = max_sample_label_chars,
        cluster_rows = is.null(row_order)
      )
      fig
    },
    #' @description
    #' heatmap of sample correlations.
    #'
    #' The Spearman correlation among all samples
    #' is computed. Then the euclidean distance is used to compute the distances.
    #'
    #'
    #' @param max_sample_label_chars maximum displayed sample label length.
    #'   Labels keep their suffix because sample prefixes are often shared.
    #' @param method correlation method: \code{"pearson"} (default), \code{"spearman"} or \code{"kendall"} for
    #'   \code{stats::cor()} on pairwise complete observations, or a function of the matrix returning a
    #'   sample-by-sample correlation matrix
    #' @return ComplexHeatmap::Heatmap
    #'
    heatmap_cor = function(max_sample_label_chars = 20, method = "pearson") {
      wide <- self$lfq$data_wide(as.matrix = TRUE)
      fig <- .plot_heatmap_cor(
        wide$data,
        wide$annotation,
        self$lfq$factor_keys(),
        self$lfq$sample_name(),
        max_sample_label_chars = max_sample_label_chars,
        method = method
      )
      fig
    },
    #' @description
    #' PCA plot
    #'
    #' A PCA is applied and the first and second principal component are shown.
    #' Features with missing values are removed. To keep all features, impute
    #' first, e.g. with \code{\link{AggregateLimpa}} in \code{impute_only} mode.
    #'
    #'
    #' @param add_txt show sample names
    #' @param PC default c(1,2) - first and second principal component
    #' @param nudge default 0.1 nudge point labels
    #' @return ggplot
    pca = function(PC = c(1, 2), add_txt = TRUE, nudge = 0.1) {
      wide <- self$lfq$data_wide(as.matrix = TRUE)
      fig <- .plot_pca(
        wide$data,
        wide$annotation,
        self$lfq$sample_name(),
        self$lfq$factor_keys(),
        PC = PC,
        add_txt = add_txt,
        nudge = nudge
      )
      fig
    },
    #' @description
    #' pca plot
    #' @param add_txt show sample names
    #' @param PC default c(1,2) - first and second principal component
    #' @return plotly
    pca_plotly = function(PC = c(1, 2), add_txt = FALSE) {
      rlang::check_installed("plotly", reason = "for interactive plots.")
      fig <- plotly::ggplotly(self$pca(PC = PC, add_txt = add_txt), tooltip = self$lfq$sample_name())
      fig
    },
    #' @description
    #' boxplots for all proteins
    #' @param facet enable facet wrap if hierarchy_depth less then hierarchy length.
    #' @return tibble with column boxplots containing ggplot objects
    boxplots = function(facet = TRUE) {
      .plot_hierarchies_boxplot_df(
        self$lfq$data_long(),
        self$lfq,
        hierarchy = self$lfq$relevant_hierarchy_keys(),
        facet_grid_on = .boxplot_facet(self$lfq, facet)
      )
    },
    #' @description
    #' histogram of intensities given number of missing in conditions
    #' @return ggplot
    missingness_histogram = function() {
      .missingness_histogram(self$lfq)
    },
    #' @description
    #' Deprecated; use `$missingness_histogram()`.
    #' @return ggplot
    missigness_histogram = function() {
      .Deprecated(msg = "`missigness_histogram()` is deprecated; use `missingness_histogram()` instead.")
      self$missingness_histogram()
    },

    #' @description
    #' heatmap of features with missing values
    #' @return ComplexHeatmap::Heatmap
    na_heatmap = function() {
      wide <- self$lfq$data_wide(as.matrix = TRUE)
      .plot_na_heatmap(
        wide$data,
        wide$annotation,
        self$lfq$factor_keys(),
        self$lfq$sample_name()
      )
    },
    #' @description
    #' density distribution of intensities
    #' @param legend show legend TRUE, FALSE do not show, NA selects
    #'   automatically based on sample count.
    #' @param max_legend_samples maximum number of samples for automatic legend
    #'   display.
    #' @return ggplot
    intensity_distribution_density = function(legend = NA, max_legend_samples = 16) {
      .plot_intensity_distribution_density(
        self$lfq$data_long(),
        self$lfq$sample_name(),
        self$lfq$response(),
        self$lfq$is_transformed(),
        legend = legend,
        max_legend_samples = max_legend_samples
      )
    },
    #' @description
    #' Violinplot showing distribution of intensities in all samples
    #' @return ggplot
    intensity_distribution_violin = function() {
      .plot_intensity_distribution_violin(
        self$lfq$data_long(),
        self$lfq$sample_name(),
        self$lfq$response(),
        self$lfq$is_transformed()
      )
    },
    #' @description
    #' pairsplot of intensities
    #' @param max maximal number of samples to show
    #' @return NULL
    pairs_smooth = function(max = 10) {
      lfq <- self$lfq
      sample_col <- lfq$sample_name()
      samples <- dplyr::select(lfq$data_long(), dplyr::all_of(sample_col)) |>
        distinct() |>
        pull()
      if (length(samples) > max) {
        limit <- samples |> sample(max)
        lfq <- lfq$get_copy()
        lfq$set_data(dplyr::filter(lfq$data_long(), .data[[sample_col]] %in% .env$limit))
      }
      .pairs_smooth(lfq$data_wide(as.matrix = TRUE)$data)
      NULL
    },
    #' @description
    #' plot of sample correlations
    #' @param method correlation method: \code{"pearson"} (default), \code{"spearman"} or \code{"kendall"} for
    #'   \code{stats::cor()} on pairwise complete observations, or a function of the matrix returning a
    #'   sample-by-sample correlation matrix
    #' @return NULL
    sample_correlation = function(method = "pearson") {
      .plot_sample_correlation(self$lfq$data_wide(as.matrix = TRUE)$data, method = method)
    },
    #' @description
    #' upset plot based on presence absence information
    #' @return plot
    upset_missing = function() {
      rlang::check_installed("UpSetR", reason = "for UpSet plots.")
      pups <- .upset_missing_stats(self$lfq)
      res <- UpSetR::upset(pups$data, order.by = "freq", nsets = pups$nsets)
      res
    },
    #' @description
    #' write boxplots to file
    #' @param path_qc path to write to
    #' @param filename file name; one ending in `.pdf` is used as given, otherwise `boxplot.pdf` is appended
    #'   (default: the plotter prefix)
    #' @param width fig width
    #' @param height fig height
    #' @return path the file was written to (invisible)
    write_boxplots = function(path_qc, filename = NULL, width = 6, height = 6) {
      if (is.null(filename)) {
        filename <- self$prefix
      }
      fpath <- file.path(path_qc, if (.has_pdf_suffix(filename)) filename else paste0(filename, "boxplot.pdf"))
      # one plot at a time: memory does not grow with the number of proteins
      hierarchy <- self$lfq$relevant_hierarchy_keys()
      xnested <- .nest_by_hierarchy(self$lfq$data_long(), hierarchy)
      ids <- xnested[[paste(hierarchy, collapse = "+")]]
      facet_grid_on <- .boxplot_facet(self$lfq, TRUE)

      inform_prolfqua("write_boxplots: writing ", fpath)
      pb <- .make_progress(nrow(xnested))

      grDevices::pdf(fpath, width = width, height = height)
      pdf_dev <- grDevices::dev.cur()
      on.exit(.dev_close(pdf_dev), add = TRUE)
      for (i in seq_len(nrow(xnested))) {
        pb$tick()
        .render_plot_to_device(
          .plot_hierarchies_boxplot(xnested$data[[i]], ids[i], self$lfq, facet_grid_on = facet_grid_on)
        )
      }
      invisible(fpath)
    },
    #' @description
    #' write figure to pdf
    #' @param fig ggplot or ComplexHeatmap::Heatmap
    #' @param path_qc path to write to
    #' @param fig_name name of figure; `.pdf` is appended unless it ends in `.pdf`
    #' @param width figure width
    #' @param height figure height
    #' @return path the file was written to
    write_pdf = function(fig, path_qc, fig_name, width = 7, height = 7) {
      fpath <- file.path(path_qc, paste0(self$prefix, fig_name, if (!.has_pdf_suffix(fig_name)) ".pdf"))
      inform_prolfqua("write_pdf: writing ", fpath)
      grDevices::pdf(fpath, width = width, height = height)
      pdf_dev <- grDevices::dev.cur()
      on.exit(.dev_close(pdf_dev), add = TRUE)
      .render_plot_to_device(fig)
      invisible(fpath)
    }
  )
)

#' Does a file name end in .pdf (any case)
#' @param name file name
#' @noRd
#' @srrstats {G4.0} PDF writers append `.pdf` only to file names without it (`write_pdf()`, `write_boxplots()`;
#'   test-tidyMS_plotting.R).
.has_pdf_suffix <- function(name) {
  grepl("\\.pdf$", name, ignore.case = TRUE)
}
