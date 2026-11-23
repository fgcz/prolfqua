test_that("plot_hierarchies_boxplot", {
  istar <- sim_lfq_data_peptide_config()
  lfq <- LFQData$new(istar$data, istar$config)

  xnested <- lfq$data_long() |>
    dplyr::group_by(dplyr::across(dplyr::all_of(lfq$relevant_hierarchy_keys()))) |>
    tidyr::nest()

  p <- .plot_hierarchies_boxplot(
    xnested$data[[3]],
    xnested$protein_Id[[3]],
    lfq,
    facet_grid_on = tail(lfq$hierarchy_keys(), 1)
  )
  expect_true("ggplot" %in% class(p))

  p <- .plot_hierarchies_boxplot(xnested$data[[3]], xnested$protein_Id[[3]], lfq)
  expect_true("ggplot" %in% class(p))
  p <- .plot_hierarchies_boxplot(
    xnested$data[[3]],
    xnested$protein_Id[[3]],
    lfq,
    beeswarm = FALSE
  )
  expect_true("ggplot" %in% class(p))
})

test_that("density legend is suppressed automatically for many samples", {
  pdata <- expand.grid(
    sample = paste0("S", seq_len(17)),
    feature = seq_len(20),
    KEEP.OUT.ATTRS = FALSE
  )
  withr::local_seed(1)
  pdata$abundance <- stats::rnorm(nrow(pdata), mean = 10)

  expect_false(prolfqua:::.show_sample_legend(pdata, "sample", NA, 16))
  expect_true(prolfqua:::.show_sample_legend(pdata, "sample", TRUE, 16))
  expect_false(prolfqua:::.show_sample_legend(pdata, "sample", FALSE, 16))

  p <- .plot_intensity_distribution_density(
    pdata,
    sample_name = "sample",
    response = "abundance",
    is_transformed = TRUE
  )
  expect_s3_class(p, "ggplot")
  expect_equal(length(p$guides$guides), 1)
})

test_that("density plot carries sample keys for plotly hover highlighting", {
  testthat::skip_if_not_installed("plotly")

  pdata <- expand.grid(
    sample = c("S1", "S2"),
    feature = seq_len(30),
    KEEP.OUT.ATTRS = FALSE
  )
  withr::local_seed(1)
  pdata$abundance <- c(stats::rnorm(30, mean = 10), stats::rnorm(30, mean = 12))

  p <- .plot_intensity_distribution_density(
    pdata,
    sample_name = "sample",
    response = "abundance",
    is_transformed = TRUE
  )
  built <- plotly::plotly_build(plotly::ggplotly(p))
  trace_keys <- vapply(built$x$data, function(trace) trace$key[[1]], character(1))

  expect_setequal(trace_keys, c("S1", "S2"))
})

test_that("PCA sample labels repel and are not clipped", {
  withr::local_seed(1)
  matrix <- matrix(stats::rnorm(60), nrow = 10)
  colnames(matrix) <- paste0("sample_", seq_len(6))
  annotation <- data.frame(
    sample = colnames(matrix),
    group = rep(c("control", "treated"), each = 3)
  )

  p <- .plot_pca(
    matrix,
    annotation,
    sample_name = "sample",
    factor_keys = "group",
    add_txt = TRUE
  )

  geom_classes <- vapply(
    p$layers,
    function(layer) class(layer$geom)[[1]],
    character(1)
  )

  expect_s3_class(p, "ggplot")
  expect_contains(geom_classes, "GeomTextRepel")
  expect_equal(p$coordinates$clip, "off")
})

test_that("plot_pca errors when no complete features remain", {
  withr::local_seed(1)
  matrix <- matrix(stats::rnorm(60), nrow = 10)
  colnames(matrix) <- paste0("sample_", seq_len(6))
  matrix[, 1] <- NA # every feature now has at least one missing value
  annotation <- data.frame(
    sample = colnames(matrix),
    group = rep(c("control", "treated"), each = 3)
  )

  expect_error(
    suppressMessages(.plot_pca(
      matrix,
      annotation,
      sample_name = "sample",
      factor_keys = "group"
    )),
    "no features without missing values"
  )
})

test_that("plot_pca errors with too few samples instead of returning NULL", {
  withr::local_seed(1)
  matrix <- matrix(stats::rnorm(20), nrow = 10)
  colnames(matrix) <- c("s1", "s2")
  annotation <- data.frame(sample = colnames(matrix), group = c("a", "b"))

  expect_error(
    .plot_pca(matrix, annotation, sample_name = "sample", factor_keys = "group"),
    "at least 3 samples"
  )
})

test_that("plot_pca errors on duplicated sample names", {
  withr::local_seed(1)
  matrix <- matrix(stats::rnorm(60), nrow = 10)
  colnames(matrix) <- c("s1", "s1", "s2", "s3", "s4", "s5")
  annotation <- data.frame(
    sample = c("s1", "s2", "s3", "s4", "s5"),
    group = c("a", "a", "b", "b", "a")
  )

  expect_error(
    .plot_pca(matrix, annotation, sample_name = "sample", factor_keys = "group"),
    "duplicated sample names"
  )
})

test_that("plot_pca joins explicitly and handles a single factor key", {
  withr::local_seed(1)
  matrix <- matrix(stats::rnorm(60), nrow = 10)
  colnames(matrix) <- paste0("sample_", seq_len(6))
  annotation <- data.frame(
    sample = colnames(matrix),
    group = rep(c("control", "treated"), each = 3),
    stringsAsFactors = FALSE
  )

  # explicit by = sample_name => no "Joining, by" message; single factor key
  # (factor_keys[2] is NA) must not add a shape scale.
  expect_no_message(
    p <- .plot_pca(matrix, annotation, sample_name = "sample", factor_keys = "group")
  )
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), ncol(matrix))
  geom_classes <- vapply(p$layers, function(layer) class(layer$geom)[[1]], character(1))
  expect_contains(geom_classes, "GeomPoint")
})

test_that("heatmap row labels are truncated without changing matrix row names", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  withr::local_seed(1)
  matrix <- matrix(stats::rnorm(12), nrow = 3)
  original_labels <- c(
    "Feature001_very_long_metabolomics_identifier_mz101.1234_rt12.34",
    "short",
    "Feature002_very_long_metabolomics_identifier_mz201.1234_rt22.34"
  )
  rownames(matrix) <- original_labels
  colnames(matrix) <- paste0("S", seq_len(4))
  annotation <- data.frame(
    sample = colnames(matrix),
    group = rep(c("control", "treated"), each = 2)
  )

  p <- .plot_raster(
    matrix,
    annotation,
    factor_keys = "group",
    sample_name = "sample",
    show_rownames = TRUE,
    max_rownames_chars = 20
  )

  labels <- p@row_names_param$labels

  expect_s4_class(p, "Heatmap")
  expect_equal(rownames(matrix), original_labels)
  expect_true(any(grepl("[.][.][.]", labels)))
  expect_false(any(original_labels[1] %in% labels))
  expect_true(all(nchar(labels[grepl("[.][.][.]", labels)]) <= 20))
})

test_long_sample_heatmap_data <- function() {
  withr::local_seed(1)
  matrix <- matrix(stats::rnorm(30), nrow = 5)
  sample_names <- paste0(
    "20250225_0",
    seq_len(6),
    "_C37649_S8852",
    seq_len(6),
    "_neg_Dig_14_D3_",
    rep(c("A", "B"), length.out = 6)
  )
  colnames(matrix) <- sample_names
  rownames(matrix) <- paste0("feature", seq_len(5))
  annotation <- data.frame(
    sample = sample_names,
    group = rep(c("control", "treated"), each = 3)
  )
  list(matrix = matrix, sample_names = sample_names, annotation = annotation)
}

test_that("correlation heatmap sample labels keep informative suffixes", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  data <- test_long_sample_heatmap_data()

  p <- .plot_heatmap_cor(
    data$matrix,
    data$annotation,
    factor_keys = "group",
    sample_name = "sample",
    max_sample_label_chars = 20
  )

  labels <- p@column_names_param$labels
  suffixes <- prolfqua:::.suffix_plot_labels(data$sample_names, 20)

  expect_s4_class(p, "Heatmap")
  expect_true(all(suffixes %in% labels))
  expect_false(any(data$sample_names %in% labels))
  expect_true(all(nchar(suffixes) <= 20))
})

test_that("abundance heatmap sample labels keep informative suffixes", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  data <- test_long_sample_heatmap_data()

  p <- .plot_heatmap(
    data$matrix,
    data$annotation,
    factor_keys = "group",
    sample_name = "sample",
    max_sample_label_chars = 20
  )

  labels <- p@column_names_param$labels
  suffixes <- prolfqua:::.suffix_plot_labels(data$sample_names, 20)

  expect_s4_class(p, "Heatmap")
  expect_true(all(suffixes %in% labels))
  expect_false(any(data$sample_names %in% labels))
  expect_true(all(nchar(suffixes) <= 20))
})

test_that("raster sample labels keep informative suffixes", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  data <- test_long_sample_heatmap_data()

  p <- .plot_raster(
    data$matrix,
    data$annotation,
    factor_keys = "group",
    sample_name = "sample",
    max_sample_label_chars = 20
  )

  labels <- p@column_names_param$labels
  suffixes <- prolfqua:::.suffix_plot_labels(data$sample_names, 20)

  expect_s4_class(p, "Heatmap")
  expect_true(all(suffixes %in% labels))
  expect_false(any(data$sample_names %in% labels))
  expect_true(all(nchar(suffixes) <= 20))
})

test_that("abundance heatmap shows missing values in light gray", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  data <- test_long_sample_heatmap_data()
  data$matrix[1, 1] <- NA

  p <- .plot_heatmap(
    data$matrix,
    data$annotation,
    factor_keys = "group",
    sample_name = "sample"
  )

  expect_s4_class(p, "Heatmap")
  expect_equal(
    grDevices::col2rgb(p@matrix_color_mapping@na_col),
    grDevices::col2rgb(prolfqua:::.HEATMAP_NA_COL)
  )
})

test_that("abundance heatmap draws sparse one-row matrices with all-NA columns", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  matrix <- matrix(
    c(NA, 12, 13, NA, 12.5, 13.5),
    nrow = 1,
    dimnames = list("P1", c("C1", "T1", "T2", "C2", "T3", "T4"))
  )
  annotation <- data.frame(
    sample = colnames(matrix),
    group = c("control", "treated", "treated", "control", "treated", "treated")
  )

  p <- .plot_heatmap(
    matrix,
    annotation,
    factor_keys = "group",
    sample_name = "sample",
    show_rownames = TRUE
  )

  expect_s4_class(p, "Heatmap")
  withr::local_pdf(NULL)
  expect_no_error(ComplexHeatmap::draw(p))
})

test_that("abundance heatmap keeps partially observed rows and drops rows that cannot be z-scored", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  matrix <- rbind(
    complete = c(11, 12, 13, 14, 8, 9, 10, 11),
    bait_only = c(14, 15, 16, 17, NA, NA, NA, NA),
    constant = rep(10, 8),
    one_observation = c(12, rep(NA, 7)),
    all_missing = rep(NA, 8)
  )
  colnames(matrix) <- c(paste0("B", 1:4), paste0("C", 1:4))
  annotation <- data.frame(
    sample = colnames(matrix),
    group = rep(c("bait", "control"), each = 4)
  )

  p <- .plot_heatmap(
    matrix,
    annotation,
    factor_keys = "group",
    sample_name = "sample",
    na_fraction = 1,
    show_rownames = TRUE
  )

  expect_s4_class(p, "Heatmap")
  expect_setequal(rownames(p@matrix), c("complete", "bait_only"))
  withr::local_pdf(NULL)
  expect_no_warning(ComplexHeatmap::draw(p))
})

test_that("LFQDataPlotter heatmaps retain all z-scoreable rows by default", {
  expect_identical(formals(LFQDataPlotter$public_methods$heatmap)$na_fraction, 1)
})

#' @srrstats {EA5.0, EA5.0b} Default colours are colour-blind safe: blue-white-red abundance heatmaps, Okabe-Ito (or
#'   viridis beyond 8 levels) for factor levels, black and Okabe-Ito for estimate types; tested here and in
#'   test-ContrastsPlotter.R.
#' @noRd
NULL

test_that("abundance colour mapping spans blue, white and red (no red-green contrast)", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  col_fun <- prolfqua:::.abundance_col_fun(matrix(c(-2, 0, 2), nrow = 1))
  cols <- grDevices::col2rgb(col_fun(c(-2, 0, 2)))
  expect_true(cols["blue", 1] > cols["red", 1] && cols["blue", 1] > cols["green", 1])
  expect_true(all(cols[, 2] > 230))
  expect_true(cols["red", 3] > cols["blue", 3] && cols["red", 3] > cols["green", 3])
  expect_true(cols["green", 1] < 150 && cols["green", 3] < 150)
})

test_that("factor colours are Okabe-Ito up to 8 levels and viridis beyond", {
  okabe_ito <- unname(grDevices::palette.colors(palette = "Okabe-Ito"))[-1L]
  expect_identical(prolfqua:::.factor_palette(3), okabe_ito[1:3])
  expect_identical(prolfqua:::.factor_palette(8), okabe_ito)
  expect_identical(prolfqua:::.factor_palette(9), grDevices::hcl.colors(9, "viridis"))
  expect_length(prolfqua:::.factor_palette(0), 0L)
})

test_that("group colours are deterministic and shared across heatmap, na-heatmap, and PCA", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  # Regression: ComplexHeatmap assigned random colours to discrete annotations
  # (no `col` set), so the heatmap, na-heatmap, and PCA legends each used a
  # different colour for the same group. They must now share one palette.
  istar <- sim_lfq_data_protein_config()
  lfq <- LFQData$new(istar$data, istar$config)
  wide <- lfq$data_wide(as.matrix = TRUE)
  fk <- lfq$factor_keys()
  sn <- lfq$sample_name()
  key <- fk[1]

  # deterministic: same annotation -> same colours on every call
  ta1 <- prolfqua:::.heatmap_top_annotation(wide$annotation, fk, sn, colnames(wide$data))
  ta2 <- prolfqua:::.heatmap_top_annotation(wide$annotation, fk, sn, colnames(wide$data))
  heat1 <- ta1@anno_list[[key]]@color_mapping@colors
  heat2 <- ta2@anno_list[[key]]@color_mapping@colors
  expect_identical(heat1, heat2)

  # na-heatmap shares the same annotation colours
  na_h <- .plot_na_heatmap(wide$data, wide$annotation, fk, sn)
  na_colors <- na_h@top_annotation@anno_list[[key]]@color_mapping@colors
  expect_identical(na_colors[order(names(na_colors))], heat1[order(names(heat1))])

  # PCA colours each group with the same colour as the heatmap annotation
  pca <- .plot_pca(wide$data, wide$annotation, sn, fk)
  build <- ggplot2::ggplot_build(pca)
  pca_map <- tapply(build$data[[1]]$colour, as.character(build$plot$data[[key]]), unique)
  for (lev in names(pca_map)) {
    expect_equal(
      grDevices::col2rgb(pca_map[[lev]]),
      grDevices::col2rgb(heat1[[lev]])
    )
  }
})

test_that("write_pdf renders a ComplexHeatmap to a non-empty PDF", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  istar <- sim_lfq_data_peptide_config()
  lfq <- LFQData$new(istar$data, istar$config)
  pl <- lfq$get_Plotter()

  out_dir <- withr::local_tempdir()
  fpath <- pl$write_pdf(pl$heatmap(), out_dir, "heatmap_smoke")

  expect_true(file.exists(fpath))
  expect_gt(file.info(fpath)$size, 0)
})

test_that("volcano_plotly does not cap small non-zero FDR values by default", {
  skip_if_not_installed("plotly")
  skip_if_not_installed("crosstalk")
  data <- data.frame(
    fc = c(-2, 0, 2),
    BFDR = c(1e-8, 1e-3, 0),
    condition = "A",
    Prey = c("P1", "P2", "P3"),
    modelName = "model"
  )

  default_plot <- volcano_plotly(data)[[1]]
  default_y <- unlist(lapply(plotly::plotly_build(default_plot)$x$data, function(trace) trace$y))

  expect_gt(max(default_y[is.finite(default_y)]), 7.9)
  expect_true(all(is.finite(default_y)))

  capped_plot <- volcano_plotly(data, minsignificance = 1e-4)[[1]]
  capped_y <- unlist(lapply(plotly::plotly_build(capped_plot)$x$data, function(trace) trace$y))

  expect_equal(max(capped_y), 4)
})

test_that("heatmap keeps only the top_n most variable features", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  istar <- sim_lfq_data_protein_config()
  lfq <- LFQData$new(istar$data, istar$config)
  wide <- lfq$data_wide(as.matrix = TRUE)
  n_all <- nrow(wide$data)
  expect_gt(n_all, 3)

  # n_most_variable >= feature count (or NULL / Inf) keeps every feature.
  expect_equal(nrow(prolfqua:::.select_most_variable_features(lfq, wide, n_most_variable = n_all + 10)), n_all)
  expect_equal(nrow(prolfqua:::.select_most_variable_features(lfq, wide, n_most_variable = NULL)), n_all)
  expect_equal(nrow(prolfqua:::.select_most_variable_features(lfq, wide, n_most_variable = Inf)), n_all)

  # n_most_variable < feature count keeps exactly n_most_variable rows, and they are the ones with
  # the highest variability statistic (CV here — untransformed data).
  k <- 3L
  sub <- prolfqua:::.select_most_variable_features(lfq, wide, n_most_variable = k)
  expect_equal(nrow(sub), k)

  st <- lfq$get_Stats(stats = "all")
  keys <- lfq$hierarchy_keys()
  ranked <- dplyr::left_join(
    wide$rowdata,
    dplyr::select(st$stats(), dplyr::all_of(c(keys, st$stat))),
    by = keys
  )
  metric <- stats::setNames(ranked[[st$stat]], rownames(wide$data))
  kept <- metric[rownames(sub)]
  dropped <- metric[setdiff(rownames(wide$data), rownames(sub))]
  expect_gte(min(kept, na.rm = TRUE), max(dropped, na.rm = TRUE))

  # end-to-end through the plotter: no hclust blow-up, returns a Heatmap.
  p <- lfq$get_Plotter()$heatmap(n_most_variable = k)
  expect_s4_class(p, "Heatmap")
  expect_warning(p <- lfq$get_Plotter()$heatmap(top_n = k), "n_most_variable", class = "deprecatedWarning")
  expect_equal(nrow(p@matrix), k)
})

test_that("heatmap row_order shows exactly the given features in that order, without clustering", {
  istar <- sim_lfq_data_protein_config(n_proteins = 20)
  lfq <- LFQData$new(istar$data, istar$config)
  plotter <- lfq$get_Plotter()
  ids <- rownames(lfq$data_wide(as.matrix = TRUE)$data)
  wanted <- rev(ids[c(2, 5, 9, 11)])
  p <- plotter$heatmap(row_order = wanted)
  expect_s4_class(p, "Heatmap")
  expect_identical(rownames(p@matrix), wanted)

  expect_error(
    plotter$heatmap(row_order = c(wanted, "not_a_feature")),
    "not_a_feature",
    class = "prolfqua_error_bad_argument"
  )
  expect_error(plotter$heatmap(row_order = wanted, n_most_variable = 2), class = "prolfqua_error_bad_argument")
})

# Graphics devices: PDF writers must leave the caller's open device open and current.
.open_caller_device <- function(env = parent.frame()) {
  grDevices::pdf(tempfile(fileext = ".pdf"))
  own <- grDevices::dev.cur()
  withr::defer(if (own %in% grDevices::dev.list()) grDevices::dev.off(own), envir = env)
  list(cur = own, list = grDevices::dev.list())
}

test_that("write_pdf keeps the caller's graphics devices", {
  istar <- sim_lfq_data_peptide_config(n_proteins = 5)
  pl <- LFQData$new(istar$data, istar$config)$get_Plotter()
  out_dir <- withr::local_tempdir()
  before <- .open_caller_device()
  fpath <- pl$write_pdf(ggplot2::ggplot(), out_dir, "device_test")
  expect_true(file.exists(fpath))
  expect_equal(grDevices::dev.cur(), before$cur)
  expect_equal(grDevices::dev.list(), before$list)
})

test_that("write_boxplots and write_plots close their device when rendering fails", {
  local_mocked_bindings(.render_plot_to_device = function(x) stop("boom"))
  istar <- sim_lfq_data_peptide_config(n_proteins = 5)
  lfq <- LFQData$new(istar$data, istar$config)
  out_dir <- withr::local_tempdir()
  agg <- lfq$get_Transformer()$log2()$lfq$get_Aggregator("medpolish")
  agg$aggregate()
  before <- .open_caller_device()

  expect_error(lfq$get_Plotter()$write_boxplots(out_dir), "boom")
  expect_equal(grDevices::dev.list(), before$list)
  expect_equal(grDevices::dev.cur(), before$cur)

  expect_error(agg$write_plots(out_dir), "boom")
  expect_equal(grDevices::dev.list(), before$list)
  expect_equal(grDevices::dev.cur(), before$cur)
})

test_that("lr_test with a path keeps the caller's graphics devices", {
  d <- sim_lfq_data_2factor_config(n_proteins = 10, with_missing = FALSE)
  lfq <- LFQData$new(d$data, d$config)
  mod_full <- build_model(lfq, strategy_lm("abundance ~ Treatment + Background"))
  mod_reduced <- build_model(lfq, strategy_lm("abundance ~ Treatment"))
  out_dir <- withr::local_tempdir()
  before <- .open_caller_device()
  lr_test(mod_full$model_df, "full", mod_reduced$model_df, "reduced", path = out_dir)
  expect_true(length(list.files(out_dir, pattern = "^hist_LRT_")) == 1)
  expect_equal(grDevices::dev.cur(), before$cur)
  expect_equal(grDevices::dev.list(), before$list)
})

test_that("plot_lrt_diagnostics closes its device when plotting fails", {
  out_dir <- withr::local_tempdir()
  before <- .open_caller_device()
  bad <- data.frame(likelihood_ratio_test.pValue = "not a number")
  expect_error(plot_lrt_diagnostics(bad, "full", "reduced", out_dir), "numeric")
  expect_equal(grDevices::dev.cur(), before$cur)
  expect_equal(grDevices::dev.list(), before$list)
})

test_that("default density and box plots work without plotly and ggbeeswarm", {
  local_mocked_bindings(.has_package = function(pkg) FALSE)
  istar <- sim_lfq_data_peptide_config(n_proteins = 10, with_missing = FALSE, seed = 1)
  lfq <- LFQData$new(istar$data, istar$config)
  pl <- lfq$get_Plotter()
  density <- pl$intensity_distribution_density()
  expect_s3_class(density, "ggplot")
  expect_no_error(ggplot2::ggplot_build(density))
  boxplot <- pl$boxplots()$boxplot[[1]]
  expect_no_error(ggplot2::ggplot_build(boxplot))
  expect_true(any(vapply(boxplot$layers, function(l) inherits(l$position, "PositionJitterdodge"), logical(1))))
})

#' @srrstats {EA5.0a, EA5.1} Text sizes follow the ggplot2 and corrplot defaults; the remaining override, volcano
#'   labels, defaults to 3 mm (about 8.5 pt) instead of 1 mm and stays an argument, as do corrplot `tl.cex` and
#'   `number.cex` through `...`.
#' @noRd
NULL

test_that("volcano labels and correlation plot text use readable default sizes", {
  dat <- data.frame(
    fc = c(-3, 3, 0.1),
    p.adjust = c(0.001, 0.001, 0.9),
    condition = "a",
    colour = "observed",
    name = c("x", "y", "z")
  )
  p <- .multigroup_volcano_labelled(dat, colour = NULL, label = "name")
  text_layer <- Filter(function(l) inherits(l$geom, "GeomTextRepel"), p$layers)
  expect_length(text_layer, 1L)
  expect_identical(text_layer[[1]]$aes_params$size, 3)
  expect_false(any(c("tl.cex", "number.cex") %in% names(formals(.plot_sample_correlation))))
  expect_true("..." %in% names(formals(.plot_sample_correlation)))
})

test_that("PDF writers append .pdf only when the name lacks it", {
  withr::local_options(prolfqua.verbose = FALSE)
  istar <- sim_lfq_data_protein_config(n_proteins = 3, with_missing = FALSE)
  pl <- LFQData$new(istar$data, istar$config)$get_Plotter()
  out_dir <- withr::local_tempdir()
  expect_identical(basename(pl$write_pdf(ggplot2::ggplot(), out_dir, "fig")), paste0(pl$prefix, "fig.pdf"))
  expect_identical(basename(pl$write_pdf(ggplot2::ggplot(), out_dir, "fig.pdf")), paste0(pl$prefix, "fig.pdf"))
  expect_identical(basename(pl$write_pdf(ggplot2::ggplot(), out_dir, "fig.PDF")), paste0(pl$prefix, "fig.PDF"))
  expect_identical(basename(pl$write_boxplots(out_dir, "qc_")), "qc_boxplot.pdf")
  expect_identical(basename(pl$write_boxplots(out_dir, "boxes.pdf")), "boxes.pdf")
  expect_false(any(grepl("\\.pdf\\.pdf$", list.files(out_dir), ignore.case = TRUE)))
})

test_that("sample correlations accept a cor() method or a function", {
  m <- matrix(c(1, 2, 3, 4, 2, 4, 5, 9, 4, 3, 2, 1), ncol = 3)
  expect_equal(prolfqua:::.sample_correlation(m, "pearson"), stats::cor(m))
  expect_equal(prolfqua:::.sample_correlation(m, "spearman"), stats::cor(m, method = "spearman"))
  custom <- function(x) stats::cor(x, method = "kendall")
  expect_equal(prolfqua:::.sample_correlation(m, custom), stats::cor(m, method = "kendall"))
  expect_error(prolfqua:::.sample_correlation(m, "robust"), class = "prolfqua_error_bad_argument")
  expect_error(prolfqua:::.sample_correlation(m, function(x) 1), class = "prolfqua_error_bad_argument")
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  istar <- sim_lfq_data_protein_config(n_proteins = 10, with_missing = FALSE)
  pl <- LFQData$new(istar$data, istar$config)$get_Plotter()
  expect_true(methods::is(pl$heatmap_cor(method = "spearman"), "Heatmap"))
})

test_that("plot_pca draws only valid shapes for a second factor with many levels (review B10)", {
  n <- 30
  annot <- data.frame(
    sample = paste0("s", seq_len(n)),
    group = rep(c("A", "B"), n / 2),
    subject = paste0("p", seq_len(n))
  )
  m <- withr::with_seed(1, matrix(rnorm(20 * n), nrow = 20, dimnames = list(NULL, annot$sample)))
  p <- .plot_pca(m, annot, "sample", c("group", "subject"))
  shapes <- ggplot2::layer_data(p, 1)$shape
  expect_true(is.null(shapes) || all(shapes %in% 0:25))

  few <- annot[1:6, ]
  few$subject <- rep(c("p1", "p2", "p3"), 2)
  p <- .plot_pca(m[, few$sample], few, "sample", c("group", "subject"))
  expect_setequal(unique(ggplot2::layer_data(p, 1)$shape), c(16, 17, 15))
})

test_that("plot_heatmap_cor maps `color` from the lowest to the highest correlation", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  data <- test_long_sample_heatmap_data()
  heatmap <- function(...) {
    .plot_heatmap_cor(data$matrix, data$annotation, factor_keys = "group", sample_name = "sample", ...)
  }
  rgb_at <- function(p, value) grDevices::col2rgb(p@matrix_color_mapping@col_fun(value))
  p <- heatmap(color = c("blue", "yellow"))
  expect_equal(rgb_at(p, min(p@matrix)), grDevices::col2rgb("blue"))
  expect_equal(rgb_at(p, max(p@matrix)), grDevices::col2rgb("yellow"))
  default <- heatmap()
  expect_equal(rgb_at(default, min(default@matrix)), grDevices::col2rgb("white"))
  expect_equal(rgb_at(default, max(default@matrix)), grDevices::col2rgb("red"))
})

test_that("pca_plotly plots the requested principal components", {
  skip_if_not_installed("plotly")
  istar <- sim_lfq_data_protein_config(with_missing = FALSE)
  lfq <- LFQData$new(istar$data, istar$config)
  layout <- plotly::plotly_build(lfq$get_Plotter()$pca_plotly(PC = c(2, 3)))$x$layout
  expect_match(layout$xaxis$title$text, "PC2")
  expect_match(layout$yaxis$title$text, "PC3")
})
