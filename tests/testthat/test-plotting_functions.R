test_that("plot_hierarchies_boxplot", {
  istar <- sim_lfq_data_peptide_config()
  lfq <- LFQData$new(istar$data, istar$config)

  xnested <- lfq$data_long() |>
    dplyr::group_by_at(lfq$relevant_hierarchy_keys()) |>
    tidyr::nest()

  p <- plot_hierarchies_boxplot(
    xnested$data[[3]],
    xnested$protein_Id[[3]],
    lfq,
    facet_grid_on = tail(lfq$hierarchy_keys(), 1)
  )
  expect_true("ggplot" %in% class(p))

  p <- plot_hierarchies_boxplot(xnested$data[[3]], xnested$protein_Id[[3]], lfq)
  expect_true("ggplot" %in% class(p))
  p <- plot_hierarchies_boxplot(
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
  pdata$abundance <- stats::rnorm(nrow(pdata), mean = 10)

  expect_false(prolfqua:::.show_sample_legend(pdata, "sample", NA, 16))
  expect_true(prolfqua:::.show_sample_legend(pdata, "sample", TRUE, 16))
  expect_false(prolfqua:::.show_sample_legend(pdata, "sample", FALSE, 16))

  p <- plot_intensity_distribution_density(
    pdata,
    sample_name = "sample",
    response = "abundance",
    is_transformed = TRUE
  )
  expect_s3_class(p, "ggplot")
  expect_equal(length(p$guides$guides), 1)
})

test_that("PCA sample labels repel and are not clipped", {
  matrix <- matrix(stats::rnorm(60), nrow = 10)
  colnames(matrix) <- paste0("sample_", seq_len(6))
  annotation <- data.frame(
    sample = colnames(matrix),
    group = rep(c("control", "treated"), each = 3)
  )

  p <- plot_pca(
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
  matrix <- matrix(stats::rnorm(60), nrow = 10)
  colnames(matrix) <- paste0("sample_", seq_len(6))
  matrix[, 1] <- NA # every feature now has at least one missing value
  annotation <- data.frame(
    sample = colnames(matrix),
    group = rep(c("control", "treated"), each = 3)
  )

  expect_error(
    suppressMessages(plot_pca(
      matrix,
      annotation,
      sample_name = "sample",
      factor_keys = "group"
    )),
    "no features without missing values"
  )
})

test_that("plot_pca errors with too few samples instead of returning NULL", {
  matrix <- matrix(stats::rnorm(20), nrow = 10)
  colnames(matrix) <- c("s1", "s2")
  annotation <- data.frame(sample = colnames(matrix), group = c("a", "b"))

  expect_error(
    plot_pca(matrix, annotation, sample_name = "sample", factor_keys = "group"),
    "at least 3 samples"
  )
})

test_that("plot_pca errors on duplicated sample names", {
  matrix <- matrix(stats::rnorm(60), nrow = 10)
  colnames(matrix) <- c("s1", "s1", "s2", "s3", "s4", "s5")
  annotation <- data.frame(
    sample = c("s1", "s2", "s3", "s4", "s5"),
    group = c("a", "a", "b", "b", "a")
  )

  expect_error(
    plot_pca(matrix, annotation, sample_name = "sample", factor_keys = "group"),
    "duplicated sample names"
  )
})

test_that("plot_pca joins explicitly and handles a single factor key", {
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
    p <- plot_pca(matrix, annotation, sample_name = "sample", factor_keys = "group")
  )
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), ncol(matrix))
  geom_classes <- vapply(p$layers, function(layer) class(layer$geom)[[1]], character(1))
  expect_contains(geom_classes, "GeomPoint")
})

test_that("heatmap row labels are truncated without changing matrix row names", {
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

  p <- plot_raster(
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
  data <- test_long_sample_heatmap_data()

  p <- plot_heatmap_cor(
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
  data <- test_long_sample_heatmap_data()

  p <- plot_heatmap(
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
  data <- test_long_sample_heatmap_data()

  p <- plot_raster(
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
  data <- test_long_sample_heatmap_data()
  data$matrix[1, 1] <- NA

  p <- plot_heatmap(
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
  matrix <- matrix(
    c(NA, 12, 13, NA, 12.5, 13.5),
    nrow = 1,
    dimnames = list("P1", c("C1", "T1", "T2", "C2", "T3", "T4"))
  )
  annotation <- data.frame(
    sample = colnames(matrix),
    group = c("control", "treated", "treated", "control", "treated", "treated")
  )

  p <- plot_heatmap(
    matrix,
    annotation,
    factor_keys = "group",
    sample_name = "sample",
    show_rownames = TRUE
  )

  expect_s4_class(p, "Heatmap")
  expect_no_error(ComplexHeatmap::draw(p))
})

test_that("abundance color mapping spans green, black, and red", {
  col_fun <- prolfqua:::.abundance_col_fun(matrix(c(-2, 0, 2), nrow = 1))
  cols <- grDevices::col2rgb(col_fun(c(-2, 0, 2)))

  # min -> green, center -> black, max -> red
  expect_true(cols["green", 1] > 200 && cols["red", 1] < 50 && cols["blue", 1] < 50)
  expect_true(all(cols[, 2] < 50))
  expect_true(cols["red", 3] > 200 && cols["green", 3] < 50 && cols["blue", 3] < 50)
})

test_that("group colours are deterministic and shared across heatmap, na-heatmap, and PCA", {
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
  na_h <- plot_na_heatmap(wide$data, wide$annotation, fk, sn)
  na_colors <- na_h@top_annotation@anno_list[[key]]@color_mapping@colors
  expect_identical(na_colors[order(names(na_colors))], heat1[order(names(heat1))])

  # PCA colours each group with the same colour as the heatmap annotation
  pca <- plot_pca(wide$data, wide$annotation, sn, fk)
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
  istar <- sim_lfq_data_peptide_config()
  lfq <- LFQData$new(istar$data, istar$config)
  pl <- lfq$get_Plotter()

  out_dir <- withr::local_tempdir()
  fpath <- pl$write_pdf(pl$heatmap(), out_dir, "heatmap_smoke")

  expect_true(file.exists(fpath))
  expect_gt(file.info(fpath)$size, 0)
})

test_that("volcano_plotly does not cap small non-zero FDR values by default", {
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
  istar <- sim_lfq_data_protein_config()
  lfq <- LFQData$new(istar$data, istar$config)
  wide <- lfq$data_wide(as.matrix = TRUE)
  n_all <- nrow(wide$data)
  expect_gt(n_all, 3)

  # top_n >= feature count (or NULL / Inf) keeps every feature.
  expect_equal(nrow(prolfqua:::.select_most_variable_features(lfq, wide, top_n = n_all + 10)), n_all)
  expect_equal(nrow(prolfqua:::.select_most_variable_features(lfq, wide, top_n = NULL)), n_all)
  expect_equal(nrow(prolfqua:::.select_most_variable_features(lfq, wide, top_n = Inf)), n_all)

  # top_n < feature count keeps exactly top_n rows, and they are the ones with
  # the highest variability statistic (CV here — untransformed data).
  k <- 3L
  sub <- prolfqua:::.select_most_variable_features(lfq, wide, top_n = k)
  expect_equal(nrow(sub), k)

  st <- lfq$get_Stats(stats = "all")
  keys <- c(lfq$hierarchy_keys(), lfq$isotope_label())
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
  p <- lfq$get_Plotter()$heatmap(top_n = k)
  expect_s4_class(p, "Heatmap")
})
