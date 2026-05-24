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

  text_grobs <- Filter(function(grob) inherits(grob, "text"), p$gtable$grobs)
  labels <- unlist(lapply(text_grobs, function(grob) grob$label), use.names = FALSE)

  expect_s3_class(p, "pheatmap")
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

  text_grobs <- Filter(function(grob) inherits(grob, "text"), p$gtable$grobs)
  labels <- unlist(lapply(text_grobs, function(grob) grob$label), use.names = FALSE)
  suffixes <- prolfqua:::.suffix_plot_labels(data$sample_names, 20)

  expect_s3_class(p, "pheatmap")
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

  text_grobs <- Filter(function(grob) inherits(grob, "text"), p$gtable$grobs)
  labels <- unlist(lapply(text_grobs, function(grob) grob$label), use.names = FALSE)
  suffixes <- prolfqua:::.suffix_plot_labels(data$sample_names, 20)

  expect_s3_class(p, "pheatmap")
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

  text_grobs <- Filter(function(grob) inherits(grob, "text"), p$gtable$grobs)
  labels <- unlist(lapply(text_grobs, function(grob) grob$label), use.names = FALSE)
  suffixes <- prolfqua:::.suffix_plot_labels(data$sample_names, 20)

  expect_s3_class(p, "pheatmap")
  expect_true(all(suffixes %in% labels))
  expect_false(any(data$sample_names %in% labels))
  expect_true(all(nchar(suffixes) <= 20))
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
