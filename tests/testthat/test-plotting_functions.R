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
