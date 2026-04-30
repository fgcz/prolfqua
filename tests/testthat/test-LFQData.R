# Test LFQData R6 class and decorator classes

# Shared setup
istar <- sim_lfq_data_peptide_config()
lfqdata <- LFQData$new(istar$data, istar$config)

test_that("LFQData creation and basic accessors", {
  expect_s3_class(lfqdata, "LFQData")
  expect_true("data.frame" %in% class(lfqdata$data_long()))
  expect_s3_class(lfqdata$get_config(), "AnalysisConfiguration")
  expect_type(lfqdata$subject_id(), "character")
  expect_s3_class(lfqdata$hierarchy(), "data.frame")
  expect_s3_class(lfqdata$factors(), "data.frame")
  expect_type(lfqdata$response(), "character")
  expect_false(lfqdata$is_transformed())
})

test_that("LFQData filtering and subsetting", {
  lfq_copy <- lfqdata$get_copy()
  orig_hier <- lfq_copy$hierarchy()
  lfq_copy$filter_proteins_by_peptide_count()
  # After filtering, should have same or fewer hierarchy entries
  expect_true(nrow(lfq_copy$hierarchy()) <= nrow(orig_hier))

  f1 <- lfqdata$omit_na(nr_na = 0)
  expect_true(nrow(f1$hierarchy()) <= nrow(lfqdata$hierarchy()))

  f2 <- lfqdata$omit_na(factor_depth = 0)
  expect_true(nrow(f2$hierarchy()) <= nrow(lfqdata$hierarchy()))

  cc <- lfqdata$get_copy()
  cc$complete_cases()
  expect_s3_class(cc, "LFQData")

  rsi <- lfqdata$get_copy()
  rsi$remove_small_intensities()
  expect_s3_class(rsi, "LFQData")

  samp <- lfqdata$get_sample(5, seed = 4)
  expect_s3_class(samp, "LFQData")
  expect_equal(nrow(samp$hierarchy()), 5)

  sub <- lfqdata$get_subset(head(lfqdata$hierarchy(), 3))
  expect_s3_class(sub, "LFQData")

  lfq2 <- lfqdata$get_copy()
  lfq2$set_data(lfq2$data_long()[1:100, ])
  res <- lfqdata$filter_difference(lfq2)
  expect_equal(nrow(res$data_long()), nrow(lfqdata$data_long()) - 100)
})

test_that("LFQData data_wide conversion", {
  tmp <- lfqdata$data_wide()
  expect_equal(nrow(tmp$data), nrow(tmp$rowdata))
  expect_equal(ncol(tmp$data), nrow(tmp$annotation) + ncol(tmp$rowdata))
  expect_true("data.frame" %in% class(tmp$data))

  tmp_mat <- lfqdata$data_wide(as.matrix = TRUE)
  expect_true("matrix" %in% class(tmp_mat$data))
})

test_that("LFQData get_copy and rename", {
  copy <- lfqdata$get_copy()
  expect_s3_class(copy, "LFQData")
  # Verify it's independent
  old_response <- copy$response()
  copy$rename_response("test_intensity")
  expect_equal(copy$response(), "test_intensity")
  # Original unchanged
  expect_equal(lfqdata$response(), old_response)

  hc <- lfqdata$hierarchy_counts()
  expect_s3_class(hc, "data.frame")

  sh <- lfqdata$summarize_hierarchy()
  expect_s3_class(sh, "data.frame")
})

test_that("LFQData decorator factory methods", {
  expect_s3_class(lfqdata$get_Transformer(), "LFQDataTransformer")
  expect_s3_class(lfqdata$get_Stats(), "LFQDataStats")
  expect_s3_class(lfqdata$get_Summariser(), "LFQDataSummariser")
  expect_s3_class(lfqdata$get_Plotter(), "LFQDataPlotter")
  expect_s3_class(lfqdata$get_Aggregator("medpolish"), "AggregateMedpolish")
})

test_that("LFQDataTransformer", {
  lfq_copy <- lfqdata$get_copy()
  tr <- lfq_copy$get_Transformer()
  expect_false(tr$lfq$is_transformed())

  tr$log2()
  expect_true(tr$lfq$is_transformed())

  scales <- tr$get_scales()
  expect_true("medians" %in% names(scales))
  expect_true("mads" %in% names(scales))

  before_scales <- tr$get_scales()
  tr$robscale()
  after_scales <- tr$get_scales()
  # robscale with preserveMean=TRUE keeps the original mean
  expect_true(abs(mean(before_scales$medians) - mean(after_scales$medians)) < 1e-8)

  # Test intensity_array
  lfq_copy2 <- lfqdata$get_copy()
  tr2 <- lfq_copy2$get_Transformer()
  tr2$intensity_array(log2)
  expect_true(tr2$lfq$is_transformed())
})

test_that("LFQDataStats", {
  lfq_copy <- lfqdata$get_copy()
  stats <- lfq_copy$get_Stats()
  expect_equal(stats$stat, "CV")

  sdf <- stats$stats()
  expect_s3_class(sdf, "data.frame")
  expect_true("sd" %in% colnames(sdf) || "CV" %in% colnames(sdf))

  sw <- stats$stats_wide()
  expect_s3_class(sw, "data.frame")

  sq <- stats$stats_quantiles()
  expect_true(all(c("long", "wide") %in% names(sq)))

  # With transformed data, stat should be "sd"
  lfq_t <- lfqdata$get_copy()
  lfq_t$is_transformed(TRUE)
  stats_t <- lfq_t$get_Stats()
  expect_equal(stats_t$stat, "sd")
})

test_that("AggregateMedpolish", {
  lfq_copy <- lfqdata$get_copy()
  tr <- lfq_copy$get_Transformer()
  tr$log2()
  tr$robscale()

  agg <- AggregateMedpolish$new(tr$lfq, "protein")
  result <- agg$aggregate()
  expect_s3_class(result, "LFQData")
  expect_true(nrow(result$data_long()) < nrow(lfqdata$data_long()))

  pl <- agg$plot()
  expect_true("plots" %in% names(pl))
})

test_that("AggregateRlm", {
  lfq_copy <- lfqdata$get_copy()
  tr <- lfq_copy$get_Transformer()
  tr$log2()
  tr$robscale()

  agg <- AggregateRlm$new(tr$lfq$get_copy(), "protein")
  result <- agg$aggregate()
  expect_s3_class(result, "LFQData")
})

test_that("AggregateTopN", {
  lfq_raw <- lfqdata$get_copy()
  agg <- AggregateTopN$new(lfq_raw, "protein", N = 3, func = "sum")
  result <- agg$aggregate()
  expect_s3_class(result, "LFQData")

  pl <- agg$plot()
  expect_true("plots" %in% names(pl))

  agg_mean <- AggregateTopN$new(lfqdata$get_copy(), "protein", func = "mean")
  result_mean <- agg_mean$aggregate()
  expect_s3_class(result_mean, "LFQData")
})
