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
})

test_that("filter_proteins_by_peptide_count counts the second hierarchy level within the first", {
  kept_proteins <- function(lfq) sort(unique(lfq$data_long()$protein_Id))

  by_default <- lfqdata$get_copy()$filter_proteins_by_peptide_count()
  explicit <- lfqdata$get_copy()$filter_proteins_by_peptide_count(
    threshold = 2,
    parent = "protein_Id",
    child = "peptide_Id"
  )
  expect_equal(by_default$data_long(), explicit$data_long())
  peptides_per_protein <- dplyr::summarize(
    by_default$data_long(),
    n = dplyr::n_distinct(peptide_Id),
    .by = "protein_Id"
  )
  expect_true(all(peptides_per_protein$n >= 2))

  # peptide-level modelling (hierarchy_depth 2) used to make this a no-op
  depth2 <- lfqdata$get_copy()
  depth2$set_config_value("hierarchy_depth", 2)
  depth2$filter_proteins_by_peptide_count()
  expect_equal(kept_proteins(depth2), kept_proteins(by_default))

  stricter <- lfqdata$get_copy()$filter_proteins_by_peptide_count(threshold = 3)
  expect_true(all(kept_proteins(stricter) %in% kept_proteins(by_default)))

  prot <- sim_lfq_data_protein_config()
  lfq_prot <- LFQData$new(prot$data, prot$config)
  expect_error(lfq_prot$filter_proteins_by_peptide_count(), class = "prolfqua_error_bad_argument")
  expect_error(
    lfqdata$get_copy()$filter_proteins_by_peptide_count(child = "no_such_column"),
    class = "prolfqua_error_missing_column"
  )

  expect_gt(nrow(prolfqua_data("data_ionstar")$filtered()$data), 0)
})

test_that("rename_response and the is_transformed setter return the object", {
  lfq <- lfqdata$get_copy()
  expect_identical(lfq$is_transformed(TRUE), lfq)
  expect_true(lfq$is_transformed())

  old <- lfq$response()
  expect_identical(lfq$rename_response(old), lfq)
  expect_equal(lfq$response(), old)

  expect_error(lfq$rename_response(lfq$file_name()), class = "prolfqua_error_bad_argument")
  expect_equal(lfq$response(), old)

  expect_identical(lfq$rename_response("transformedIntensity"), lfq)
  expect_equal(lfq$response(), "transformedIntensity")
  expect_true("transformedIntensity" %in% colnames(lfq$data_long()))
  expect_false(old %in% colnames(lfq$data_long()))
})

test_that("remove_small_intensities re-establishes explicit missing rows", {
  istar_missing <- sim_lfq_data_peptide_config(Nprot = 5, weight_missing = 1, seed = 1)
  lfq_missing <- LFQData$new(istar_missing$data, istar_missing$config)

  filtered <- lfq_missing$get_copy()
  filtered$remove_small_intensities()

  response <- filtered$response()
  expect_gt(sum(is.na(filtered$data_long()[[response]])), 0)
  expect_equal(nrow(filtered$data_long()), nrow(prolfqua::complete_cases(filtered)))
})

test_that("LFQData data_wide conversion", {
  tmp <- lfqdata$data_wide()
  expect_equal(nrow(tmp$data), nrow(tmp$rowdata))
  expect_equal(ncol(tmp$data), nrow(tmp$annotation) + ncol(tmp$rowdata))
  expect_true("data.frame" %in% class(tmp$data))

  tmp_mat <- lfqdata$data_wide(as.matrix = TRUE)
  expect_true("matrix" %in% class(tmp_mat$data))
})

test_that("LFQData data_wide aligns matrix columns and annotation rows", {
  shuffled <- lfqdata$get_copy()
  sample_order <- unique(lfqdata$data_long()$sampleName)
  interleaved <- sample_order[c(9, 2, 10, 1, 11, 3, 5, 12, 4, 6, 7, 8)]
  shuffled$set_data(
    lfqdata$data_long() |>
      dplyr::mutate(.sample_order = match(.data$sampleName, interleaved)) |>
      dplyr::arrange(.data$.sample_order) |>
      dplyr::select(-".sample_order")
  )

  tmp <- shuffled$data_wide()
  sample_cols <- setdiff(colnames(tmp$data), colnames(tmp$rowdata))
  expect_equal(sample_cols, tmp$annotation[[shuffled$sample_name()]])

  tmp_mat <- shuffled$data_wide(as.matrix = TRUE)
  expect_equal(colnames(tmp_mat$data), tmp_mat$annotation[[shuffled$sample_name()]])
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
  expect_warning(
    expect_s3_class(lfqdata$get_Aggregator("medpolish"), "AggregateMedpolish"),
    "did not transform the intensities"
  )
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

test_that("sample offsets adjust the log2 response by sample ID", {
  cfg <- AnalysisConfiguration$new()
  cfg$file_name <- "run"
  cfg$sample_name <- "sample_id"
  cfg$hierarchy <- list(feature = "feature")
  cfg$set_response("log2_intensity")
  cfg$is_response_transformed <- TRUE
  data <- data.frame(
    feature = rep(c("p1", "p2", "p3"), each = 2),
    run = rep(c("a.raw", "b.raw"), 3),
    sample_id = factor(rep(c("A", "B"), 3)),
    log2_intensity = c(9, 5, 13, NA, 7, 10),
    log2_intensity_offset_adjusted = 0
  )
  lfq <- LFQData$new(data, cfg, prefix = "protein_")
  tr <- lfq$get_Transformer()

  result <- withVisible(tr$subtract_sample_offsets(c(B = -2, A = 1)))
  expect_identical(result$value, tr)
  expect_false(result$visible)
  expect_identical(tr$lfq$response(), "log2_intensity_offset_adjusted.1")
  adjusted <- tr$lfq$data_long()
  expect_equal(adjusted[[tr$lfq$response()]], c(8, 7, 12, NA, 6, 12))
  expect_equal(diff(adjusted[[tr$lfq$response()]][c(1, 3, 5)]), diff(data$log2_intensity[c(1, 3, 5)]))
  expect_identical(adjusted[names(data)], data)
  expect_true(tr$lfq$is_transformed())
  expect_identical(tr$lfq$prefix, "protein_")
  expect_identical(lfq$data_long(), data)
  expect_identical(lfq$get_config()$work_intensity, "log2_intensity")
})

test_that("sample offsets require one finite numeric value per sample ID", {
  logged <- lfqdata$get_Transformer()$log2()$lfq
  tr <- logged$get_Transformer()
  sample_ids <- unique(logged$data_long()[[logged$sample_name()]])
  offsets <- setNames(numeric(length(sample_ids)), sample_ids)
  invalid <- list(
    as.character(offsets),
    as.matrix(offsets),
    unname(offsets),
    offsets[-1],
    setNames(offsets, replace(sample_ids, 1, sample_ids[2])),
    setNames(offsets, replace(sample_ids, 1, "")),
    setNames(offsets, replace(sample_ids, 1, NA_character_)),
    c(offsets, unexpected_sample = 0),
    replace(offsets, 1, NA_real_),
    replace(offsets, 1, Inf)
  )
  for (value in invalid) {
    expect_error(tr$subtract_sample_offsets(value), class = "prolfqua_error_bad_argument")
  }
  expect_identical(tr$lfq$data_long(), logged$data_long())
  expect_identical(tr$lfq$response(), logged$response())
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
  expect_equal(nrow(result$data_long()), nrow(result$factors()) * nrow(result$hierarchy()))

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

  # invalid func must error via match.arg, not be silently accepted
  expect_error(AggregateTopN$new(lfqdata$get_copy(), "protein", func = "median"), "should be one of")
})

test_that("Transformer and Summariser deep-clone their input", {
  lfq_copy <- lfqdata$get_copy()
  expect_false(lfq_copy$is_transformed())

  tr <- lfq_copy$get_Transformer()
  tr$log2()
  # decorator's in-place ops must not mutate the caller's LFQData
  expect_false(lfq_copy$is_transformed())
  expect_true(tr$lfq$is_transformed())

  sm <- lfq_copy$get_Summariser()
  expect_false(identical(sm$lfq, lfq_copy))
})

test_that("tidy_to_wide_config takes the row keys", {
  default <- tidy_to_wide_config(lfqdata, as.matrix = TRUE)
  no_label <- tidy_to_wide_config(lfqdata, as.matrix = TRUE, row_keys = lfqdata$hierarchy_keys())
  expect_equal(dim(no_label$data), dim(default$data))
  expect_false(any(grepl(lfqdata$data_long()[[lfqdata$isotope_label()]][1], rownames(no_label$data), fixed = TRUE)))
  expect_equal(colnames(no_label$rowdata), lfqdata$hierarchy_keys())

  # two isotope labels: the hierarchy keys alone do not identify a row
  heavy <- dplyr::mutate(lfqdata$data_long(), !!lfqdata$isotope_label() := "heavy")
  two_labels <- LFQData$new(dplyr::bind_rows(lfqdata$data_long(), heavy), lfqdata$get_config())
  expect_error(
    tidy_to_wide_config(two_labels, row_keys = two_labels$hierarchy_keys()),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("get_Aggregator offers limpa", {
  skip_if_not_installed("limpa")
  lfq_log2 <- lfqdata$get_Transformer()$log2()$lfq
  expect_s3_class(lfq_log2$get_Aggregator("limpa"), "AggregateLimpa")
  expect_error(lfq_log2$get_Aggregator("nope"), class = "prolfqua_error_bad_argument")
})

test_that("LFQDataToSummarizedExperiment returns a SummarizedExperiment", {
  skip_if_not_installed("SummarizedExperiment")
  se <- LFQDataToSummarizedExperiment(lfqdata)
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(ncol(se), nrow(lfqdata$data_wide()$annotation))
})

test_that("lfqdata_from_wide inverts tidy_to_wide_config", {
  istar <- sim_lfq_data_peptide_config(Nprot = 10)
  lfq <- LFQData$new(istar$data, istar$config)
  wide <- tidy_to_wide_config(lfq, as.matrix = TRUE)
  nr <- tidy_to_wide_config(lfq, as.matrix = TRUE, value = lfq$nr_children_col())
  values <- stats::setNames(list(wide$data, nr$data), c(lfq$response(), lfq$nr_children_col()))
  back <- lfqdata_from_wide(values, wide$rowdata, wide$annotation, lfq$get_config())
  expect_s3_class(back, "LFQData")

  keys <- c(lfq$hierarchy_keys(), lfq$sample_name())
  cols <- c(keys, lfq$isotope_label(), lfq$file_name(), lfq$factor_keys(), lfq$response(), lfq$nr_children_col())
  tidy <- function(x) {
    x <- dplyr::select(x$data_long(), dplyr::all_of(cols))
    as.data.frame(dplyr::arrange(x, dplyr::across(dplyr::all_of(keys))))
  }
  expect_equal(tidy(back), tidy(lfq))

  single <- lfqdata_from_wide(wide$data, wide$rowdata, wide$annotation, lfq$get_config())
  expect_equal(single$response(), lfq$response())

  expect_error(
    lfqdata_from_wide(wide$data, wide$rowdata[-1, ], wide$annotation, lfq$get_config()),
    class = "prolfqua_error_bad_argument"
  )
  expect_error(
    lfqdata_from_wide(wide$data, wide$rowdata, wide$annotation[-1, ], lfq$get_config()),
    class = "prolfqua_error_bad_argument"
  )
  expect_error(
    lfqdata_from_wide(wide$data, wide$rowdata[, "protein_Id", drop = FALSE], wide$annotation, lfq$get_config()),
    class = "prolfqua_error_missing_column"
  )
})

test_that("intensity_array can mark a back-transform as untransformed", {
  logged <- lfqdata$get_Transformer()$intensity_array(log)$lfq
  expect_true(logged$is_transformed())
  back <- logged$get_Transformer()$intensity_array(exp, force = TRUE, transformed = FALSE)$lfq
  expect_false(back$is_transformed())
  still <- logged$get_Transformer()$intensity_array(exp, force = TRUE)$lfq
  expect_true(still$is_transformed())
})

test_that("is_aggregatable tells whether there is a level below the modelled one", {
  expect_true(lfqdata$is_aggregatable())
  prot <- sim_lfq_data_protein_config()
  expect_false(LFQData$new(prot$data, prot$config)$is_aggregatable())
  depth2 <- lfqdata$get_copy()
  depth2$set_config_value("hierarchy_depth", 2)
  expect_false(depth2$is_aggregatable())
})
