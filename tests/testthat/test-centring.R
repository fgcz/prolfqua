# Effect of centring on contrasts (srr RE2.3).

centring_fixture <- function() {
  sim <- sim_lfq_data_protein_config(n_proteins = 10, with_missing = FALSE, seed = 42)
  suppressMessages(LFQData$new(sim$data, sim$config)$get_Transformer()$log2()$lfq)
}

contrast_table <- function(lfq) {
  fit <- build_contrast_analysis(lfq, "~ group_", c(AvsC = "group_A - group_Ctrl", BvsC = "group_B - group_Ctrl"), "lm")
  dplyr::arrange(fit$get_contrasts(), .data$protein_Id, .data$contrast)
}

test_that("preserve_mean changes the average abundance, not the contrasts", {
  withr::local_options(prolfqua.verbose = FALSE)
  lfq <- centring_fixture()
  kept <- contrast_table(lfq$get_Transformer()$robscale(preserve_mean = TRUE)$lfq)
  zero <- contrast_table(lfq$get_Transformer()$robscale(preserve_mean = FALSE)$lfq)
  for (column in c("diff", "std.error", "statistic", "p.value")) {
    expect_equal(kept[[column]], zero[[column]], tolerance = 1e-10, label = column)
  }
  expect_gt(min(abs(kept$avgAbd - zero$avgAbd)), 1)
})

test_that("a shift common to all samples leaves contrasts unchanged; a per-group shift moves them by its difference", {
  withr::local_options(prolfqua.verbose = FALSE)
  lfq <- centring_fixture()
  base <- contrast_table(lfq)
  samples <- lfq$factors()
  common <- stats::setNames(rep(3, nrow(samples)), samples$sampleName)
  shifted <- contrast_table(lfq$get_Transformer()$subtract_sample_offsets(common)$lfq)
  expect_equal(shifted$diff, base$diff, tolerance = 1e-10)
  expect_equal(shifted$p.value, base$p.value, tolerance = 1e-10)
  expect_equal(shifted$avgAbd, base$avgAbd - 3, tolerance = 1e-10)
  by_group <- stats::setNames(ifelse(samples$group_ == "A", 0.5, 0), samples$sampleName)
  moved <- contrast_table(lfq$get_Transformer()$subtract_sample_offsets(by_group)$lfq)
  is_avsc <- moved$contrast == "AvsC"
  expect_equal(moved$diff[is_avsc], base$diff[is_avsc] - 0.5, tolerance = 1e-10)
  expect_equal(moved$diff[!is_avsc], base$diff[!is_avsc], tolerance = 1e-10)
})

# robscale_lowcv ----

lowcv_fixture <- function(with_missing = FALSE) {
  sim <- sim_lfq_data_protein_config(n_proteins = 100, with_missing = with_missing, seed = 42)
  suppressMessages(LFQData$new(sim$data, sim$config)$get_Transformer()$log2()$lfq)
}

response_table <- function(lfq) {
  data <- lfq$data_long()
  keys <- c(lfq$hierarchy_keys(), lfq$sample_name())
  out <- data[keys]
  out$response <- data[[lfq$response()]]
  dplyr::arrange(out, dplyr::across(dplyr::all_of(keys)))
}

test_that("robscale_lowcv removes per-sample shifts and keeps the overall level", {
  withr::local_options(prolfqua.verbose = FALSE)
  lfq <- lowcv_fixture()
  samples <- lfq$factors()$sampleName
  shift <- stats::setNames(seq(0, 2, length.out = length(samples)), samples)
  shifted <- lfq$get_Transformer()$subtract_sample_offsets(-shift)$lfq
  base <- lfq$get_Transformer()$robscale_lowcv()
  moved <- shifted$get_Transformer()$robscale_lowcv()
  expect_identical(moved$lowcv_features, base$lowcv_features)
  difference <- response_table(moved$lfq)$response - response_table(base$lfq)$response
  expect_equal(difference, rep(mean(shift), length(difference)), tolerance = 1e-10)
})

test_that("low-CV features are the requested share of the features observed in every sample", {
  withr::local_options(prolfqua.verbose = FALSE)
  lfq <- lowcv_fixture(with_missing = TRUE)
  data <- lfq$data_long()
  n_samples <- dplyr::n_distinct(data$sampleName)
  observed <- tapply(!is.na(data[[lfq$response()]]), data$protein_Id, sum)
  complete <- names(observed)[observed == n_samples]
  expect_gt(length(complete), 0)
  expect_lt(length(complete), length(observed))
  tr <- lfq$get_Transformer()$robscale_lowcv(fraction = 0.25)
  expect_equal(nrow(tr$lowcv_features), round(0.25 * length(complete)))
  expect_true(all(tr$lowcv_features$protein_Id %in% complete))
  expect_equal(nrow(lfq$get_Transformer()$robscale_lowcv(fraction = 1)$lowcv_features), length(complete))
})

test_that("features that differ between groups are not low-CV features", {
  withr::local_options(prolfqua.verbose = FALSE)
  lfq <- lowcv_fixture()
  data <- lfq$data_long()
  proteins <- unique(data$protein_Id)
  differential <- proteins[seq(1, length(proteins), by = 3)]
  in_a <- data$protein_Id %in% differential & data$group_ == "A"
  data[[lfq$response()]][in_a] <- data[[lfq$response()]][in_a] + 3
  perturbed <- lfq$get_copy()
  perturbed$set_data(data)
  tr <- perturbed$get_Transformer()$robscale_lowcv(fraction = 0.4)
  expect_equal(nrow(tr$lowcv_features), 40)
  expect_false(any(tr$lowcv_features$protein_Id %in% differential))
})

test_that("robscale_lowcv rejects untransformed data, incomplete features and a bad fraction", {
  withr::local_options(prolfqua.verbose = FALSE)
  sim <- sim_lfq_data_protein_config(n_proteins = 20, with_missing = FALSE, seed = 42)
  raw <- suppressMessages(LFQData$new(sim$data, sim$config))
  expect_error(raw$get_Transformer()$robscale_lowcv(), "log-transformed", class = "prolfqua_error_bad_argument")

  lfq <- lowcv_fixture()
  expect_error(lfq$get_Transformer()$robscale_lowcv(fraction = 0), class = "prolfqua_error_bad_argument")
  expect_error(lfq$get_Transformer()$robscale_lowcv(fraction = 1.5), class = "prolfqua_error_bad_argument")

  data <- lfq$data_long()
  samples <- unique(data$sampleName)
  protein_index <- match(data$protein_Id, unique(data$protein_Id))
  missing_sample <- samples[protein_index %% length(samples) + 1]
  data[[lfq$response()]][data$sampleName == missing_sample] <- NA
  incomplete <- lfq$get_copy()
  incomplete$set_data(data)
  expect_error(
    incomplete$get_Transformer()$robscale_lowcv(),
    "observed in every sample",
    class = "prolfqua_error_bad_argument"
  )
})
