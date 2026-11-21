## Tests for the injectable progress reporter (.make_progress + wiring)

# Aggregated protein-level data (the "aggregated" facades, e.g. lm, need this).
make_progress_protein_data <- function(n_proteins = 30) {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = n_proteins, seed = 42)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  lfqdata
}

CONTRASTS_P <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
MODELSTR_P <- "~ group_"

# Record (i, total, label) tuples, keeping the max i and the total per label.
make_recorder <- function() {
  rec <- new.env(parent = emptyenv())
  rec$seen <- list()
  rec$callback <- function(i, total, label) {
    prev_i <- if (is.null(rec$seen[[label]])) 0 else rec$seen[[label]][["max_i"]]
    rec$seen[[label]] <- list(max_i = max(prev_i, i), total = total)
    invisible(NULL)
  }
  rec
}

# ---- .make_progress unit behaviour ----

test_that(".make_progress yields a safe no-op reporter for empty loops", {
  for (total in list(0, -1, NA_integer_, integer(0))) {
    pb <- prolfqua:::.make_progress(total)
    expect_silent(pb$tick())
    expect_null(pb$tick(2))
  }
})

test_that(".make_progress (NULL reporter) returns the legacy progress_bar", {
  skip_if_not_installed("progress")
  pb <- prolfqua:::.make_progress(3, reporter = NULL)
  expect_s3_class(pb, "progress_bar")
})

test_that("function reporter gets monotonic ticks reaching total with the label", {
  rec <- make_recorder()
  is_seen <- integer(0)
  cb <- function(i, total, label) is_seen[[length(is_seen) + 1]] <<- i
  pb <- prolfqua:::.make_progress(5, label = "unit", reporter = cb)
  for (k in 1:5) {
    pb$tick()
  }
  expect_true(all(diff(is_seen) > 0)) # monotonic
  expect_equal(is_seen[length(is_seen)], 5) # final i == total

  rec2 <- make_recorder()
  pb2 <- prolfqua:::.make_progress(4, label = "unit2", reporter = rec2$callback)
  for (k in 1:4) {
    pb2$tick()
  }
  expect_equal(rec2$seen[["unit2"]][["max_i"]], 4)
  expect_equal(rec2$seen[["unit2"]][["total"]], 4)
})

test_that("unknown reporter is rejected", {
  expect_error(prolfqua:::.make_progress(3, reporter = "nope"), "prolfqua.progress")
})

test_that("'message' reporter emits a heartbeat in a non-interactive session", {
  pb <- prolfqua:::.make_progress(3, label = "msg", reporter = "message")
  expect_message(pb$tick(), "msg: 1/3")
})

# ---- Wiring: the option reaches the real pipeline ----

test_that("prolfqua.progress callback fires through build_contrast_analysis (lm)", {
  lfqdata <- make_progress_protein_data()
  rec <- make_recorder()
  withr::local_options(prolfqua.progress = rec$callback)

  fa <- prolfqua::build_contrast_analysis(lfqdata, MODELSTR_P, CONTRASTS_P, method = "lm")
  res <- fa$get_contrasts()

  expect_true(is.data.frame(res))
  expect_true(length(rec$seen) > 0) # progress fired at least once
  # Every pass that emitted ran to completion (final tick always emits).
  for (lab in names(rec$seen)) {
    expect_equal(rec$seen[[lab]][["max_i"]], rec$seen[[lab]][["total"]], info = lab)
  }
})

test_that("firth pass labels appear in progress output (firth)", {
  skip_if_not_installed("logistf")
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 30, with_missing = TRUE, weight_missing = 0.5, seed = 42)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  rec <- make_recorder()
  withr::local_options(prolfqua.progress = rec$callback)

  fa <- prolfqua::build_contrast_analysis(lfqdata, MODELSTR_P, CONTRASTS_P, method = "firth")
  res <- fa$get_contrasts()

  expect_true(is.data.frame(res))
  expect_true(any(grepl("^firth", names(rec$seen)))) # distinct firth labels seen
})

# ---- Regression: the reporter never alters results ----

test_that("setting prolfqua.progress does not change get_contrasts() results", {
  res_default <- prolfqua::build_contrast_analysis(
    make_progress_protein_data(),
    MODELSTR_P,
    CONTRASTS_P,
    method = "lm"
  )$get_contrasts()

  withr::local_options(prolfqua.progress = function(i, total, label) invisible(NULL))
  res_opt <- prolfqua::build_contrast_analysis(
    make_progress_protein_data(),
    MODELSTR_P,
    CONTRASTS_P,
    method = "lm"
  )$get_contrasts()

  expect_equal(res_default, res_opt)
})

test_that(".make_progress falls back to a no-op reporter without the progress package", {
  local_mocked_bindings(.has_package = function(pkg) FALSE)
  pb <- .make_progress(3)
  expect_silent(pb$tick())
})
