# Provenance: each object records the problems and dropped items that happen in it; one summary per step.

prov_peptide <- function(n_proteins = 10, seed = 42) {
  sim <- sim_lfq_data_peptide_config(n_proteins = n_proteins, with_missing = FALSE, seed = seed)
  LFQData$new(sim$data, sim$config)
}

test_that("Provenance records entries, skips empty ones and reports each entry once", {
  log <- Provenance$new()
  log$record("step_a", "rows removed", c("p1", "p2", "p2"), "decoy")
  log$record("step_a", "nothing", character(0), "none")
  log$record("step_b", "fit failed", "p3", "singular design", failure = TRUE)
  entries <- log$entries()
  expect_equal(entries$step, c("step_a", "step_b"))
  expect_equal(entries$n, c(2L, 1L))
  expect_equal(entries$ids[[1]], c("p1", "p2"))
  expect_equal(entries$failure, c(FALSE, TRUE))
  expect_message(
    expect_warning(log$report(), "step_b: 1 fit failed .*singular design", class = "prolfqua_fit_failed"),
    "step_a: 2 rows removed \\(e.g. p1, p2\\): decoy",
    class = "prolfqua_message"
  )
  expect_silent(log$report())
})

test_that("a copied log is independent of the original", {
  log <- Provenance$new()
  log$record("a", "x", "p1", "r")
  copy <- log$clone(deep = TRUE)
  copy$record("b", "y", "p2", "r")
  expect_equal(nrow(log$entries()), 1L)
  expect_equal(nrow(copy$entries()), 2L)
})

test_that("remove_small_intensities records values set to NA and features removed", {
  lfq <- prov_peptide()
  data <- lfq$data_long()
  response <- lfq$response()
  features <- unique(data$peptide_Id)
  data[[response]][data$peptide_Id == features[1]] <- 1
  data[[response]][data$peptide_Id == features[2]][1] <- 1
  lfq <- LFQData$new(data, lfq$get_config())
  expect_message(lfq$remove_small_intensities(threshold = 4), "values set to NA", class = "prolfqua_message")
  entries <- lfq$provenance()
  set_na <- entries[entries$event == "values set to NA", ]
  removed <- entries[entries$event == "features removed", ]
  n_samples <- length(unique(data$sampleName))
  expect_equal(set_na$n, n_samples + 1L)
  expect_equal(removed$n, 1L)
  expect_match(removed$ids[[1]], features[1], fixed = TRUE)
  expect_equal(set_na$response, response)
})

test_that("remove_decoys and remove_contaminants continue the log and leave the original unchanged", {
  lfq <- prov_peptide()
  data <- lfq$data_long()
  proteins <- unique(data$protein_Id)
  data$protein_Id[data$protein_Id == proteins[1]] <- paste0("REV_", proteins[1])
  data$protein_Id[data$protein_Id == proteins[2]] <- paste0("zz|CON_", proteins[2])
  lfq <- LFQData$new(data, lfq$get_config())
  withr::local_options(prolfqua.verbose = FALSE)
  targets <- lfq$remove_decoys()$remove_contaminants()
  entries <- targets$provenance()
  expect_equal(entries$step, c("remove_decoys", "remove_contaminants"))
  expect_equal(entries$ids, list(paste0("REV_", proteins[1]), paste0("zz|CON_", proteins[2])))
  expect_equal(nrow(lfq$provenance()), 0L)
})

test_that("filter_proteins_by_peptide_count records the proteins removed", {
  lfq <- prov_peptide()
  withr::local_options(prolfqua.verbose = FALSE)
  counts <- table(unique(lfq$data_long()[c("protein_Id", "peptide_Id")])$protein_Id)
  lfq$filter_proteins_by_peptide_count(threshold = 2)
  entries <- lfq$provenance()
  expect_equal(entries$step, "filter_proteins_by_peptide_count")
  expect_setequal(entries$ids[[1]], names(counts)[counts < 2])
})

test_that("a log transform records zeros set to NA; the log continues through transforms and subsets", {
  lfq <- prov_peptide()
  data <- lfq$data_long()
  data[[lfq$response()]][1:3] <- 0
  lfq <- LFQData$new(data, lfq$get_config())
  expect_message(
    transformed <- lfq$get_Transformer()$log2()$robscale()$lfq,
    "log transform: 3 values set to NA",
    class = "prolfqua_message"
  )
  entries <- transformed$provenance()
  expect_equal(entries$n, 3L)
  expect_equal(entries$response, lfq$response())
  expect_equal(nrow(lfq$provenance()), 0L)
  subset <- transformed$get_subset(transformed$data_long()[1, ])
  expect_equal(subset$provenance(), entries)
  expect_equal(transformed$get_copy()$provenance(), entries)
})

test_that("prolfqua.verbose = FALSE silences the provenance message, not the record", {
  lfq <- prov_peptide()
  data <- lfq$data_long()
  data[[lfq$response()]][1] <- 0
  lfq <- LFQData$new(data, lfq$get_config())
  withr::local_options(prolfqua.verbose = FALSE)
  expect_silent(transformed <- lfq$get_Transformer()$log2()$lfq)
  expect_equal(nrow(transformed$provenance()), 1L)
})

test_that("rlm aggregation gives one warning naming the non-converged proteins and records them", {
  sim <- sim_lfq_data_peptide_config(n_proteins = 20, seed = 42)
  withr::local_options(prolfqua.verbose = FALSE)
  lfq <- LFQData$new(sim$data, sim$config)$get_Transformer()$log2()$lfq
  n_warnings <- 0L
  protein <- withCallingHandlers(
    lfq$get_Aggregator("rlm")$aggregate(),
    prolfqua_fit_failed = function(w) {
      n_warnings <<- n_warnings + 1L
      expect_match(conditionMessage(w), "'rlm' failed to converge")
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(n_warnings, 1L)
  entries <- protein$provenance()
  expect_equal(entries$step, "aggregate (rlm)")
  expect_true(entries$failure)
  expect_true(all(entries$ids[[1]] %in% protein$data_long()$protein_Id))
  expect_equal(nrow(lfq$provenance()), 0L)
})

# protein-level data where the first protein is observed in group Ctrl only
prov_one_group <- function() {
  sim <- sim_lfq_data_protein_config(n_proteins = 10, with_missing = FALSE, seed = 42)
  data <- sim$data
  target <- unique(data$protein_Id)[1]
  data$abundance[data$protein_Id == target & data$group_ != "Ctrl"] <- NA
  lfq <- LFQData$new(data, sim$config)$get_Transformer()$log2()$lfq
  list(lfq = lfq, target = target)
}

prov_fit_warnings <- function(expr) {
  messages <- character(0)
  result <- withCallingHandlers(
    expr,
    prolfqua_fit_failed = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(result = result, messages = messages)
}

test_that("a failed fit gives one warning per step and a provenance entry with its reason", {
  withr::local_options(prolfqua.verbose = FALSE)
  fixture <- prov_one_group()
  res <- prov_fit_warnings(build_model(fixture$lfq, strategy_lm(paste(fixture$lfq$response(), "~ group_"))))
  expect_length(res$messages, 1L)
  expect_match(res$messages, "contrasts can be applied only to factors with 2 or more levels")
  entries <- res$result$provenance()
  failed <- entries[entries$event == "fits failed", ]
  expect_equal(failed$ids, list(fixture$target))
  expect_true(failed$failure)
  dropped <- entries[entries$event == "missing responses dropped", ]
  expect_equal(nrow(dropped), 0L)
})

test_that("missing responses dropped per fit are counted (G2.14b)", {
  withr::local_options(prolfqua.verbose = FALSE)
  fixture <- prov_one_group()
  data <- fixture$lfq$data_long()
  second <- unique(data$protein_Id)[2]
  data[[fixture$lfq$response()]][data$protein_Id == second][1:2] <- NA
  lfq <- LFQData$new(data, fixture$lfq$get_config())
  n_target <- sum(is.na(data[[fixture$lfq$response()]][data$protein_Id == fixture$target]))
  # lm: the one-group protein fails to fit, so none of its responses count as dropped; limma fits it with NA
  # coefficients and drops its missing responses
  expected <- list(lm = list(n = 2L, ids = second), limma = list(n = 2L + n_target, ids = c(fixture$target, second)))
  for (method in names(expected)) {
    facade <- suppressWarnings(build_contrast_analysis(lfq, "~ group_", c(AvsC = "group_A - group_Ctrl"), method))
    dropped <- facade$model$provenance()
    dropped <- dropped[dropped$event == "missing responses dropped", ]
    expect_equal(dropped$n, expected[[method]]$n, label = method)
    expect_setequal(dropped$ids[[1]], expected[[method]]$ids)
  }
})

test_that("limma records non-estimable coefficients as singular fits; lm_impute records the LOD refit", {
  withr::local_options(prolfqua.verbose = FALSE)
  fixture <- prov_one_group()
  contrasts <- c(AvsC = "group_A - group_Ctrl")
  limma <- prov_fit_warnings(build_contrast_analysis(fixture$lfq, "~ group_", contrasts, "limma"))
  expect_length(limma$messages, 0L)
  singular <- limma$result$model$provenance()
  expect_equal(singular$ids[singular$event == "singular fits"], list(fixture$target))
  for (method in c("lm_impute", "limma_impute")) {
    imputed <- prov_fit_warnings(build_contrast_analysis(fixture$lfq, "~ group_", contrasts, method))
    expect_length(imputed$messages, 0L)
    entries <- imputed$result$model$provenance()
    expect_equal(entries$ids[entries$event == "refit at the limit of detection"], list(fixture$target), label = method)
  }
})

test_that("rlm non-convergence and lmer failures each give one warning for the fitting step", {
  withr::local_options(prolfqua.verbose = FALSE)
  sim <- sim_lfq_data_peptide_config(n_proteins = 20, seed = 42)
  peptide <- LFQData$new(sim$data, sim$config)$get_Transformer()$log2()$lfq
  protein <- suppressWarnings(peptide$get_Aggregator()$aggregate())
  rlm <- prov_fit_warnings(build_contrast_analysis(protein, "~ group_", c(AvsC = "group_A - group_Ctrl"), "rlm"))
  expect_length(rlm$messages, 1L)
  expect_match(rlm$messages, "fits with a warning.*'rlm' failed to converge")
  lmer <- suppressMessages(prov_fit_warnings(build_contrast_analysis(
    peptide,
    "~ group_ + (1 | peptide_Id) + (1 | sampleName)",
    c(AvsC = "group_A - group_Ctrl"),
    "lmer_nested"
  )))
  expect_length(lmer$messages, 1L)
  expect_match(lmer$messages, "fits failed")
})

test_that("an analysis lists a protein observed in one group only as absent, with the reason", {
  withr::local_options(prolfqua.verbose = FALSE)
  fixture <- prov_one_group()
  contrasts <- c(AvsC = "group_A - group_Ctrl")
  lm <- prov_fit_warnings(build_contrast_analysis(fixture$lfq, "~ group_", contrasts, "lm"))$result
  entries <- lm$provenance()
  absent <- entries[entries$source == "contrasts", ]
  expect_equal(absent$ids, list(fixture$target))
  expect_match(absent$reason, "^fit failed: .*contrasts can be applied only to factors")
  expect_true(fixture$target %in% entries$ids[[which(entries$source == "model" & entries$event == "fits failed")]])
  limma <- build_contrast_analysis(fixture$lfq, "~ group_", contrasts, "limma")$provenance()
  absent <- limma[limma$source == "contrasts", ]
  expect_equal(absent$ids, list(fixture$target))
  expect_equal(absent$reason, "contrast not estimable")
})

test_that("an analysis reads the input data record, e.g. zeros set to NA by log2", {
  sim <- sim_lfq_data_protein_config(n_proteins = 10, with_missing = FALSE, seed = 42)
  sim$data$abundance[1:2] <- 0
  withr::local_options(prolfqua.verbose = FALSE)
  lfq <- LFQData$new(sim$data, sim$config)$get_Transformer()$log2()$lfq
  facade <- build_contrast_analysis(lfq, "~ group_", c(AvsC = "group_A - group_Ctrl"), "lm")
  entries <- facade$provenance()
  zeros <- entries[entries$source == "input data", ]
  expect_equal(zeros$step, "log transform")
  expect_equal(zeros$n, 2L)
  expect_setequal(zeros$ids[[1]], sim$data$protein_Id[1:2])
})

test_that("hurdle_nested reads the input, the aggregated data and both component models", {
  sim <- sim_lfq_data_peptide_config(n_proteins = 20, weight_missing = 0.4, seed = 42)
  withr::local_options(prolfqua.verbose = FALSE)
  lfq <- LFQData$new(sim$data, sim$config)$get_Transformer()$log2()$lfq
  facade <- suppressWarnings(build_contrast_analysis(
    lfq,
    "~ group_",
    c(AvsC = "group_A - group_Ctrl"),
    "hurdle_nested"
  ))
  sources <- unique(facade$.recorded_provenance()$source)
  expect_true(all(sources %in% c("input data", "aggregated data", "model (intensity)", "model (detection)")))
  expect_true("model (intensity)" %in% sources)
})

test_that("print() and summary() show a one-line provenance summary", {
  withr::local_options(prolfqua.verbose = FALSE)
  fixture <- prov_one_group()
  facade <- prov_fit_warnings(
    build_contrast_analysis(fixture$lfq, "~ group_", c(AvsC = "group_A - group_Ctrl"), "lm")
  )$result
  expect_output(print(facade), "Provenance: 1 failure and [0-9]+ other entries; see \\$provenance\\(\\)")
  expect_output(print(summary(facade)), "Provenance: 1 failure and [0-9]+ other entries")
})
