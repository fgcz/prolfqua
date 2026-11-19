# Provenance: each object records the problems and dropped items that happen in it; one summary per step.

prov_peptide <- function(Nprot = 10, seed = 42) {
  sim <- sim_lfq_data_peptide_config(Nprot = Nprot, with_missing = FALSE, seed = seed)
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
