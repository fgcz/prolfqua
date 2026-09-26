make_complete_lfqdata <- function(Nprot = 10, paired = FALSE) {
  istar <- sim_lfq_data_protein_config(Nprot = Nprot, with_missing = FALSE, paired = paired, seed = 42)
  LFQData$new(istar$data, istar$config)
}

set_missing <- function(lfqdata, protein, keep) {
  data <- lfqdata$data_long()
  drop <- data$protein_Id == protein & !keep(data)
  data$abundance[drop] <- NA
  LFQData$new(data, lfqdata$get_config())
}

imputed_cells <- function(before, after) {
  dplyr::inner_join(
    before$data_long(),
    after$data_long(),
    by = c("protein_Id", "sampleName"),
    suffix = c(".before", ".after")
  )
}

test_that("a ~ group_ fit fills a missing cell with the group mean of the observed values", {
  lfq <- make_complete_lfqdata()
  protein <- lfq$data_long()$protein_Id[1]
  missing_sample <- lfq$data_long() |>
    dplyr::filter(protein_Id == protein, group_ == "A") |>
    dplyr::pull(sampleName) |>
    dplyr::first()
  lfq_na <- set_missing(lfq, protein, function(d) d$sampleName != missing_sample)
  mod <- build_model(lfq_na, strategy_lm("abundance ~ group_"))

  res <- impute_from_model(mod, lfq_na)

  observed_a <- lfq_na$data_long() |>
    dplyr::filter(protein_Id == protein, group_ == "A", !is.na(abundance)) |>
    dplyr::pull(abundance)
  filled <- res$lfqdata$data_long() |>
    dplyr::filter(protein_Id == protein, sampleName == missing_sample) |>
    dplyr::pull(abundance)
  expect_equal(filled, mean(observed_a))
  row <- res$summary[res$summary$protein_Id == protein, ]
  expect_equal(row$route, "fitted")
  expect_equal(row$n_imputed, 1L)
  expect_equal(row$n_observed, nrow(lfq_na$factors()) - 1L)
})

test_that("an additive paired design fills a missing cell with group plus pair effect", {
  lfq <- make_complete_lfqdata(paired = TRUE)
  protein <- lfq$data_long()$protein_Id[1]
  group_effect <- c(A = 0, B = 1, Ctrl = 2)
  pair_effect <- c(P1 = 0, P2 = 0.5, P3 = -0.5, P4 = 1)
  data <- lfq$data_long()
  exact <- data$protein_Id == protein
  data$abundance[exact] <- 20 + group_effect[data$group_[exact]] + pair_effect[data$subject_[exact]]
  lfq <- LFQData$new(data, lfq$get_config())
  lfq_na <- set_missing(lfq, protein, function(d) !(d$group_ == "B" & d$subject_ == "P4"))
  mod <- build_model(lfq_na, strategy_lm("abundance ~ group_ + subject_"))

  res <- impute_from_model(mod, lfq_na)

  filled <- res$lfqdata$data_long() |>
    dplyr::filter(protein_Id == protein, group_ == "B", subject_ == "P4") |>
    dplyr::pull(abundance)
  expect_equal(filled, 20 + 1 + 1)
  expect_equal(res$summary$route[res$summary$protein_Id == protein], "fitted")
})

test_that("a subject missing a whole group gets finite predictions through lod_refit", {
  lfq <- make_complete_lfqdata()
  protein <- lfq$data_long()$protein_Id[1]
  lfq_na <- set_missing(lfq, protein, function(d) d$group_ != "A")
  lod <- min(lfq_na$data_long()$abundance, na.rm = TRUE)
  mod <- build_model_impute(lfq_na, strategy_lm("abundance ~ group_"), lod = lod)
  lod_fit <- mod$model_df$linear_model[[which(mod$model_df$protein_Id == protein)]]
  expect_s3_class(lod_fit, "imputed_model")
  expect_s3_class(lod_fit, "lm")

  res <- impute_from_model(mod, lfq_na)

  filled <- res$lfqdata$data_long() |>
    dplyr::filter(protein_Id == protein, group_ == "A") |>
    dplyr::pull(abundance)
  expect_length(filled, 4)
  expect_true(all(is.finite(filled)))
  row <- res$summary[res$summary$protein_Id == protein, ]
  expect_equal(row$route, "lod_refit")
  expect_equal(row$n_imputed, 4L)
})

test_that("nearly complete data still has a LOD, so a subject missing a whole group is refitted", {
  lfq <- make_complete_lfqdata()
  protein <- lfq$data_long()$protein_Id[1]
  lfq_na <- set_missing(lfq, protein, function(d) d$group_ != "A")
  mod <- build_model_impute(lfq_na, strategy_lm("abundance ~ group_"))

  res <- impute_from_model(mod, lfq_na)

  expect_equal(res$summary$route[res$summary$protein_Id == protein], "lod_refit")
  expect_false(anyNA(res$lfqdata$data_long()$abundance))
})

test_that("a subject missing a whole pair is refitted at the LOD and its pair cells are predicted", {
  lfq <- make_complete_lfqdata(paired = TRUE)
  protein <- lfq$data_long()$protein_Id[1]
  lfq_na <- set_missing(lfq, protein, function(d) d$subject_ != "P1")
  lod <- min(lfq_na$data_long()$abundance, na.rm = TRUE)
  mod <- build_model_impute(lfq_na, strategy_lm("abundance ~ group_ + subject_"), lod = lod)

  res <- impute_from_model(mod, lfq_na)

  filled <- res$lfqdata$data_long() |>
    dplyr::filter(protein_Id == protein, subject_ == "P1") |>
    dplyr::pull(abundance)
  expect_length(filled, 3)
  expect_true(all(is.finite(filled)))
  expect_equal(res$summary$route[res$summary$protein_Id == protein], "lod_refit")
})

test_that("the weighted lm_impute facade model fills rows absent from the long table", {
  istar <- sim_lfq_data_protein_config(Nprot = 30, weight_missing = 0.5, seed = 7)
  data <- istar$data[!is.na(istar$data$abundance), ]
  lfq <- LFQData$new(data, istar$config)
  facade <- ContrastsLMImputeFacade$new(lfq, "~ group_", c(AvsCtrl = "group_A - group_Ctrl"))
  expect_equal(facade$model$model_strategy$weights, "nrPeptides")

  res <- impute_from_model(facade$model, lfq)

  expect_equal(nrow(res$lfqdata$data_long()), nrow(lfq$factors()) * dplyr::n_distinct(data$protein_Id))
  filled <- res$summary$route %in% c("fitted", "lod_refit")
  expect_true(any(filled))
  expect_equal(res$summary$n_observed + res$summary$n_imputed, rep(nrow(lfq$factors()), nrow(res$summary)))
})

test_that("observed values are unchanged and every subject x sample row is present", {
  istar <- sim_lfq_data_protein_config(Nprot = 30, weight_missing = 0.5, seed = 42)
  lfq <- LFQData$new(istar$data, istar$config)
  data <- lfq$data_long()
  lfq <- LFQData$new(data[!is.na(data$abundance) | seq_len(nrow(data)) %% 2 == 0, ], lfq$get_config())
  mod <- build_model_impute(lfq, strategy_lm("abundance ~ group_"))

  res <- impute_from_model(mod, lfq)

  n_samples <- nrow(lfq$factors())
  n_proteins <- dplyr::n_distinct(lfq$data_long()$protein_Id)
  expect_equal(nrow(res$lfqdata$data_long()), n_samples * n_proteins)
  cells <- imputed_cells(lfq, res$lfqdata)
  observed <- !is.na(cells$abundance.before)
  expect_equal(cells$abundance.after[observed], cells$abundance.before[observed])
  expect_equal(sum(res$summary$n_observed), sum(!is.na(lfq$data_long()$abundance)))
  expect_false(identical(res$lfqdata$get_config(), lfq$get_config()))
})

test_that("a complete subject gets route complete and no imputed cells", {
  lfq <- make_complete_lfqdata()
  mod <- build_model(lfq, strategy_lm("abundance ~ group_"))

  res <- impute_from_model(mod, lfq)

  expect_true(all(res$summary$route == "complete"))
  expect_true(all(res$summary$n_imputed == 0L))
  cells <- imputed_cells(lfq, res$lfqdata)
  expect_equal(nrow(cells), nrow(lfq$data_long()))
  expect_equal(cells$abundance.after, cells$abundance.before)
})

test_that("build_model fits that cannot predict the missing cells get route none", {
  lfq <- make_complete_lfqdata(paired = TRUE)
  proteins <- unique(lfq$data_long()$protein_Id)
  lfq_na <- set_missing(lfq, proteins[1], function(d) d$group_ != "A")
  lfq_na <- set_missing(lfq_na, proteins[2], function(d) d$subject_ != "P1")
  mod <- build_model(lfq_na, strategy_lm("abundance ~ group_ + subject_"))

  res <- impute_from_model(mod, lfq_na)

  none <- res$summary |> dplyr::filter(protein_Id %in% proteins[1:2])
  expect_equal(none$route, c("none", "none"))
  expect_equal(none$n_imputed, c(0L, 0L))
  expect_equal(none$n_observed, c(8L, 9L))
  still_na <- res$lfqdata$data_long() |> dplyr::filter(protein_Id %in% proteins[1:2], is.na(abundance))
  expect_equal(nrow(still_na), 4L + 3L)
})

test_that("a model without lm fits is an error", {
  lfq <- make_complete_lfqdata(paired = TRUE)
  mod <- build_model(lfq, strategy_lmer("abundance ~ group_ + (1 | subject_)"))

  expect_error(impute_from_model(mod, lfq), "needs lm fits")
})
