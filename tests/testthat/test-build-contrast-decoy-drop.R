# build_contrast_analysis drops decoys before the fit only when opted in.

make_data_with_decoys <- function(n_decoy = 5, Nprot = 30) {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = Nprot)
  data <- istar$data
  config <- istar$config
  top <- config$hierarchy_keys()[1]
  prots <- unique(data[[top]])
  data[[top]] <- ifelse(
    data[[top]] %in% prots[seq_len(n_decoy)],
    paste0("REV_", data[[top]]),
    data[[top]]
  )
  list(data = data, config = config, top = top)
}

test_that("pattern_decoys unset -> decoys remain in contrasts (unchanged behaviour)", {
  s <- make_data_with_decoys()
  lfq <- prolfqua::LFQData$new(s$data, s$config)
  lfq$rename_response("transformedIntensity")
  fa <- prolfqua::build_contrast_analysis(
    lfq,
    "~ group_",
    c("A_vs_Ctrl" = "group_A - group_Ctrl"),
    method = "lm"
  )
  res <- fa$get_contrasts()
  expect_true(any(prolfqua::is_decoy(unique(res[[s$top]]))))
})

test_that("pattern_decoys set -> no decoy in contrasts (targets-only fit)", {
  s <- make_data_with_decoys()
  config <- s$config$clone(deep = TRUE)
  config$pattern_decoys <- "" # opt in with defaults only
  lfq <- prolfqua::LFQData$new(s$data, config)
  lfq$rename_response("transformedIntensity")
  fa <- prolfqua::build_contrast_analysis(
    lfq,
    "~ group_",
    c("A_vs_Ctrl" = "group_A - group_Ctrl"),
    method = "lm"
  )
  res <- fa$get_contrasts()
  expect_false(any(prolfqua::is_decoy(unique(res[[s$top]]))))
  # every modelled protein is a target (a subset of the non-decoy ids)
  all_ids <- unique(s$data[[s$top]])
  target_ids <- all_ids[!prolfqua::is_decoy(all_ids)]
  expect_true(all(unique(res[[s$top]]) %in% target_ids))
})

test_that("targets-only fit also applies to a moderated (deqms) facade", {
  s <- make_data_with_decoys()
  config <- s$config$clone(deep = TRUE)
  config$pattern_decoys <- ""
  lfq <- prolfqua::LFQData$new(s$data, config)
  lfq$rename_response("transformedIntensity")
  fa <- prolfqua::build_contrast_analysis(
    lfq,
    "~ group_",
    c("A_vs_Ctrl" = "group_A - group_Ctrl"),
    method = "deqms"
  )
  res <- fa$get_contrasts()
  expect_false(any(prolfqua::is_decoy(unique(res[[s$top]]))))
})
