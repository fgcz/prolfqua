# Quant-side decoy / contaminant handling on LFQData.

make_lfq_with_prefixes <- function(decoy = character(), cont = character()) {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 20)
  data <- istar$data
  config <- istar$config
  top <- config$hierarchy_keys()[1]
  prots <- unique(data[[top]])
  relabel <- function(ids, targets, prefix) {
    ifelse(ids %in% targets, paste0(prefix, ids), ids)
  }
  data[[top]] <- relabel(data[[top]], prots[decoy], "REV_")
  data[[top]] <- relabel(data[[top]], prots[cont], "zz|CON|")
  list(lfq = prolfqua::LFQData$new(data, config), top = top, n_prot = length(prots))
}

test_that("remove_decoys drops decoy proteins using built-in defaults (no pattern)", {
  s <- make_lfq_with_prefixes(decoy = 1:4)
  expect_null(s$lfq$get_config()$pattern_decoys)
  clean <- s$lfq$remove_decoys()
  ids <- unique(clean$data_long()[[s$top]])
  expect_false(any(prolfqua::is_decoy(ids)))
  expect_equal(length(ids), s$n_prot - 4)
  # self is not mutated
  expect_true(any(prolfqua::is_decoy(unique(s$lfq$data_long()[[s$top]]))))
})

test_that("decoy_proportion is nr decoy subjects / nr subjects", {
  s <- make_lfq_with_prefixes(decoy = 1:4)
  expect_equal(s$lfq$decoy_proportion(), 4 / 20)
  s0 <- make_lfq_with_prefixes()
  expect_equal(s0$lfq$decoy_proportion(), 0)
})

test_that("a configured pattern_decoys unions with the defaults", {
  s <- make_lfq_with_prefixes(decoy = 1:3)
  data <- s$lfq$data_long()
  top <- s$top
  prots <- unique(data[[top]])
  # relabel one more protein with a non-default, custom decoy prefix
  target <- prots[!prolfqua::is_decoy(prots)][1]
  data[[top]] <- ifelse(data[[top]] == target, paste0("shuffled_", data[[top]]), data[[top]])
  config <- s$lfq$get_config()
  config$pattern_decoys <- "^shuffled_"
  lfq <- prolfqua::LFQData$new(data, config)
  # 3 REV_ (defaults) + 1 shuffled_ (configured) = 4 decoys
  expect_equal(lfq$decoy_proportion(), 4 / 20)
  clean <- lfq$remove_decoys()
  expect_equal(length(unique(clean$data_long()[[top]])), 20 - 4)
})

test_that("contaminants only act when pattern_contaminants is configured", {
  s <- make_lfq_with_prefixes(cont = 1:5)
  # no pattern -> not identified, unchanged, proportion 0
  expect_null(s$lfq$get_config()$pattern_contaminants)
  expect_equal(s$lfq$contaminant_proportion(), 0)
  same <- s$lfq$remove_contaminants()
  expect_equal(length(unique(same$data_long()[[s$top]])), s$n_prot)

  # opt in with a pattern (union with defaults catches the zz|CON| prefix)
  config <- s$lfq$get_config()
  config$pattern_contaminants <- "^KRT"
  lfq <- prolfqua::LFQData$new(s$lfq$data_long(), config)
  expect_equal(lfq$contaminant_proportion(), 5 / 20)
  clean <- lfq$remove_contaminants()
  expect_equal(length(unique(clean$data_long()[[s$top]])), s$n_prot - 5)
})

test_that("decoy removal propagates from protein to peptide rows", {
  istar <- prolfqua::sim_lfq_data_peptide_config()
  data <- istar$data
  config <- istar$config
  top <- config$hierarchy_keys()[1]
  prots <- unique(data[[top]])
  decoy <- prots[1:2]
  data[[top]] <- ifelse(data[[top]] %in% decoy, paste0("REV_", data[[top]]), data[[top]])
  lfq <- prolfqua::LFQData$new(data, config)
  clean <- lfq$remove_decoys()
  # every peptide belonging to a decoy protein is gone
  expect_false(any(prolfqua::is_decoy(unique(clean$data_long()[[top]]))))
  expect_lt(nrow(clean$data_long()), nrow(lfq$data_long()))
})
