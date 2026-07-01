# Shared decoy / contaminant identifier detection.

test_that("is_decoy flags default decoy prefixes and not normal ids", {
  ids <- c(
    "REV_sp|P1|X",
    "rev_P2",
    "DECOY_P3",
    "decoy_P4",
    "XXX_P5",
    "reverse_P6",
    "##P7",
    "sp|P8|X",
    "tr|P9|X",
    "normalProtein"
  )
  expect_equal(
    is_decoy(ids),
    c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
  )
})

test_that("is_decoy: NULL / empty / 'a^' pattern use defaults only (never match-all)", {
  ids <- c("REV_P1", "P2")
  for (pat in list(NULL, "", "a^")) {
    expect_equal(is_decoy(ids, pattern = pat), c(TRUE, FALSE))
  }
})

test_that("is_decoy unions a configured pattern with the defaults", {
  ids <- c("shuffled_P1", "REV_P2", "P3")
  expect_equal(is_decoy(ids, pattern = "^shuffled_"), c(TRUE, TRUE, FALSE))
})

test_that("is_contaminant flags default contaminant prefixes", {
  ids <- c("zz|Cont00001|X", "CON__ALBU", "CON_KRT", "Cont_trypsin", "contam_1", "sp|P1|X")
  expect_equal(is_contaminant(ids), c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE))
})

test_that("is_contaminant unions a configured pattern with the defaults", {
  ids <- c("KERATIN_1", "zz|C1|X", "sp|P1|X")
  expect_equal(is_contaminant(ids, pattern = "^KERATIN_"), c(TRUE, TRUE, FALSE))
})

test_that("effective_*_pattern exposes the applied regex", {
  expect_true(grepl("REV_", effective_decoy_pattern(), fixed = TRUE))
  expect_true(grepl("shuffled", effective_decoy_pattern("^shuffled_"), fixed = TRUE))
  expect_identical(effective_decoy_pattern(""), effective_decoy_pattern(NULL))
  expect_true(grepl("zz", effective_contaminant_pattern(), fixed = TRUE))
})
