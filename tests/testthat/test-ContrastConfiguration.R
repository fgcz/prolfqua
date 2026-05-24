test_that("ContrastConfiguration constructor uses LM-flavoured defaults", {
  cfg <- ContrastConfiguration$new(subject_id = "protein_Id")
  expect_equal(cfg$subject_id, "protein_Id")
  expect_equal(cfg$contrast_col, "contrast")
  expect_equal(cfg$effect_col, "diff")
  expect_equal(cfg$score_col, "statistic")
  expect_equal(cfg$pvalue_col, "p.value")
  expect_equal(cfg$fdr_col, "FDR")
  expect_equal(cfg$avg_abundance_col, "avgAbd")
  expect_true(cfg$supports_dea_qc)
  expect_false(cfg$needs_saint_annotation)
})

test_that("ContrastConfiguration fields are the canonical accessors", {
  cfg <- ContrastConfiguration$new(subject_id = c("p", "q"))
  expect_equal(cfg$subject_id, c("p", "q"))
  expect_equal(cfg$contrast_col, "contrast")
  expect_equal(cfg$effect_col, "diff")
  expect_equal(cfg$score_col, "statistic")
  expect_equal(cfg$pvalue_col, "p.value")
  expect_equal(cfg$fdr_col, "FDR")
  expect_equal(cfg$model_name_col, "modelName")
  expect_equal(cfg$avg_abundance_col, "avgAbd")
})

test_that("has_pvalue reflects pvalue_col presence", {
  cfg <- ContrastConfiguration$new(subject_id = "p")
  expect_true(cfg$has_pvalue())

  cfg_na <- ContrastConfiguration$new(subject_id = "p", pvalue_col = NA_character_)
  expect_false(cfg_na$has_pvalue())

  cfg_empty <- ContrastConfiguration$new(subject_id = "p", pvalue_col = "")
  expect_false(cfg_empty$has_pvalue())
})

test_that("SAINT-flavoured config carries non-default column names and flags", {
  saint <- ContrastConfiguration$new(
    subject_id = "protein_Id",
    contrast_col = "Bait",
    effect_col = "log2_EFCs",
    score_col = "SaintScore",
    pvalue_col = NA_character_,
    fdr_col = "BFDR",
    supports_dea_qc = FALSE,
    needs_saint_annotation = TRUE
  )
  expect_equal(saint$contrast_col, "Bait")
  expect_equal(saint$effect_col, "log2_EFCs")
  expect_equal(saint$score_col, "SaintScore")
  expect_equal(saint$fdr_col, "BFDR")
  expect_false(saint$has_pvalue())
  expect_false(saint$supports_dea_qc)
  expect_true(saint$needs_saint_annotation)
})

test_that("deep clone is independent", {
  cfg <- ContrastConfiguration$new(subject_id = "p")
  cfg2 <- cfg$clone(deep = TRUE)
  cfg2$contrast_col <- "Bait"
  expect_equal(cfg$contrast_col, "contrast")
  expect_equal(cfg2$contrast_col, "Bait")
})
