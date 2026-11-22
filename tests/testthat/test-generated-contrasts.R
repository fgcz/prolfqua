test_that("main-effect contrasts expose a top-level subtraction", {
  contrasts <- main_effect_contrasts(
    primary_levels = c("WT", "KO"),
    secondary_levels = c("Early", "Late")
  )

  expect_identical(
    unname(contrasts$KO_vs_WT),
    "(G_KO_Early + G_KO_Late)/2 - (G_WT_Early + G_WT_Late)/2"
  )
})

test_that("generated contrasts give the intended differences of group means", {
  # group means: interaction 2, effect of B 3 at x and 5 at y
  means <- list(G_A_x = 1, G_A_y = 2, G_B_x = 4, G_B_y = 7)
  value <- function(contrast) eval(parse(text = contrast), envir = means)
  res <- generate_contrasts(c("A", "B"), c("x", "y"))
  expect_identical(res$ContrastName, c("B_vs_A", "B_vs_A_at_x", "B_vs_A_at_y", "interaction_B_vs_A_at_y_vs_x"))
  expect_equal(vapply(res$Contrast, value, numeric(1), USE.NAMES = FALSE), c(4, 3, 5, 2))
  no_interaction <- generate_contrasts(c("A", "B"), c("x", "y"), interactions = FALSE)
  expect_identical(no_interaction$ContrastName, res$ContrastName[1:3])
})

test_that("contrast generators cover every pair of levels", {
  expect_identical(
    generate_contrasts_for_factor(c("a", "b", "c")),
    list(b_vs_a = "b - a", c_vs_a = "c - a", c_vs_b = "c - b")
  )
  expect_identical(
    level_specific_contrasts(c("A", "B"), c("x", "y")),
    list(B_vs_A_at_x = "G_B_x - G_A_x", B_vs_A_at_y = "G_B_y - G_A_y")
  )
  expect_identical(
    interaction_contrasts(c("A", "B", "C"), c("x", "y")),
    list(
      interaction_B_vs_A_at_y_vs_x = "(G_B_y - G_A_y) - (G_B_x - G_A_x)",
      interaction_C_vs_A_at_y_vs_x = "(G_C_y - G_A_y) - (G_C_x - G_A_x)",
      interaction_C_vs_B_at_y_vs_x = "(G_C_y - G_B_y) - (G_C_x - G_B_x)"
    )
  )
})

test_that("annotation_add_contrasts adds the group column and pads the contrast columns", {
  df <- data.frame(
    Sample = paste0("s", 1:6),
    genotype = rep(c("WT", "KO"), each = 3),
    time = rep(c("t1", "t2", "t1"), 2)
  )
  res <- annotation_add_contrasts(df, "genotype", "time", dataset_id = "d1")
  expect_identical(res$name, "DEA_G__d1.csv")
  expect_identical(res$annot$Group, paste(df$genotype, df$time, sep = "_"))
  expect_identical(res$annot$Sample, df$Sample)
  expect_false(any(c("genotype", "time") %in% colnames(res$annot)))
  # levels sorted: KO before WT, t1 before t2
  expect_identical(
    res$annot$ContrastName,
    c("WT_vs_KO", "WT_vs_KO_at_t1", "WT_vs_KO_at_t2", "interaction_WT_vs_KO_at_t2_vs_t1", NA, NA)
  )
  expect_identical(res$annot$Contrast[2], "G_WT_t1 - G_KO_t1")
})
