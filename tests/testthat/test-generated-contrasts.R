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
