hurdle_components_fixture <- function() {
  istar <- sim_lfq_data_peptide_config(Nprot = 60, weight_missing = 0.6, seed = 7)
  lfq_peptide <- LFQData$new(istar$data, istar$config)$get_Transformer()$log2()$lfq
  lfq_protein <- suppressWarnings(lfq_peptide$get_Aggregator()$aggregate())
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl", "B_vs_Ctrl" = "group_B - group_Ctrl")
  fit <- function(lfqdata, method) {
    suppressWarnings(build_contrast_analysis(lfqdata, "~ group_", contrasts, method = method))
  }
  list(
    lfq_peptide = lfq_peptide,
    contrasts = contrasts,
    lm = fit(lfq_protein, "lm"),
    lm_impute = fit(lfq_protein, "lm_impute"),
    firth_nested = fit(lfq_peptide, "firth_nested"),
    binomial_nested = fit(lfq_peptide, "binomial_nested"),
    binomial_firth_nested = fit(lfq_peptide, "binomial_firth_nested")
  )
}

fx <- hurdle_components_fixture()

# Component get_contrasts() warns about zero sample variances on the simulated data.
new_hurdle <- function(intensity, detection) {
  suppressWarnings(ContrastsHurdleFacade$new(intensity, detection))
}

test_that("hurdle p-values are Fisher's combination of the component p-values", {
  res <- new_hurdle(fx$lm_impute, fx$firth_nested)$get_contrasts()
  both <- res[res$hurdle_components == "both", ]
  expect_true(nrow(both) > 0)
  chisq <- -2 * (log(both$p.value.intensity) + log(both$p.value.detection))
  expect_equal(both$chisq.hurdle, chisq)
  expect_equal(both$df.hurdle, rep(4, nrow(both)))
  expect_equal(both$p.value, stats::pchisq(chisq, 4, lower.tail = FALSE))

  intensity_only <- res[res$hurdle_components == "intensity", ]
  expect_true(nrow(intensity_only) > 0)
  expect_equal(intensity_only$df.hurdle, rep(2, nrow(intensity_only)))
  expect_equal(intensity_only$p.value, intensity_only$p.value.intensity)
})

test_that("proteins absent from lm become detection-only rows", {
  res <- new_hurdle(fx$lm, fx$firth_nested)$get_contrasts()
  lm_rows <- fx$lm$get_contrasts()
  absent <- dplyr::anti_join(res, lm_rows, by = c("protein_Id", "contrast"))
  expect_true(nrow(absent) > 0)
  expect_equal(nrow(res), nrow(fx$firth_nested$get_contrasts()))
  expect_true(all(absent$hurdle_components == "detection"))
  expect_equal(absent$df.hurdle, rep(2, nrow(absent)))
  expect_equal(absent$p.value, absent$p.value.detection)
  expect_true(all(is.na(absent$diff)))
})

test_that("fully detected proteins keep their intensity p-value", {
  completed <- fx$lfq_peptide$get_copy()
  completed$complete_cases()
  fully_detected <- completed$data_long() |>
    dplyr::group_by(.data$protein_Id) |>
    dplyr::summarise(full = !anyNA(.data[[completed$response()]]), .groups = "drop") |>
    dplyr::filter(.data$full)
  expect_true(nrow(fully_detected) > 0)

  for (detection in list(fx$firth_nested, fx$binomial_nested, fx$binomial_firth_nested)) {
    res <- new_hurdle(fx$lm_impute, detection)$get_contrasts()
    full_rows <- res[res$protein_Id %in% fully_detected$protein_Id, ]
    expect_true(all(full_rows$hurdle_components == "intensity"))
    expect_equal(full_rows$p.value, full_rows$p.value.intensity)
    expect_true(all(res$hurdle_components[!res$protein_Id %in% fully_detected$protein_Id] == "both"))
  }
})

test_that("FDR is Benjamini-Hochberg within contrast", {
  res <- new_hurdle(fx$lm, fx$binomial_nested)$get_contrasts()
  expect_setequal(unique(res$contrast), names(fx$contrasts))
  for (ct in names(fx$contrasts)) {
    sub <- res[res$contrast == ct, ]
    expect_equal(sub$FDR, stats::p.adjust(sub$p.value, method = "BH"))
  }
})

test_that("all intensity x detection pairings return the contrast schema", {
  for (intensity in c("lm", "lm_impute")) {
    for (detection in c("firth_nested", "binomial_nested", "binomial_firth_nested")) {
      hurdle <- new_hurdle(fx[[intensity]], fx[[detection]])
      res <- hurdle$get_contrasts()
      expect_true(all(hurdle$column_description()$column_name %in% colnames(res)))
      expect_true(all(
        c(
          "p.value.intensity",
          "FDR.intensity",
          "diff.detection",
          "p.value.detection",
          "FDR.detection",
          "chisq.hurdle",
          "df.hurdle",
          "hurdle_components"
        ) %in%
          colnames(res)
      ))
      expect_setequal(unique(res$modelName), "hurdle")
      expect_false(dplyr::is_grouped_df(res))
      expect_equal(anyDuplicated(res[, c("protein_Id", "contrast")]), 0)
    }
  }
  hurdle <- new_hurdle(fx$lm_impute, fx$firth_nested)
  expect_equal(hurdle$get_config()$subject_id, "protein_Id")
  wide <- hurdle$to_wide()
  expect_equal(nrow(wide), length(unique(hurdle$get_contrasts()$protein_Id)))
  expect_r6_class(hurdle$get_Plotter(), "ContrastsPlotter")
  expect_true(nrow(hurdle$filter_significant(FDR_threshold = 0.5, diff_threshold = 0)) > 0)
})

test_that("hurdle supports a two-factor design", {
  istar <- sim_lfq_data_peptide_config(Nprot = 20, weight_missing = 0.8, seed = 11)
  istar$data <- istar$data |>
    dplyr::mutate(
      Treatment = dplyr::if_else(.data$group_ == "B", "B", "A"),
      Background = dplyr::if_else(grepl("_V[34]$", .data$sampleName), "Z", "X")
    )
  istar$config$factors <- list(Treatment = "Treatment", Background = "Background")
  istar$config$factor_depth <- 2
  lfq_peptide <- LFQData$new(istar$data, istar$config)$get_Transformer()$log2()$lfq
  lfq_protein <- suppressWarnings(lfq_peptide$get_Aggregator()$aggregate())
  contrasts <- c(
    "treatment" = "TreatmentB - TreatmentA",
    "background" = "BackgroundZ - BackgroundX"
  )
  modelstr <- "~ Treatment * Background"
  intensity <- suppressWarnings(build_contrast_analysis(lfq_protein, modelstr, contrasts, method = "lm_impute"))
  detection <- suppressWarnings(build_contrast_analysis(lfq_peptide, modelstr, contrasts, method = "binomial_nested"))
  res <- new_hurdle(intensity, detection)$get_contrasts()
  expect_setequal(unique(res$contrast), names(contrasts))
  for (ct in names(contrasts)) {
    sub <- res[res$contrast == ct, ]
    expect_equal(sub$FDR, stats::p.adjust(sub$p.value, method = "BH"))
  }
})

test_that("mismatched components raise errors", {
  one_contrast <- suppressWarnings(build_contrast_analysis(
    fx$lfq_peptide,
    "~ group_",
    fx$contrasts[1],
    method = "firth_nested"
  ))
  expect_error(new_hurdle(fx$lm, one_contrast), "same contrasts")

  renamed <- dplyr::rename(fx$firth_nested$get_contrasts(), protein = "protein_Id")
  other_subject <- ContrastsTable$new(renamed, subject_id = "protein", model_name = "firth_nested")
  expect_error(new_hurdle(fx$lm, other_subject), "same subject_id")

  expect_error(new_hurdle(fx$lm, fx$lm_impute), "binomial_firth_nested facade")
})

test_that("hurdle_nested builds both components from peptide input", {
  facade <- suppressWarnings(build_contrast_analysis(
    fx$lfq_peptide,
    "~ group_",
    fx$contrasts,
    method = "hurdle_nested"
  ))
  expect_r6_class(facade, "ContrastsHurdleNestedFacade")
  res <- facade$get_contrasts()
  expect_setequal(unique(res$modelName), "hurdle_nested")
  expect_true(all(ContrastsInterface$new()$column_description()$column_name %in% colnames(res)))
  expect_true(all(c("p.value.intensity", "p.value.detection", "hurdle_components") %in% colnames(res)))

  manual <- new_hurdle(fx$lm_impute, fx$binomial_firth_nested)$get_contrasts()
  expect_equal(nrow(res), nrow(manual))
  joined <- dplyr::inner_join(res, manual, by = c("protein_Id", "contrast"), suffix = c("", ".manual"))
  expect_equal(nrow(joined), nrow(manual))
  expect_equal(joined$diff, joined$diff.manual)
  expect_equal(joined$p.value, joined$p.value.manual)
  expect_equal(joined$FDR, joined$FDR.manual)

  expect_s3_class(facade$get_missing(), "data.frame")
  expect_s3_class(facade$to_wide(), "data.frame")
  expect_r6_class(facade$get_Plotter(), "ContrastsPlotter")
  expect_true(nrow(facade$filter_significant(FDR_threshold = 0.5, diff_threshold = 0)) > 0)

  entry <- lookup_facade("hurdle_nested")
  expect_equal(entry$class, "ContrastsHurdleNestedFacade")
  expect_equal(entry$needs, "nested")
  expect_error(
    ContrastsHurdleNestedFacade$new(fx$lm_impute$.lfqdata, "~ group_", fx$contrasts),
    "requires LFQData with additional hierarchy"
  )
})
