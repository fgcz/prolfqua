#' Statistical condition checks
#'
#' @srrstats {G5.8b} Validation tests reject character, complex, logical and factor responses with classed errors.
#' @noRd

validation_fixture <- function() {
  x <- sim_lfq_data_protein_config(Nprot = 5, seed = 42)
  list(
    data = x$data,
    config = x$config,
    lfq = LFQData$new(x$data, x$config),
    contrasts = c(A_vs_Ctrl = "group_A - group_Ctrl")
  )
}

test_that("analysis requires LFQData before inspecting its configuration", {
  expect_error(build_contrast_analysis(list(), "~ group_", c(a = "A - B")), class = "prolfqua_error_bad_argument")
})

test_that("analysis method keys match exactly and are case sensitive", {
  x <- validation_fixture()
  for (method in list("ropecaV", "LM", NA_character_, c("lm", "limma"))) {
    expect_error(
      build_contrast_analysis(x$lfq, "~ group_", x$contrasts, method = method),
      class = "prolfqua_error_bad_argument"
    )
  }
})

test_that("facades require nonempty named character contrasts", {
  x <- validation_fixture()
  for (contrasts in list(
    "group_A - group_Ctrl",
    c(a = 1),
    character(),
    setNames("group_A", ""),
    setNames("group_A", NA_character_)
  )) {
    expect_error(ContrastsLMFacade$new(x$lfq, "~ group_", contrasts), class = "prolfqua_error_bad_argument")
  }
})

test_that("facades require one RHS-only formula string", {
  x <- validation_fixture()
  for (formula in list("abundance ~ group_", "group_", c("~ group_", "~ 1"), NA_character_, "~", ~group_)) {
    expect_error(ContrastsLMFacade$new(x$lfq, formula, x$contrasts), class = "prolfqua_error_bad_argument")
  }
})

test_that("LFQData rejects data of the wrong class", {
  x <- validation_fixture()
  expect_error(LFQData$new(list(), x$config), class = "prolfqua_error_bad_argument")
})

test_that("LFQData rejects configuration of the wrong class", {
  expect_error(LFQData$new(data.frame(), list()), class = "prolfqua_error_invalid_configuration")
})

test_that("LFQData requires the configured response column", {
  x <- validation_fixture()
  x$data[[x$config$get_response()]] <- NULL
  expect_error(LFQData$new(x$data, x$config), class = "prolfqua_error_missing_column")
})

test_that("LFQData requires a numeric response", {
  x <- validation_fixture()
  for (value in list("1", 1 + 1i, TRUE, factor("1"))) {
    x$data[[x$config$get_response()]] <- value
    expect_error(LFQData$new(x$data, x$config), class = "prolfqua_error_bad_argument")
  }
})

test_that("setup_analysis rejects non-data-frame input", {
  x <- validation_fixture()
  expect_error(setup_analysis(list(), x$config), class = "prolfqua_error_bad_argument")
})

test_that("setup_analysis rejects invalid configuration", {
  expect_error(setup_analysis(data.frame(), list()), class = "prolfqua_error_invalid_configuration")
})

test_that("setup_analysis checks source hierarchy factor and response columns", {
  x <- validation_fixture()
  for (role in c("hierarchy", "factors", "response")) {
    cfg <- x$config$clone(deep = TRUE)
    cfg$hierarchy <- list(protein_Id = "protein_Id")
    cfg$factors <- list(group_ = "group_")
    if (role == "response") {
      cfg$set_response("absent")
    } else {
      cfg[[role]][[1]] <- "absent"
    }
    expect_error(setup_analysis(x$data, cfg), class = "prolfqua_error_missing_column")
  }
})

test_that("setup_analysis rejects a character response", {
  x <- validation_fixture()
  x$data[[x$config$get_response()]] <- "1"
  expect_error(setup_analysis(x$data, x$config), class = "prolfqua_error_bad_argument")
})

test_that("matrix strategies enforce the same formula contract as StrategyBase", {
  for (strategy in list(strategy_limma, strategy_limpa)) {
    for (formula in list("abundance group_", c("a ~ b", "c ~ d"), NA_character_, "a ~")) {
      expect_error(strategy(formula), class = "prolfqua_error_bad_argument")
    }
    expect_s3_class(strategy(abundance ~ group_)$formula, "formula")
  }
})

test_that("limma builders reject misspelled weights instead of silently fitting unweighted", {
  x <- validation_fixture()
  strat <- strategy_limma("abundance ~ group_", weights = "misspelled")
  for (builder in list(
    build_model_limma,
    build_model_limma_impute,
    build_model_limma_voom,
    build_model_limma_voom_impute
  )) {
    expect_error(builder(x$lfq, strat), class = "prolfqua_error_missing_column")
  }
})

test_that("per-protein strategies require their weights column", {
  x <- validation_fixture()
  expect_error(
    build_model(x$lfq, strategy_lm("abundance ~ group_", weights = "misspelled")),
    class = "prolfqua_error_missing_column"
  )
})

test_that("limma builders require the configured formula response", {
  x <- validation_fixture()
  for (formula in c("other ~ group_", "~ group_", "log(abundance) ~ group_")) {
    expect_error(build_model_limma(x$lfq, strategy_limma(formula)), class = "prolfqua_error_bad_argument")
  }
})

test_that("imputation requires LFQData", {
  expect_error(build_model_impute(list(), strategy_lm("abundance ~ group_")), class = "prolfqua_error_bad_argument")
})

test_that("imputation requires one finite numeric LOD or NULL", {
  x <- validation_fixture()
  for (lod in list(NA_real_, Inf, NaN, c(1, 2), "1", numeric())) {
    expect_error(
      build_model_impute(x$lfq, strategy_lm("abundance ~ group_"), lod = lod),
      class = "prolfqua_error_bad_argument"
    )
  }
})

test_that("DEqMS requires a data frame of counts", {
  expect_error(
    ContrastsModeratedDEqMS$new(ContrastsInterface$new(), list(), "count"),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("DEqMS requires the counts column", {
  expect_error(
    ContrastsModeratedDEqMS$new(ContrastsInterface$new(), data.frame(), "count"),
    class = "prolfqua_error_missing_column"
  )
})

test_that("DEqMS requires numeric counts", {
  expect_error(
    ContrastsModeratedDEqMS$new(ContrastsInterface$new(), data.frame(count = "1"), "count"),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("DEqMS requires one loess span in (0, 1]", {
  for (span in list(0, -1, 1.1, Inf, NA_real_, c(0.5, 0.75), "0.75")) {
    expect_error(
      ContrastsModeratedDEqMS$new(ContrastsInterface$new(), data.frame(count = 1), "count", loess_span = span),
      class = "prolfqua_error_bad_argument"
    )
  }
})

test_that("Contrasts requires a Model", {
  expect_error(Contrasts$new(list(), c(a = "group_A - group_B")), class = "prolfqua_error_bad_argument")
})

test_that("Contrasts requires named contrasts", {
  x <- validation_fixture()
  model <- build_model(x$lfq, strategy_lm("abundance ~ group_"))
  expect_error(Contrasts$new(model, "group_A - group_B"), class = "prolfqua_error_bad_argument")
})

test_that("moderation requires a contrast object", {
  expect_error(ContrastsModerated$new(list()), class = "prolfqua_error_bad_argument")
})

test_that("simulators require a positive whole number of proteins", {
  for (simulator in list(sim_lfq_data_protein_config, sim_lfq_data_peptide_config, sim_lfq_data_2factor_config)) {
    for (n in list(0, -1, 1.5, NA_real_, Inf, c(1, 2), "5")) {
      expect_error(simulator(Nprot = n), class = "prolfqua_error_bad_argument")
    }
  }
})

test_that("simulators require a missingness weight in [0, 1]", {
  for (simulator in list(sim_lfq_data_protein_config, sim_lfq_data_peptide_config, sim_lfq_data_2factor_config)) {
    for (weight in list(-0.1, 1.1, NA_real_, Inf, c(0, 1), "0.2")) {
      expect_error(simulator(weight_missing = weight), class = "prolfqua_error_bad_argument")
    }
  }
})

test_that("ANOVA test choices match exactly", {
  expect_error(get_anova_df("bad"), class = "prolfqua_error_bad_argument")
})

test_that("Firth histogram choices match exactly", {
  model <- ModelFirth$new(list(model_df = data.frame()))
  expect_error(model$anova_histogram("p"), class = "prolfqua_error_bad_argument")
})

test_that("robscale_subset aborts on mismatched transformation state", {
  x <- validation_fixture()
  ref <- x$lfq$get_copy()$get_Transformer()$log2()$lfq
  expect_error(x$lfq$get_Transformer()$robscale_subset(ref), class = "prolfqua_error_bad_argument")
})
