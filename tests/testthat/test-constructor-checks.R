# Constructors reject malformed arguments with a typed error.

test_that("strategies need a model formula", {
  expect_error(strategy_lm("abundance group_"), class = "prolfqua_error_bad_argument")
  expect_error(strategy_lm(c("a ~ b", "c ~ d")), class = "prolfqua_error_bad_argument")
  expect_s3_class(strategy_lm("abundance ~ group_")$formula, "formula")
  expect_s3_class(strategy_lm(abundance ~ group_)$formula, "formula")
})

test_that("Model needs a model table with the subject columns and a strategy", {
  strat <- strategy_lm("abundance ~ group_")
  expect_error(
    Model$new(model_df = list(), model_strategy = strat, model_name = "m"),
    class = "prolfqua_error_bad_argument"
  )
  expect_error(
    Model$new(model_df = data.frame(x = 1), model_strategy = strat, model_name = "m"),
    class = "prolfqua_error_missing_column"
  )
  expect_error(
    Model$new(model_df = data.frame(protein_Id = "P1"), model_strategy = list(), model_name = "m"),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("ModelFirth needs models1, models2 and hkey", {
  expect_error(ModelFirth$new(list(models1 = NULL)), class = "prolfqua_error_bad_argument")
  expect_s3_class(ModelFirth$new(list(models1 = NULL, models2 = NULL, hkey = NULL)), "ModelFirth")
})

test_that("ContrastsFirth needs a ModelFirth and named contrasts", {
  mod <- ModelFirth$new(list(models1 = NULL, models2 = NULL, hkey = NULL))
  expect_error(ContrastsFirth$new(list(), c(a = "group_A - group_B")), class = "prolfqua_error_bad_argument")
  expect_error(ContrastsFirth$new(mod, "group_A - group_B"), class = "prolfqua_error_bad_argument")
})
