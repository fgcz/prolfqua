context("test-tidyconfig_functions")

test_that("check config", {
  config <- create_config_Skyline()
  config$table$factors[["Time"]] = "Sampling.Time.Point"
  expect_equal(config$table$factor_keys(),"Time")
  expect_equal(config$table$hierarchy_keys(),c("protein_Id","peptide_Id","precursor_Id","fragment_Id"))
})


test_that("my_contrast_V2 works",{
  data_lm_models_to_test <- prolfqua_data('data_lm_models_to_test')
  linfct_lm <- linfct_from_model(data_lm_models_to_test$lm_complete)
  linfct_interaction <- linfct_lm$linfct_interactions

  for (model in data_lm_models_to_test) {
    print(my_contrast_V2(model, linfct_interaction))
  }
  linfct_contrasts <- linfct_factors_contrasts(data_lm_models_to_test$lm_complete)

})
