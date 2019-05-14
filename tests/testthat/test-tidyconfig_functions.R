context("test-tidyconfig_functions")

test_that("check config", {
  config <- createSkylineConfiguration()
  config$table$factors[["Time"]] = "Sampling.Time.Point"
  expect_equal(config$table$factorKeys(),"Time")
  expect_equal(config$table$hierarchyKeys(),c("protein_Id","peptide_Id","precursor_Id","fragment_Id"))
})

test_that("run script_MetaboExample_create.R",{
  source("script_MetaboExample_create.R", echo=FALSE)
})

test_that("run script_Metabo_modelling.R",{
  source("script_Metabo_modelling.R", echo=FALSE)
})

test_that("run test_resultsV1_2954_modelling.R",{
  source("script_resultsV1_2954_modelling.R", echo=FALSE)
})

test_that("run script_resultsV1_2954_create.R",{
  source("script_resultsV1_2954_create.R", echo=FALSE)
})

test_that("my_contrast_V2 works",{
  lm_models_to_test <- LFQService::lm_models_to_test
  linfct_lm <- linfct_from_model(lm_models_to_test$lm_complete)
  linfct_interaction <- linfct_lm$linfct_interactions

  for(model in lm_models_to_test){
    print(my_contrast_V2(model, linfct_interaction))
  }
  linfct_contrasts <- linfct_factors_contrasts(lm_models_to_test$lm_complete)

  for(model in lm_models_to_test){
    print(my_contrast_V2(model, linfct_contrasts))
  }

})
