rm(list=ls())
library(tidyverse)
library(LFQService)
library(multcomp)
VIS_PROT <- FALSE

results <- LFQService::results_MetaboData
results$config_pepIntensityNormalized$table$factorLevel <- 3

pepConfig<- results$config_dataTransformed
pepConfig$table$factorLevel <- 3
pepConfig$table$factorKeys()

# first model ----

modelName  <- "f_Mortality_Intervention_NRS"
modelFunction <- make_custom_model_lm("log2_rawIntensity_robust_scale  ~ Mortality + Intervention  + NRS")

modellingResult_A <- model_analyse(results$dataTransformed,
                                   modelFunction,
                                   modelName,
                                   subject_Id = pepConfig$table$hkeysLevel())



modelSummary_A <- model_analyse_summarize(modellingResult_A$modelProtein,modelName,subject_Id = pepConfig$table$hkeysLevel())
usethis::use_data(modelSummary_A, overwrite = TRUE)


visualization <- model_analyse_summarize_vis(modelSummary_A,pepConfig$table$hkeysLevel())
model_analyse_summarize_write(modelSummary_A , results$path)
model_analyse_summarize_vis_write(visualization ,path = results$path)



# Compute contrasts between main factors -----
m <- get_complete_model_fit(modellingResult_A$modelProtein)

factor_contrasts <- linfct_factors_contrasts(m)
factor_levelContrasts <- workflow_contrasts_linfct( modelSummary_A$modelProteinF,
                                                    modelSummary_A$modelName,
                                                    factor_contrasts,
                                                    subject_Id = pepConfig$table$hkeysLevel() )


wfs <- workflow_contrasts_linfct_vis(factor_levelContrasts,
                                     modellingResult_A$modelName,
                                     subject_Id = "Compound")

workflow_contrasts_linfct_vis_write(wfs, path=results$path)
workflow_contrasts_linfct_write(factor_levelContrasts,
                                modellingResult_A$modelName ,
                                path=results$path,
                                subject_Id = "Compound" )

# Compute subgroup averages ----
linfct <- linfct_from_model(m)
models_interaction_Averages <- workflow_contrasts_linfct( modelSummary_A$modelProteinF,
                                                          modelSummary_A$modelName,
                                                          linfct$linfct_factors,
                                                          subject_Id = pepConfig$table$hkeysLevel() )

workflow_contrasts_linfct_write(models_interaction_Averages,
                                modellingResult_A$modelName ,
                                prefix = "GroupAverages",
                                path=results$path,
                                subject_Id = "Compound" )

wfs <- workflow_contrasts_linfct_vis(models_interaction_Averages,
                                     modellingResult_A$modelName ,
                                     prefix = "GroupAverages",
                                     subject_Id = "Compound")
workflow_contrasts_linfct_vis_write(wfs, path=results$path)

contrastres_fun <- workflow_contrasts_linfct_ALL(modelSummary_A, linfct$linfct_factors , subject_Id = pepConfig$table$hkeysLevel(), prefix = "GroupAverages" )
contrasts <- contrastres_fun()
res <- contrastres_fun(path=results$path)



# second model For likelihood ratio test ----

results$config_pepIntensityNormalized$table$factorLevel <- 2

modelName  <- "f_Mortality_Intervention"
modelFunction <- make_custom_model_lm("log2_rawIntensity_robust_scale   ~ Mortality + Intervention ")

pepIntensity <- results$dataTransformed
pepConfig <- results$config_dataTransformed

modellingResult_B <- model_analyse(results$dataTransformed,
                                   modelFunction,
                                   modelName, subject_Id = "Compound")

res <-workflow_model_analyse(results$dataTransformed,
                             modelFunction,
                             modelName,
                             subject_Id = pepConfig$table$hkeysLevel())

## rund model comparison ----
lltest <- workflow_likelihood_ratio_test(modelSummary_A$modelProteinF,
                                         modelSummary_A$modelName,
                                         res()$summaryResult$modelProteinF,
                                         res()$summaryResult$modelName,
                                         subject_Id = "Compound",
                                         path = results$path)




