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
modellingResult_A <- workflow_model_analyse(results$dataTransformed,
                                            results$config_dataTransformed,
                                            modelFunction,
                                            modelName,
                                            isSingular = isSingular_lm)


visualization <- workflow_model_analyse_vis(modellingResult_A,pepConfig ,  modelName)
workflow_model_analyse_vis_write(visualization,modelName,path = results$path)


## Compute contrasts between main factors
m <- modellingResult_A$modelProteinF$lmer_f_Mortality_Intervention_NRS[[1]]
#LFQService:::.contrast_tukey_multcomp(m,"Intervention")
factor_contrasts <- linfct_factors_contrasts(m)
factor_levelContrasts <- workflow_contrasts_linfct( modellingResult_A$modelProteinF,
                                                    modellingResult_A$modelName,
                                                    factor_contrasts,
                                                    subject_Id = pepConfig$table$hkeysLevel() )

workflow_contrasts_linfct_write(factor_levelContrasts,
                                modellingResult_A$modelName ,
                                path=results$path,
                                subject_Id = "Compound" )

wfs <- workflow_contrasts_linfct_vis(factor_levelContrasts,
                                     modellingResult_A$modelName,
                                     subject_Id = "Compound")
wfs$histogram_coeff_p.values
wfs$VolcanoPlot
workflow_contrasts_linfct_vis_write(wfs, path=results$path)

# Compute subgroup averages
linfct <- linfct_from_model(m)
models_interaction_Averages <- workflow_contrasts_linfct( modellingResult_A$modelProteinF,
                                                          modellingResult_A$modelName,
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




# second model For likelihood ratio test ----

results$config_pepIntensityNormalized$table$factorLevel <- 2

modelName  <- "f_Mortality_Intervention"
modelFunction <- make_custom_model_lm("log2_rawIntensity_robust_scale   ~ Mortality + Intervention ")

pepIntensity <- results$dataTransformed
pepConfig <- results$config_dataTransformed

modellingResult_B <- workflow_model_analyse(results$dataTransformed,
                                            results$config_dataTransformed,
                                            modelFunction,
                                            modelName,
                                            isSingular = isSingular_lm)


## rund model comparison ----

lltest <- workflow_likelihood_ratio_test(modellingResult_A$modelProteinF,
                                         modellingResult_A$modelName,
                                         modellingResult_B$modelProteinF,
                                         modellingResult_B$modelName,
                                         subject_Id = "Compound",
                                         path = results$path)




