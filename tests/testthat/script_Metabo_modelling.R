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
modelFunction <- make_custom_model_lm("log2_rawIntensity_robust_scale  ~ Mortality + Intervention  + NRS", model_name =  modelName)

modellingResult_A <- model_analyse(results$dataTransformed,
                                   modelFunction,
                                   modelName,
                                   subject_Id = pepConfig$table$hkeysLevel())
#usethis::use_data(modellingResult_A, overwrite = TRUE)

modelSummary_A <- model_analyse_summarize(modellingResult_A$modelProtein,modelName,subject_Id = pepConfig$table$hkeysLevel())



visualization <- model_analyse_summarize_vis(modelSummary_A,pepConfig$table$hkeysLevel())
model_analyse_summarize_write(modelSummary_A , results$path)
model_analyse_summarize_vis_write(visualization ,path = results$path)



# Compute contrasts between main factors -----
m <- get_complete_model_fit(modellingResult_A$modelProtein)

factor_contrasts <- linfct_factors_contrasts(m$lmer_f_Mortality_Intervention_NRS[[1]])
factor_levelContrasts <- contrasts_linfct( m,
                                           modellingResult_A$modelName,
                                           factor_contrasts,
                                           subject_Id = pepConfig$table$hkeysLevel(),
                                           contrastfun = LFQService::my_contrast_V2)


wfs <- contrasts_linfct_vis(factor_levelContrasts,
                            modellingResult_A$modelName,
                            subject_Id = "Compound")

contrasts_linfct_vis_write(wfs, path=results$path)
contrasts_linfct_write(factor_levelContrasts,
                       modellingResult_A$modelName ,
                       path=results$path,
                       subject_Id = "Compound" )

# Compute subgroup averages ----
linfct <- linfct_from_model(m$lmer_f_Mortality_Intervention_NRS[[1]])
models_interaction_Averages <- contrasts_linfct( m,
                                                 modelSummary_A$modelName,
                                                 linfct$linfct_factors,
                                                 subject_Id = pepConfig$table$hkeysLevel(),
                                                contrastfun = LFQService::my_contrast_V2,
                                                )

contrasts_linfct_write(models_interaction_Averages,
                       modellingResult_A$modelName ,
                       prefix = "GroupAverages",
                       path=results$path,
                       subject_Id = "Compound" )

wfs <- contrasts_linfct_vis(models_interaction_Averages,
                            modellingResult_A$modelName ,
                            prefix = "GroupAverages",
                            subject_Id = "Compound")
contrasts_linfct_vis_write(wfs, path=results$path)

contrastres_fun <- workflow_contrasts_linfct(m,
                                             modellingResult_A$modelName,
                                             linfct$linfct_factors ,
                                             subject_Id = pepConfig$table$hkeysLevel(),
                                             prefix = "GroupAverages",
                                             contrastfun = LFQService::my_contrast_V2)

contrasts <- contrastres_fun()
res <- contrastres_fun(path=results$path)



# second model For likelihood ratio test ----

results$config_pepIntensityNormalized$table$factorLevel <- 2

modelName  <- "f_Mortality_Intervention"
modelFunction <- make_custom_model_lm("log2_rawIntensity_robust_scale   ~ Mortality + Intervention ", modelName)

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
dim(m)
dim(modellingResult_B$modelProtein)

anova(m$lmer_f_Mortality_Intervention_NRS[[1]],modellingResult_B$modelProtein$lmer_f_Mortality_Intervention[[1]] )
lltest <- workflow_likelihood_ratio_test(m,
                                         modelSummary_A$modelName,
                                         modellingResult_B$modelProtein,
                                         modellingResult_B$modelName,
                                         subject_Id = "Compound",
                                         path = NULL)

lltest <- workflow_likelihood_ratio_test(m,
                                         modelSummary_A$modelName,
                                         modellingResult_B$modelProtein,
                                         modellingResult_B$modelName,
                                         subject_Id = "Compound",
                                         path = results$path)



