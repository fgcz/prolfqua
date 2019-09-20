rm(list=ls())
library(tidyverse)
library(LFQService)
library(conflicted)

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
factor_contrasts <- linfct_factors_contrasts(m$linear_model[[1]])


# Compute subgroup averages ----
linfct <- linfct_from_model(m$linear_model[[1]], as_list = F)
linfct_A <- rbind(linfct, factor_contrasts)

models_interaction_Averages <- contrasts_linfct( m,
                                                 linfct,
                                                 subject_Id = pepConfig$table$hkeysLevel(),
                                                 contrastfun = modelFunction$contrast_fun
)


head(models_interaction_Averages)

modelProteinF <- modellingResult_A$modelProtein %>%
  dplyr::filter(exists_lmer == TRUE)
res_contrasts <- workflow_contrasts_linfct(modelProteinF,
                                           linfct_A,
                                           pepConfig,
                                           prefix =  "Contrasts",
                                           contrastfun = modelFunction$contrast_fun)
tmp <- res_contrasts()$contrast_result

xx <- tmp %>% dplyr::select(pepConfig$table$hierarchyKeys(), "lhs", "estimate")
xx <- xx %>% pivot_wider(names_from = "lhs", values_from = "estimate")


