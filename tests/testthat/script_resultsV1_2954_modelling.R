rm(list=ls())
library(tidyverse)
library(LFQService)
library(lme4)
library(conflicted)


VIS_PROT <- FALSE
Phoenix <- TRUE

results <- LFQService::resultsV12954

results$path <- file.path(tempdir(), results$path)

if(!dir.exists(results$path)){
  dir.create(results$path)
}


length(unique(results$resDataStart$protein_Id))
length(unique(results$filteredPep$protein_Id))


results$config_pepIntensityNormalized$table$factorDepth <- 1
pepConfig<- results$config_pepIntensityNormalized

pepConfig$table$factorKeys()

# Model 1
modelName  <- "f_Condition_r_peptide"

formula_randomPeptide <- make_custom_model_lmer("transformedIntensity  ~ Condition + (1 | peptide_Id)", model_name = modelName)
models_base <- model_analyse(results$pepIntensityNormalized, formula_randomPeptide, modelName)
summary_base <- model_analyse_summarize(models_base$modelProtein,models_base$modelName)




#model_analyse_write(models_base, modelName, results$path)
reslist <- model_analyse_summarize_vis(summary_base)
model_analyse_summarize_vis_write(reslist, path = results$path)



# Model2
modelName  <- "f_Condition_r_peptid_r_patient"
formula_randomPatient <- make_custom_model_lmer("transformedIntensity  ~ Condition +  (1|patient_id) + (1|peptide_Id)", model_name=modelName)

models_interaction <- model_analyse(results$pepIntensityNormalized,
                                    formula_randomPatient,
                                    modelName)


if(TRUE){
  m <- get_complete_model_fit(models_interaction$modelProtein)
  m$linear_model[[1]]
  factor_contrasts <- linfct_factors_contrasts( m$linear_model[[1]])
  factor_levelContrasts <- contrasts_linfct( m,
                                             factor_contrasts,
                                             subject_Id = "protein_Id")
  head(factor_levelContrasts)
  m$linear_model[[1]]
  my_contest(m$linear_model[[1]],factor_contrasts )
}










# usethis::use_data(models_interaction, overwrite = TRUE)

summary_interaction <- model_analyse_summarize(models_interaction$modelProtein,
                                               models_interaction$modelName)

model_analyse_summarize_write(summary_interaction,  results$path)
reslist <- model_analyse_summarize_vis(summary_interaction)
model_analyse_summarize_vis_write(reslist,  path = results$path)

# saveRDS(models_interaction, file="models_interaction.rda")
# saveRDS(res_cond_r_pep_r_pat,file=paste0(modelName,".rda"))
# res_cond_r_pep_r_pat <- readRDS(paste0(modelName,".rda"))

# modelName  <- "f_Condition_r_patient.peptpid"
# formula_randomPatient <- make_custom_model_lmer("transformedIntensity  ~ Condition + (1 | patient_id/peptide_Id)", model_name = modelName)
# names(rres_cond_r_pat.pep)
# saveRDS(rres_cond_r_pat.pep,file=paste0(modelName,".rda"))


lltest <- workflow_likelihood_ratio_test(models_base$modelProtein,
                                         models_base$modelName,
                                         models_interaction$modelProtein,
                                         models_interaction$modelName,
                                         subject_Id = pepConfig$table$hkeysDepth(),
                                         path = results$path)


models <- get_complete_model_fit(models_interaction$modelProtein)
m <- linfct_from_model(models$linear_model[[1]])
m$linfct_interactions
# Group averages for one of the models
models_interaction_Averages <- contrasts_linfct(models,
                                                 m$linfct_interactions,
                                                 subject_Id = pepConfig$table$hkeysDepth() )

contrasts_linfct_write(models_interaction_Averages, pepConfig, models_base$modelName ,  path=results$path )

wfs <- contrasts_linfct_vis(models_interaction_Averages,models_base$modelName )
contrasts_linfct_vis_write(wfs, path=results$path)


all_linfct <- LFQService:::.linfct_all_possible_contrasts(m$linfct_interactions)

models_allContrasts <- contrasts_linfct( models,
                                         all_linfct,
                                         subject_Id = pepConfig$table$hkeysDepth() )

wfs <- contrasts_linfct_vis(models_allContrasts,models_interaction$modelName )
contrasts_linfct_vis_write(wfs, path=results$path)
pivot_model_contrasts_2_Wide(models_allContrasts)

