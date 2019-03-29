rm(list=ls())
library(tidyverse)
library(LFQService)
library(multcomp)
library(lme4)

VIS_PROT <- FALSE

Phoenix <- TRUE


results <- LFQService::resultsV12954
dir.create(results$path)
length(unique(results$resDataStart$protein_Id))
length(unique(results$filteredPep$protein_Id))


results$config_pepIntensityNormalized$table$factorLevel <- 1
pepConfig<- results$config_pepIntensityNormalized
pepConfig$table$factorKeys()

pepConfig$table$factorKeys()
modelName  <- "f_Condition_r_peptide"

formula_randomPeptide <- make_custom_model("transformedIntensity  ~ Condition + (1 | peptide_Id)")
res_cond_r_pep <- workflow_no_interaction_modelling(results,formula_randomPeptide, modelName)
#  saveRDS(res_cond_r_pep,file=paste0(modelName,".rda"))
#res_cond_r_pep <- readRDS(paste0(modelName,".rda"))

modelName  <- "f_Condition_r_peptid_r_patient"
formula_randomPatient <- make_custom_model("transformedIntensity  ~ Condition + (1 | peptide_Id) + (1|patient_id)")
res_cond_r_pep_r_pat <- workflow_no_interaction_modelling(results,formula_randomPatient, modelName)


#  saveRDS(res_cond_r_pep_r_pat,file=paste0(modelName,".rda"))

#res_cond_r_pep_r_pat <- readRDS(paste0(modelName,".rda"))

#modelName  <- "f_Condition_r_patient.peptpid"
#formula_randomPatient <- make_custom_model("transformedIntensity  ~ Condition + (1 | patient_id/peptide_Id)")
#rres_cond_r_pat.pep <- workflow_no_interaction_modelling(results,formula_randomPatient, modelName)
#names(rres_cond_r_pat.pep)
#saveRDS(rres_cond_r_pat.pep,file=paste0(modelName,".rda"))


lltest <- workflow_likelihood_ratio_test(res_cond_r_pep$modelProteinF,
                                         res_cond_r_pep$modelName,
                                         res_cond_r_pep_r_pat$modelProteinF,
                                         res_cond_r_pep_r_pat$modelName)

pdf(file = file.path(results$path, "LLratioTest.pdf"))
hist(lltest$likelihood_ratio_test.pValue, breaks=20, main="Testing the significance of the random effect")
dev.off()

res_cond_r_pep$TwoFactorModelFactor2 <- inner_join(res_cond_r_pep$TwoFactorModelFactor2,lltest , by = "protein_Id")
res_cond_r_pep_r_pat$TwoFactorModelFactor2 <- inner_join(res_cond_r_pep_r_pat$TwoFactorModelFactor2,lltest , by = "protein_Id")

#readr::write_csv(res_cond_r_pep$TwoFactorModelFactor2, path = file.path(results$path, paste0("Contrasts_Auto_SignificanceValues_",  res_cond_r_pep$modelName, ".csv")))
#readr::write_csv(res_cond_r_pep_r_pat$TwoFactorModelFactor2, path = file.path(results$path, paste0("Contrasts_Auto_SignificanceValues_",  res_cond_r_pep_r_pat$modelName, ".csv")))


res_cond_r_pep_Pivot <- pivot_model_contrasts_2_Wide(res_cond_r_pep$TwoFactorModelFactor2)
#write_csv(res_cond_r_pep_Pivot, path=file.path(results$path, paste0("Contrasts_SignificanceValues_", res_cond_r_pep$modelName, "_PIVOT.csv")))

res_cond_r_pep_r_pat_Pivot <- pivot_model_contrasts_2_Wide(res_cond_r_pep_r_pat$TwoFactorModelFactor2)
#write_csv(res_cond_r_pep_r_pat_Pivot,
#          path=file.path(results$path, paste0("Contrasts_SignificanceValues_", res_cond_r_pep_r_pat$modelName, "_PIVOT.csv")))

tmp <- inner_join(res_cond_r_pep$TwoFactorModelFactor2,  res_cond_r_pep_r_pat$TwoFactorModelFactor2, by=c("protein_Id", "factor", 'isSingular', 'peptide_Id_n', 'lhs'))
ggplot(tmp, aes(x = estimate.x, y = estimate.y, )) + geom_point() + facet_wrap(~lhs)
ggplot(tmp, aes(x = p.value.x, y = p.value.y, )) + geom_point() + facet_wrap(~lhs)



# add group averages -----

m <- res_cond_r_pep$modelProteinF$lmer_f_Condition_r_peptide[[1]]
linfct <- lmer4_linfct_from_model(m)
res_cond_r_pep_grA <- workflow_group_averages(res_cond_r_pep$modelProteinF,res_cond_r_pep$modelName,
                                              results$path, linfct$linfct_interactions, lltest )
res_cond_r_pep_r_pat_grA <- workflow_group_averages(res_cond_r_pep_r_pat$modelProteinF,res_cond_r_pep_r_pat$modelName,
                                                    results$path, linfct$linfct_interactions,lltest )

#write_csv(res_cond_r_pep_grA, path=file.path(results$path, paste0("GroupAverages_", res_cond_r_pep$modelName, ".csv")))
#write_csv(res_cond_r_pep_r_pat_grA, path=file.path(results$path, paste0("GroupAverages_", res_cond_r_pep_r_pat$modelName, ".csv")))


# Doing protein plots with prediction -----
res_cond_r_pep$modelProteinF <- res_cond_r_pep$modelProteinF %>%
  mutate(nRandom_Interaction_plot = map2(!!sym(paste0("lmer_", res_cond_r_pep$modelName)), protein_Id, plot_model_and_data))


library(lme4)
#plot_model_and_data_TWO(res_cond_r_pep_r_pat$modelProteinF$lmer_f_Condition_r_peptid_r_patient[[1]],"A",legend.position = "bottom", firstlast = TRUE)
#plot_model_and_data_TWO(res_cond_r_pep_r_pat$modelProteinF$lmer_f_Condition_r_peptid_r_patient[[1]],"A",legend.position = "bottom", firstlast = FALSE)

res_cond_r_pep_r_pat$modelProteinF <- res_cond_r_pep_r_pat$modelProteinF %>%
  mutate(nRandom_Interaction_plot = map2(!!sym(paste0("lmer_", res_cond_r_pep_r_pat$modelName)), protein_Id, plot_model_and_data_TWO, firstlast = TRUE))

tmp <- inner_join(dplyr::select(res_cond_r_pep$modelProteinF, protein_Id, nRandom_Interaction_plot),
                  dplyr::select(res_cond_r_pep_r_pat$modelProteinF, protein_Id, nRandom_Interaction_plot), by="protein_Id" )

tmp <- inner_join(tmp, lltest)

tmp <- tmp %>% dplyr::filter(likelihood_ratio_test.pValue < 0.000001)
pdf(file.path(results$path,"Significant_Random_Patient.pdf"), width =8, height =8)
for(i in 1:nrow(tmp)){
  print(gridExtra::grid.arrange(tmp$nRandom_Interaction_plot.x[[i]], tmp$nRandom_Interaction_plot.y[[i]]))
}
dev.off()
