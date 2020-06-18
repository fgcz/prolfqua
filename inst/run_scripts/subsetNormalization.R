# Run mixed models on benchmark dataset.

rm(list = ls())

library(conflicted)
library(LFQService)
library(tidyverse)
library(dplyr)

conflict_prefer("filter", "dplyr")



inputMQfile <-
  "C:/Users/wewol/MAS_WEW/LFQServiceAnalysisTemplate/inst/benchmarkData/MQ_Ionstar2018_PXD003881.zip"
outpath <- "results_modelling_all"
#outpath <- "results_modelling_WHO_noSex"

inputAnnotation <- "C:\\Users\\wewol\\MAS_WEW\\LFQServiceAnalysisTemplate\\inst\\MQ_Ionstar2018_PXD003881/annotationIonstar.xlsx"

ps <- ProjectStructure$new(outpath = ".",
                     project_Id = 3000,
                     order_Id = "IonStar" ,
                     workunit_Id = "IonStar",
                     inputData = inputMQfile,
                     inputAnnotation = inputAnnotation)

mqdata <- tidyMQ_Peptides_Config(ps$inputData)

annotation <- readxl::read_xlsx(ps$inputAnnotation)

mqdata$config$table$factors[["dilution."]] = "sample"
mqdata$config$table$factors[["run_ID"]] = "run_ID"
mqdata$config$table$factorDepth <- 1




# specify model definition

res <- application_add_annotation(
  mqdata$data,
  annotation
)

mqdata$data <- setup_analysis(res, mqdata$config)

dd <-
  LFQService:::filter_proteins_by_peptide_count(mqdata$data, mqdata$config )

pepIntensityNormalized <-
  transform_work_intensity(dd$data, mqdata$config, log2)

subset <-
  pepIntensityNormalized %>% dplyr::filter(grepl("HUMAN", pepIntensityNormalized$protein_Id))

pepIntensityNormalized <-
  scale_with_subset(pepIntensityNormalized, subset, mqdata$config)

hist(pepIntensityNormalized$log2_peptide.intensity)
hist(pepIntensityNormalized$log2_peptide.intensity_subset_scaled)

pepIntensityNormalized <- pepIntensityNormalized %>%
  dplyr::rename(transformedIntensity = mqdata$config$table$getWorkIntensity())
mqdata$config$table$popWorkIntensity()
mqdata$config$table$setWorkIntensity("transformedIntensity")

res <-
  list(
    pepIntensityNormalized = pepIntensityNormalized,
    config_pepIntensityNormalized = mqdata$config
  )

protL <- medpolish_protein_quants(res$pepIntensityNormalized, res$config_pepIntensityNormalized)
#protL("plot")

resProt <- list(
  protIntensities = protL("unnest")$data,
  config_protIntensities = protL("unnest")$config
)

dataIonstarSubsetNorm_V2 <- list()
dataIonstarSubsetNorm_V2$resultsPep <- res
dataIonstarSubsetNorm_V2$resultsProt <- resProt
usethis::use_data(dataIonstarSubsetNorm_V2, overwrite = TRUE)
