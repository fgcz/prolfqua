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

inputAnntation <- "C:\\Users\\wewol\\MAS_WEW\\LFQServiceAnalysisTemplate\\inst\\MQ_Ionstar2018_PXD003881/annotationIonstar.xlsx"
assign("lfq_write_format", "xlsx", envir = .GlobalEnv)



# creates default configuration
config <- LFQService::create_MQ_peptide_Configuration()

annotation <- readxl::read_xlsx(inputAnntation)

config$table$factors[["dilution."]] = "sample"
config$table$factors[["run_ID"]] = "run_ID"


config$table$factorLevel <- 1

config$order_Id = "IonStar"
config$project_Id = "p3000"
config$workunit_Id = "IonStar"

# specify model definition

res <- application_set_up_MQ_run(
  outpath = outpath,
  inputMQfile = inputMQfile,
  inputAnnotation = inputAnntation,
  config = config$clone(deep = TRUE)
)


dd <-
  LFQService:::.workflow_MQ_filter_peptides_V3(res$data, res$config)
pepIntensityNormalized <-
  transform_work_intensity(dd$data, dd$config, log2)


subset <-
  pepIntensityNormalized %>% dplyr::filter(grepl("HUMAN", pepIntensityNormalized$protein_Id))
pepIntensityNormalized <-
  scale_with_subset(pepIntensityNormalized, subset, dd$config)

hist(pepIntensityNormalized$log2_peptide.intensity)
hist(pepIntensityNormalized$log2_peptide.intensity_subset_scaled)

pepIntensityNormalized <- pepIntensityNormalized %>%
  dplyr::rename(transformedIntensity = dd$config$table$getWorkIntensity())
dd$config$table$popWorkIntensity()
dd$config$table$setWorkIntensity("transformedIntensity")

res <-
  list(
    pepIntensityNormalized = pepIntensityNormalized,
    config_pepIntensityNormalized = dd$config
  )

protL <- medpolish_protein_quants(res$pepIntensityNormalized, res$config_pepIntensityNormalized)
protL("plot")

resProt <- list(
  protIntensities = protL("unnest")$data,
  config_protIntensities = protL("unnest")$config
)

dataIonstarSubsetNorm_V2 <- list()
dataIonstarSubsetNorm_V2$resultsPep <- res
dataIonstarSubsetNorm_V2$resultsProt <- resProt
usethis::use_data(dataIonstarSubsetNorm_V2, overwrite = TRUE)
