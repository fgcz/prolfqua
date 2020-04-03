# Run mixed models on benchmark dataset.

rm(list = ls())


library(LFQService)
library(tidyverse)
library(dplyr)

inputMQfile <- "../samples/test_MQ_IonStar2018_PXD003881.zip"

inputMQfile <- "C:\\Users\\wolski\\MAS_WEW/LFQServiceAnalysisTemplate//inst/benchmarkData/MQ_Ionstar2018_PXD003881.zip"
outpath <- "results_modelling_all"
#outpath <- "results_modelling_WHO_noSex"

inputAnntation <- "../samples/annotationIonstar.xlsx"
assign("lfq_write_format", "xlsx", envir = .GlobalEnv)


# MQPeptides<- "D:/Dropbox/DataAnalysis/p2109_PEPTIDE_Analysis/data/MQWorkunit.zip"
# unz(MQPeptides,"modificationSpecificPeptides.txt")
# read.csv(unz(MQPeptides,"modificationSpecificPeptides.txt"),
#         header = TRUE, sep="\t", stringsAsFactors = FALSE)


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

#modelName  <- "Model"
memodel <- "~ dilution +  (1|peptide_Id) + (1|sampleName)"
rlmpep <- "~ dilution. +  peptide_Id"
lmmodel <- "~ dilution."


DEBUG <- FALSE

Contrasts <- c(
  "dilution_(9/3)_3" =   "dilution.e - dilution.a",
  "dilution_(9/4.5)_2" =   "dilution.e - dilution.b",
  "dilution_(9/6)_1.5" =   "dilution.e - dilution.c",
  "dilution_(9/7.5)_1.2" =   "dilution.e - dilution.d",

  "dilution_(7.5/3)_2.5" =   "dilution.d - dilution.a",
  "dilution_(7.5/4.5)_1.6(6)" =   "dilution.d - dilution.b",
  "dilution_(7.5/6)_1.25" =   "dilution.d - dilution.c",

  "dilution_(6/3)_2" =   "dilution.c - dilution.a",
  "dilution_(6/4.5)_1.3(3)" =   "dilution.c - dilution.b",

  "dilution_(4.5/3)_1.5" =   "dilution.b - dilution.a"
)


if (TRUE) {
  assign("lfq_write_format", "xlsx", envir = .GlobalEnv)
  #source("c:/Users/wolski/prog/LFQService/R/tidyMS_application.R")
  res <- application_set_up_MQ_run(outpath = outpath,
                                   inputMQfile = inputMQfile,
                                   inputAnnotation = inputAnntation,
                                   config = config)


  summarised <- data_pep_to_prot(res$data,
                                 res$config,
                                 res$qc_path)

  summarised <- summarised(DEBUG = TRUE)
  data_c <- summarised$results$pepIntensityNormalized
  config_c <- summarised$results$config_pepIntensityNormalized$clone(deep = TRUE)

  #usethis::use_data(data_c)
  #usethis::use_data(config_c)
  #dataIonstar <- summarised
  ##usethis::use_data(dataIonstar)

} else {

  res <- application_set_up_MQ_run(outpath = outpath,
                                   inputMQfile = inputMQfile,
                                   inputAnnotation = inputAnntation,
                                   config = config$clone(deep = TRUE))


  dd <- LFQService:::.workflow_MQ_filter_peptides_V3(res$data, res$config)
  pepIntensityNormalized <- transform_work_intensity(dd$data, dd$config, log2)


  subset <- pepIntensityNormalized %>% filter(grepl("HUMAN", pepIntensityNormalized$protein_Id))
  dd$config$table$workIntensity
  pepIntensityNormalized <- scale_with_subset(pepIntensityNormalized,subset, dd$config)

  dd$config$table$workIntensity

  hist(pepIntensityNormalized$log2_peptide.intensity)
  hist(pepIntensityNormalized$log2_peptide.intensity_subset_scaled)

  pepIntensityNormalized <-
    pepIntensityNormalized <- pepIntensityNormalized %>%
    dplyr::rename(transformedIntensity = dd$config$table$getWorkIntensity())
  dd$config$table$popWorkIntensity()
  dd$config$table$setWorkIntensity("transformedIntensity")

  res <- list(pepIntensityNormalized = pepIntensityNormalized, config_pepIntensityNormalized = dd$config)
  dataIonstarSubsetNorm <- list()
  dataIonstarSubsetNorm$results <- res
  usethis::use_data(dataIonstarSubsetNorm, overwrite = TRUE)

 }

message("######################## fit mixed #######################")

mean(is.na(summarised$results$pepIntensityNormalized$transformedIntensity))

memodel_full <- paste0(config_c$table$getWorkIntensity() , memodel)
modelFunction <- make_custom_model_lmer( memodel_full, model_name = "meModel")
reportColumns <- c(
  "statistic",
  "p.value",
  "p.value.adjusted")

#foo
#source("c:/Users/wolski/prog/LFQService/R/tidyMS_application.R")
if (TRUE) {
  resXXmixmodel <- application_run_modelling_V2(
    outpath = outpath,
    data = data_c,
    pepConfig = config_c,
    modelFunction = modelFunction,
    contrasts = Contrasts,
    modelling_dir = "modelling_results_peptide")
  resXXmixmodel <- resXXmixmodel(do = "result")
}


message("################## fit medpolish ######################")
prot <- summarised$protintensity_fun("unnest")
model <- paste0(prot$config$table$getWorkIntensity() , lmmodel)

modelFunction <- make_custom_model_lm( model, model_name = "Model")
if (TRUE) {
  resXXmedpolish <- application_run_modelling_V2(
    outpath = outpath,
    data = prot$data,
    pepConfig = prot$config,
    modelFunction = modelFunction,
    contrasts = Contrasts,
    modelling_dir = "modelling_results_peptide")
  resXXmedpolish <- resXXmedpolish(do = "result")
}


message("###################### fit ROPECA #######################")
model <- paste0(summarised$results$config_pepIntensityNormalized$table$getWorkIntensity()  , lmmodel)
summarised$results$config_pepIntensityNormalized$table$hierarchyLevel <- 2
modelFunction <- make_custom_model_lm( model, model_name = "pepModel")


if (TRUE) {
  resXXRopeca <- application_run_modelling_V2(
    outpath = outpath,
    data = summarised$results$pepIntensityNormalized,
    pepConfig = summarised$results$config_pepIntensityNormalized,
    modelFunction = modelFunction,
    contrasts = Contrasts,
    modelling_dir = "modelling_results_peptide")
  resXXRopeca <- resXXRopeca(do = "result")
}

detach("package:LFQService",unload = TRUE)
library(LFQService)
ropeca_P <- summary_ROPECA_median_p.scaled(resXXRopeca,contrast = "contrast")


allresults <- list(ropeca_P = ropeca_P, resXXmedpolish = resXXmedpolish, resXXmixmodel = resXXmixmodel )
#saveRDS(allresults, file = "allresults.Rds")

tmp <- contrasts_linfct_vis(ropeca_P,columns = c("beta.based.significance"),
                            estimate = "median.estimate",
                            contrast = "contrast",
                            modelName = "protModelRopeca")
contrasts_linfct_vis_write(tmp, path = file.path(outpath,"modelling_results_peptide"))

