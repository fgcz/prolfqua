# Run mixed models on benchmark dataset.

rm(list = ls())
library(LFQService)
library(tidyverse)
library(dplyr)

conflicted::conflict_prefer("filter", "dplyr")

inputMQfile <-
  "C:/Users/wewol/MAS_WEW/LFQServiceAnalysisTemplate/inst/benchmarkData/MQ_Ionstar2018_PXD003881.zip"
outpath <- "results_modelling_all"
outpath <- "results_modelling_WHO_noSex"
inputAnnotation <- "C:\\Users\\wewol\\MAS_WEW\\LFQServiceAnalysisTemplate\\inst\\MQ_Ionstar2018_PXD003881/annotationIonstar.xlsx"

pStruct <- ProjectStructure$new(outpath,
                                project_Id = 1,
                                order_Id = 1,
                                workunit_Id = 1,
                                inputAnnotation = inputAnnotation,
                                inputData = inputMQfile)
pStruct$create()



# MQPeptides<- "D:/Dropbox/DataAnalysis/p2109_PEPTIDE_Analysis/data/MQWorkunit.zip"
# unz(MQPeptides,"modificationSpecificPeptides.txt")
# read.csv(unz(MQPeptides,"modificationSpecificPeptides.txt"),
#         header = TRUE, sep="\t", stringsAsFactors = FALSE)


mqdata <- tidyMQ_Peptides_Config(pStruct$inputData)
annotation <- readxl::read_xlsx(inputAnnotation)
head(annotation)

res <- application_add_annotation(
  mqdata$data,
  inputAnnotation = annotation,
  fileName  = mqdata$config$table$fileName
)

# creates default configuration

config <- mqdata$config
config$table$factors[["dilution."]] = "sample"
config$table$factors[["fPairing."]] = "fakePair_Id"
config$table$factors[["run_Id"]] = "run_Id"
config$table$factorDepth <- 1
mqdata$config <- config


# specify model definition

#modelName  <- "Model"
memodel <- "~ dilution. +  (1|peptide_Id) + (1|sampleName)"
#rlmpep <- "~ dilution. +  peptide_Id"
lmmodel <- "~ dilution."

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


data <- setup_analysis(res, mqdata$config)

filteredData <- LFQService::filter_proteins_by_peptide_count(data, mqdata$config)



normalizedData <- normalize_log2_robscale(filteredData$data,
                                          mqdata$config)


protintensity_fun <- medpolish_protein_quants( normalizedData$data,
                                               normalizedData$config )




dataIonstarProtein <- protintensity_fun("unnest")
usethis::use_data(dataIonstarProtein, overwrite = TRUE )


memodel_full <- paste0(dataIonstarNormalizedPep$config$table$getWorkIntensity() , memodel)
modelFunction <- make_custom_model_lmer(memodel_full, model_name = "meModel")
reportColumns <- c("statistic",
                   "p.value",
                   "p.value.adjusted")

#foo
#source("c:/Users/wolski/prog/LFQService/R/tidyMS_application.R")
if (TRUE) {
  resXXmixmodel <- application_run_modelling_V2(
    data = dataIonstarNormalizedPep$data,
    config = dataIonstarNormalizedPep$config,
    modelFunction = modelFunction,
    contrasts = Contrasts,
    modelling_dir = "modelling_results_peptide"
  )
  resXXmixmodel <- resXXmixmodel(do = "result")
}


message("################## fit medpolish ######################")
model <- paste0(dataIonstarProtein$config$table$getWorkIntensity() , lmmodel)

modelFunction <- make_custom_model_lm(model, model_name = "Model")
if (TRUE) {
  resXXmedpolish <- application_run_modelling_V2(
    data = dataIonstarProtein$data,
    config = dataIonstarProtein$config,
    modelFunction = modelFunction,
    contrasts = Contrasts,
    modelling_dir = "modelling_results_peptide"
  )
  resXXmedpolish <- resXXmedpolish(do = "result")
}


message("###################### fit ROPECA #######################")
model <-
  paste0(
    dataIonstarNormalizedPep$config$table$getWorkIntensity()  ,
    lmmodel
  )
  dataIonstarNormalizedPep$config$table$hierarchyDepth <- 2
  modelFunction <- make_custom_model_lm(model, model_name = "pepModel")


if (TRUE) {
  resXXRopeca <- application_run_modelling_V2(
    data = dataIonstarNormalizedPep$data,
    config = dataIonstarNormalizedPep$config,
    modelFunction = modelFunction,
    contrasts = Contrasts,
    modelling_dir = "modelling_results_peptide"
  )
  resXXRopeca <- resXXRopeca(do = "result")
}

ropeca_P <-
  summary_ROPECA_median_p.scaled(resXXRopeca,
                                 contrast = "contrast")

allresults <-
  list(
    ropeca_P = ropeca_P,
    resXXmedpolish = resXXmedpolish,
    resXXmixmodel = resXXmixmodel
  )
#saveRDS(allresults, file = "allresults.Rds")

tmp <-
  contrasts_linfct_vis(
    ropeca_P,
    columns = c("beta.based.significance"),
    estimate = "estimate",
    contrast = "contrast",
    modelName = "protModelRopeca"
  )
contrasts_linfct_vis_write(tmp, path = file.path(outpath, "modelling_results_peptide"))
