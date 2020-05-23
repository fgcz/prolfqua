rm(list = ls())

library(LFQService)
library(tidyverse)
library(dplyr)

inputMQfile <- "../samples/test_MQ_IonStar2018_PXD003881.zip"
outpath <- "results_modelling"
#outpath <- "results_modelling_WHO_noSex"

inputAnntation <- "../samples/annotationIonstar.xlsx"
assign("lfq_write_format", "xlsx", envir = .GlobalEnv)


# MQPeptides<- "D:/Dropbox/DataAnalysis/p2109_PEPTIDE_Analysis/data/MQWorkunit.zip"
# unz(MQPeptides,"modificationSpecificPeptides.txt")
# read.csv(unz(MQPeptides,"modificationSpecificPeptides.txt"),
#         header=TRUE, sep="\t", stringsAsFactors = FALSE)


# creates default configuration
config <- LFQService::create_MQ_peptide_Configuration()
annotation <- readxl::read_xlsx(inputAnntation)

config$table$factors[["dilution."]] = "sample"
config$table$factors[["run_ID"]] = "run_ID"

config$table$factorDepth <- 1

config$order_Id = "IonStar"
config$project_Id = "p3000"
config$workunit_Id = "IonStar"

# specify model definition

modelName  <- "Model"
memodel <- "~ dilution. +  (1|peptide_Id) + (1|sampleName)"
memodel_trunc <- "|trunc() ~ dilution. +  (1|peptide_Id) + (1|sampleName)"

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
  saveRDS(summarised, file = "aaa_summarized.RDA")

}else{
  summarised <- readRDS("aaa_summarized.RDA")
}


message("######################## fit mixed #######################")
memodel_full <- paste0(summarised$results$config_pepIntensityNormalized$table$getWorkIntensity() , memodel)
modelFunction <- make_custom_model_lmer( memodel_full, model_name = "meModel")

#source("c:/Users/wolski/prog/LFQService/R/tidyMS_application.R")
if (TRUE) {
  resXXmixmodel <- application_run_modelling_V2(
    outpath = outpath,
    data = summarised$results$pepIntensityNormalized,
    pepConfig = summarised$results$config_pepIntensityNormalized,
    modelFunction = modelFunction,
    contrasts = Contrasts,
    modelling_dir = "modelling_results_peptide")

  tmp <- resXXmixmodel(DEBUG = TRUE)
  dd <- tmp$res_contrasts(DEBUG = TRUE)
}

# Work on brms code

library(brms)
data <- summarised$results$pepIntensityNormalized
config <- summarised$results$config_pepIntensityNormalized


nested <- data %>% group_by_at(config$table$hkeysDepth()) %>% nest()

mdata2 <- nested$data[[2]]
mdata1 <- nested$data[[1]]
mdata26 <- nested$data[[26]]

startmodel <- brms::brm(memodel_full, mdata2, cores = 6, refresh = 0)
tmp <- update(startmodel, newdata = mdata26)

source("../../R/tidyMS_stanr.R")
library(MCMCvis)

res <- ms_brms_model(mdata = mdata2,
              memodel = startmodel,
              fixef = config$table$fkeysDepth(),
              linfct_A = dd$linfct_A)

if (TRUE) {
res <- nested %>% mutate(summary =
                           purrr::map( data, ms_brms_model,
                                       startmodel,
                                       config$table$fkeysDepth(),
                                       dd$linfct_A)
                         )
saveRDS(res, file = "rstandSimpleMixed.RDS")
}

if (FALSE) {
  memodel_trunc_full <- paste0(summarised$results$config_pepIntensityNormalized$table$getWorkIntensity() , memodel_trunc)
  startmodel <- brms::brm(memodel_trunc_full, mdata2, cores = 6, refresh = 0)
  res <- nested %>% mutate(summary =
                             purrr::map( data, ms_brms_model,
                                         startmodel,
                                         config$table$fkeysDepth(),
                                         dd$linfct_A)
  )
  saveRDS(res, file = "rstandTruncMixed.RDS")
}
