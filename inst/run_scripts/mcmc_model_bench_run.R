rm(list=ls())
library(LFQService)
library(tidyverse)

LFQServiceAnalysisTemplate::setglobalLocation()
thisdir <- file.path(projectLocation,
                     thisPackageName,
                     "inst/MQ_Ionstar2018_PXD003881")
setwd(thisdir)

inputMQfile <- "../samples/p2558_o5748_peptides.zip"


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
annotation

config$table$factors[["dilution_"]] = "sample"
config$table$factors[["run_ID"]] = "run_ID"


config$table$factorLevel <- 1

config$order_Id = "IonStar"
config$project_Id = "p3000"
config$workunit_Id = "IonStar"

# specify model definition

modelName  <- "Model"
memodel <- "~ dilution_ +  (1|peptide_Id)"
lmmodel <- "~ dilution_"

DEBUG <- FALSE

Contrasts <- c(
  "dilution_(9/3)_3" =   "dilution_e - dilution_a",
  "dilution_(9/4.5)_2" =   "dilution_e - dilution_b",
  "dilution_(9/6)_1.5" =   "dilution_e - dilution_c",
  "dilution_(9/7.5)_1.2" =   "dilution_e - dilution_d",

  "dilution_(7.5/3)_2.5" =   "dilution_d - dilution_a",
  "dilution_(7.5/4.5)_1.6(6)" =   "dilution_d - dilution_b",
  "dilution_(7.5/6)_1.25" =   "dilution_d - dilution_c",

  "dilution_(6/3)_2" =   "dilution_c - dilution_a",
  "dilution_(6/4.5)_1.3(3)" =   "dilution_c - dilution_b",

  "dilution_(4.5/3)_1.5" =   "dilution_b - dilution_a"
)


if(TRUE){
  assign("lfq_write_format", "xlsx", envir = .GlobalEnv)
  #source("c:/Users/wolski/prog/LFQService/R/tidyMS_application.R")
  res <- application_set_up_MQ_run(outpath = outpath,
                                   inputMQfile = inputMQfile,
                                   inputAnnotation = inputAnntation,
                                   config = config)

  summarised <- application_summarize_data_pep_to_prot(res$data,
                                                       res$config,
                                                       res$qc_path,
                                                       DEBUG = DEBUG, WRITE_PROTS=FALSE)
  saveRDS(summarised, file="aaa_summarized.RDA")
}else{
  summarised <- readRDS("aaa_summarized.RDA")
}

message("######################## fit mixed #######################")
memodel <- paste0(summarised$results$config_pepIntensityNormalized$table$getWorkIntensity() , memodel)
modelFunction <- make_custom_model_lmer( memodel, model_name = "meModel")
reportColumns <- c("p.value",
                   "p.value.adjusted")


#source("c:/Users/wolski/prog/LFQService/R/tidyMS_application.R")
if(TRUE){
  resXXmixmodel <- application_run_modelling_V2(
    outpath = outpath,
    data = summarised$results$pepIntensityNormalized,
    pepConfig = summarised$results$config_pepIntensityNormalized,
    modelFunction = modelFunction,
    contrasts = Contrasts,
    modelling_dir = "modelling_results_peptide")
  resXXmixmodel <- resXXmixmodel(do="result")
  #saveRDS(resXX, file="resXX.rda")
}



message("################## fit medpolish ######################")
prot <- summarised$protintensity_fun("unnest")
model <- paste0(prot$config$table$getWorkIntensity() , lmmodel)

modelFunction <- make_custom_model_lm( model, model_name = "Model")
if(TRUE){
  resXXmedpolish <- application_run_modelling_V2(
    outpath = outpath,
    data = prot$data,
    pepConfig = prot$config,
    modelFunction = modelFunction,
    contrasts = Contrasts,
    modelling_dir = "modelling_results_peptide")
  resXXmedpolish <- resXXmedpolish(do="result")
}


message("###################### fit ROPECA #######################")
model <- paste0(summarised$results$config_pepIntensityNormalized$table$getWorkIntensity()  , lmmodel)
summarised$results$config_pepIntensityNormalized$table$hierarchyLevel <- 2
modelFunction <- make_custom_model_lm( model, model_name = "pepModel")


if(TRUE){
  resXXRopeca <- application_run_modelling_V2(
    outpath = outpath,
    data = summarised$results$pepIntensityNormalized,
    pepConfig = summarised$results$config_pepIntensityNormalized,
    modelFunction = modelFunction,
    contrasts = Contrasts,
    modelling_dir = "modelling_results_peptide")
  resXXRopeca <- resXXRopeca(do="result")

}

detach("package:LFQService",unload=TRUE)
library(LFQService)
ropeca_P <- summary_ROPECA_median_p.scaled(resXXRopeca,contrast = "contrast")


allresults <- list(ropeca_P = ropeca_P, resXXmedpolish = resXXmedpolish, resXXmixmodel = resXXmixmodel )
saveRDS(allresults, file="allresults.Rds")
tmp <- contrasts_linfct_vis(ropeca_P,columns = c("beta.based.significance"),
                            estimate = "median.estimate",
                            contrast = "contrast",
                            modelName = "protModelRopeca")
contrasts_linfct_vis_write(tmp, path = file.path(outpath,"modelling_results_peptide"))


relevantParameters <- list(outpath = outpath,
                           inputMQfile = inputMQfile,
                           modelling_dir = "modelling_results_peptide",
                           workunit_Id = config$workunit_Id,
                           annotation = annotation,
                           reportColumns = reportColumns,
                           config = config,
                           model= model,
                           Contrasts = Contrasts,
                           project_Id = config$project_Id,
                           order_Id = config$order_Id,
                           author = "Witold Wolski <wew@fgcz.ethz.ch>"
)

LFQService::copy_mixed_model_analysis_script()
rmarkdown::render("mixed_model_analysis_script_Report.Rmd",
                  params= list(pars = relevantParameters),
                  output_format = "html_document",
                  output_dir = outpath,
                  output_file = "index.html")



