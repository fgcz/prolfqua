rm(list=ls())
library(LFQService)
library(tidyverse)



outpath <- "results_modelling_testing"

inputMQfile <-  "../samples/p2558_05748_modspec.zip"
inputMQfile <-  "../samples/p2558_o5748_peptides.zip"


inputAnntation <- "../samples/p2558_05748_annotation.xlsx"
assign("lfq_write_format", "xlsx", envir = .GlobalEnv)


# creates default configuration
config <- LFQService::create_MQ_peptide_Configuration()

annotation <- readxl::read_xlsx(inputAnntation)
annotation <- annotation %>% filter(annotation$SCI != "un")

config$table$factors[["drug_"]] = "genotype"
config$table$factors[["SCI_"]] = "SCI"
config$table$factorLevel <- 2

config$order_Id = "o5748"
config$project_Id = "p2558"
config$workunit_Id = "20191120_MQ_repack.zip"


# specify model definition
modelName  <- "Model"
memodel <- "~ drug_ * SCI_  + (1|peptide_Id)"
lmmodel <- "~ drug_ * SCI_"

DEBUG <- TRUE
RUN_ALL <- TRUE

Contrasts <- c("8wk_vs_1wk" = "SCI_8wk - SCI_1wk",
               "t_vs_v" = "drug_t - drug_v",
               "t_vs_v_given_8wk" = "`drug_t:SCI_8wk` - `drug_v:SCI_8wk`",
               "t_vs_v_given_1wk" = "`drug_t:SCI_1wk` - `drug_v:SCI_1wk`",
               "interaction_construct_with_time" = "t_vs_v_given_8wk - t_vs_v_given_1wk"
)

if (TRUE) {
  assign("lfq_write_format", "xlsx", envir = .GlobalEnv)
  #source("c:/Users/wolski/prog/LFQService/R/tidyMS_application.R")

  res <- application_set_up_MQ_run(outpath = outpath,
                                   inputMQfile = inputMQfile,
                                   inputAnnotation = annotation,
                                   config = config,
                                   id_extractor = NULL,
                                   use = "peptides")


  summarised <- application_summarize_data_pep_to_prot(res$data,
                                                       res$config,
                                                       res$qc_path,
                                                       DEBUG = FALSE,
                                                       WRITE_PROTS = FALSE)


  summarised("render")
  summarised("plotprot")
  summarised("pepwrite")
  summarised("protwrite")

  saveRDS(summarised,"aaa_summarized.RDA")
}else{
  summarised <- readRDS("aaa_summarized.RDA")
}

message("######################## fit mixed #######################")
memodel <- paste0(summarised$results$config_pepIntensityNormalized$table$getWorkIntensity() , memodel)
modelFunction <- make_custom_model_lmer( memodel, model_name = "meModel")
reportColumns <- c("p.value",
                   "p.value.adjusted")


#source("c:/Users/wolski/prog/LFQService/R/tidyMS_application.R")
if (TRUE) {
  resXXmixmodel <- application_run_modelling_V2(
    outpath = outpath,
    data = summarised$results$pepIntensityNormalized,
    pepConfig = summarised$results$config_pepIntensityNormalized,
    modelFunction = modelFunction,
    contrasts = Contrasts,
    modelling_dir = "modelling_results_peptide")

  resXXmixmodel(do = "write_modelling")
  resXXmixmodel(do = "write_contrasts")

}




relevantParameters <- list(outpath = outpath,
                           inputMQfile = inputMQfile,
                           modelling_dir = "modelling_results_peptide",
                           workunit_Id = config$workunit_Id,
                           annotation = annotation,
                           reportColumns = reportColumns,
                           config = config,
                           model = memodel,
                           Contrasts = Contrasts,
                           project_Id = config$project_Id,
                           order_Id = config$order_Id,
                           author = "Witold Wolski <wew@fgcz.ethz.ch>"
)

LFQService::copy_mixed_model_analysis_script()
rmarkdown::render("mixed_model_analysis_script_Report.Rmd",
                  params = list(pars = relevantParameters),
                  output_format = "html_document",
                  output_dir = outpath,
                  output_file = "index.html")



