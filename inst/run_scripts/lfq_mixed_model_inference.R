rm(list=ls())
library(LFQService)
library(tidyverse)

DEBUG <- FALSE


outpath <- "results_modelling"
inputMQfile <-  "data/MQresults-20190805_combined_txt.zip"
inputAnntation <- "annotation.xlsx"
assign("lfq_write_format", "xlsx", envir = .GlobalEnv)


# creates default configuration
config <- LFQService::create_MQ_peptide_Configuration()

annotation <- readxl::read_xlsx(inputAnntation)

config$table$factors[["genotype_"]] = "genotype"
config$table$factors[["SCI_"]] = "SCI"
config$table$factors[["wasIstDas"]] = "wasIstDas"
config$table$factorLevel <- 2


# specify model definition
modelName  <- "f_genotype_SCI"
model <- "~ genotype_*SCI_  + (1|peptide_Id)"
is.mixed <- TRUE
DEBUG <- TRUE

Contrasts <- c("8w - nv" = "SCI_8w - SCI_nv",
               "7d - nv" = "SCI_7d - SCI_nv" ,
               "tlr4 - wt" = "genotype_tlr4 - genotype_wt",
               "tlr4_vs_wt_given_nv" = "`genotype_tlr4:SCI_nv` - `genotype_wt:SCI_nv`",
               "tlr4_vs_wt_given_7d" = "`genotype_tlr4:SCI_7d` - `genotype_wt:SCI_7d`",
               "tlr4_vs_wt_given_8w" = "`genotype_tlr4:SCI_8w` - `genotype_wt:SCI_8w`",
               "interaction_genotype_with_7d_nv" = "tlr4_vs_wt_given_7d - tlr4_vs_wt_given_nv",
               "interaction_genotype_with_8w_nv" = "tlr4_vs_wt_given_8w - tlr4_vs_wt_given_nv"
)


assign("lfq_write_format", "xlsx", envir = .GlobalEnv)

# specify model definition
modelName  <- "f_genotype_SCI"
model <- "~ genotype_*SCI_  + (1|peptide_Id)"

reportColumns <- c("p.value",
                   "p.value.adjusted")


res <- application_set_up_MQ_run(outpath = outpath,
                                 inputMQfile = inputMQfile,
                                 inputAnntation = inputAnntation)

summarised <- application_summarize_data(res$data,
                                         res$config,
                                         res$qc_path,
                                         DEBUG = DEBUG)



model <- paste0(summarised$results$config_pepIntensityNormalized$table$getWorkIntensity() ,model)
modelFunction <- make_custom_model_lmer( model, model_name = "Model")

source("c:/Users/wolski/prog/LFQService/R/tidyMS_application.R")
res <- application_run_modelling(outpath = outpath,
                                 protIntensityNormalized = summarised$results$pepIntensityNormalized,
                                 pepConfig = summarised$results$config_pepIntensityNormalized,
                                 modelFunction = modelFunction,
                                 Contrasts = Contrasts,
                                 modelling_dir = "modelling_results_peptide")


relevantParameters <- list(outpath = outpath,
                           inputMQfile = inputMQfile,
                           modelling_dir = "modelling_results_peptide",
                           workunit_name = "MaxQuant_p3147_o5997",
                           annotation = annotation,
                           reportColumns = reportColumns,
                           config = config,
                           model= model,
                           Contrasts = Contrasts,
                           project_id = "p3147",
                           order_id = "o5997",
                           author = "Witold Wolski <wew@fgcz.ethz.ch>"
)

copy_mixed_model_analysis_script()
rmarkdown::render("mixed_model_analysis_script_Report.Rmd",
                  params= list(pars = relevantParameters),
                  output_format = "html_document",
                  output_dir = outpath,
                  output_file = "index.html")



