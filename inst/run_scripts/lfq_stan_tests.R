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
annotation <- annotation %>% filter(annotation$SCI!="un")

config$table$factors[["drug"]] = "genotype"
config$table$factors[["SCI"]] = "SCI"
config$table$factorLevel <- 2

config$order_Id = "o5748"
config$project_Id = "p2558"
config$workunit_Id = "20191120_MQ_repack.zip"


# specify model definition
modelName  <- "Model"
memodel <- "~ drug * SCI  + (1|peptide_Id)"
lmmodel <- "~ drug * SCI"

DEBUG <- TRUE
RUN_ALL <- TRUE

Contrasts <- c("8wk_vs_1wk" = "SCI8wk - SCI1wk",
               "t_vs_v" = "drugt - drugv",
               "t_vs_v_given_8wk" = "`drugt:SCI8wk` - `drugv:SCI8wk`",
               "t_vs_v_given_1wk" = "`drugt:SCI1wk` - `drugv:SCI1wk`",
               "interaction_construct_with_time" = "t_vs_v_given_8wk - t_vs_v_given_1wk"
)

if(TRUE){
  assign("lfq_write_format", "xlsx", envir = .GlobalEnv)
  #source("c:/Users/wolski/prog/LFQService/R/tidyMS_application.R")

  res <- application_set_up_MQ_run(outpath = outpath,
                                   inputMQfile = inputMQfile,
                                   inputAnnotation = annotation,
                                   config=config,
                                   id_extractor = NULL,
                                   use="peptides")

  summarised <- data_pep_to_prot(res$data,
                                 res$config,
                                 res$qc_path)
  summarised <- summarised(DEBUG=TRUE)
  saveRDS(summarised,"aaa_summarized.RDA")
}else{
  summarised <- readRDS("aaa_summarized.RDA")
}

message("######################## fit mixed #######################")
memodel <- paste0(summarised$results$config_pepIntensityNormalized$table$getWorkIntensity() , memodel)
modelFunction <- make_custom_model_lmer( memodel, model_name = "meModel")
modelFunction


#source("c:/Users/wolski/prog/LFQService/R/tidyMS_application.R")
if(TRUE){
  resXXmixmodel <- application_run_modelling_V2(
    outpath = outpath,
    data = summarised$results$pepIntensityNormalized,
    pepConfig = summarised$results$config_pepIntensityNormalized,
    modelFunction = modelFunction,
    contrasts = Contrasts,
    modelling_dir = "modelling_results_peptide")

  tmp <- resXXmixmodel(DEBUG=TRUE)
  dd <- tmp$res_contrasts(DEBUG=TRUE)
}

# Work on brms code

data <- summarised$results$pepIntensityNormalized
config <- summarised$results$config_pepIntensityNormalized
nested <- data %>% group_by_at(config$table$hkeysLevel()) %>% nest()

library(brms)



#library(snow)
#cl <- makeCluster(4)
startmodel <- brms::brm(memodel, mdata, cores=6)

nestedT <- nested[1:10,]
res <- nestedT %>% mutate(summary = purrr::map( data, ms_brms_model, startmodel, dd$linfct_A))
res <- res %>% filter(!sapply(summary, is.null))
res$summary[[1]]
res %>% select(protein_Id, summary) %>% unnest()

as_tibble(res$summary[[1]], rownames="contrast")
"transformedIntensity~ drug * SCI  + (1|peptide_Id)"
