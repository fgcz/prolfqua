rm(list = ls())
library(LFQService)
library(tidyverse)

outpath <- "results_modelling_protAggregate_msstats"
inputMQfile <-  "../samples/timstof/MSstats.zip"

# massaging of the input.
inputFile <- readr::read_csv(unz(inputMQfile, filename = "MSstats.csv"))
inputFile$BioReplicate <- paste("br", inputFile$BioReplicate, sep = "")
inputFile$Condition <- make.names(inputFile$Condition)
inputFile$pep <- 0

# setup modelling
config <- create_MSstats_MSFragger_config()
config$table$factors[["Celline_"]] <- "Condition"
config$table$factors[["BioReplicate_"]] <- "BioReplicate"
config$table$hierarchyLevel <- 1 # max 2 - for plotting (heatmaps QC etc)

memodel <- "~ Celline_"
# repeated measurements
#memodel <- "~ Celline_ + (1|BioReplicate_)"
# techreps
#memodel <- "~ Celline_ + (1|BioReplicate_) + (1|Run)"
# Factorial design
# memodel <- "~ Celline_ * BioReplicate_"

# set contrasts
Contrasts <- c("RKOvsRKO.R" = "Celline_RKO - Celline_RKO.R")

# Bookkeeping
config$order_Id = "o1"
config$project_Id = "p3061"
config$workunit_Id = "Resource : 1577525 - 20200317_TK.zip"

##### Boilerplate code

annotation <- inputFile %>%
  dplyr::select(Run, Condition, BioReplicate) %>%
  distinct()

qcdir = "qc_results"
qc_path <- file.path(outpath, qcdir )
if (!dir.exists(qc_path)) {
  dir.create(qc_path, recursive = TRUE)
}

dirlayout <- list()
dirlayout$outdir <- outpath
dirlayout$qc_path = qc_path
#####

assign("lfq_write_format", "xlsx", envir = .GlobalEnv)
resDataStart <- setup_analysis( inputFile, config)
resDataStart <- remove_small_intensities( resDataStart, config, threshold = 4 ) %>%  complete_cases(config)

res <- list()
res$data <- resDataStart
res$config <- config
res$qc_path <- dirlayout$qc_path

summarised <- data_pep_to_prot(res$data,
                               res$config,
                               res$qc_path)

saveRDS(summarised , file = "summarised.Rds")


if (TRUE) {
  .Device
  summarised("render")
  .Device
  summarised("pepwrite")
  .Device
  summarised("protwrite")
  .Device
  #summarised("plotprot")
  saveRDS(summarised,"aaa_summarized.RDA")
}

summarised <- summarised(DEBUG = TRUE)
protIntensity <- summarised$protintensity_fun("unnest")

message("######################## fit mixed #######################")
#reportColumns <- c("p.value",
#                   "p.value.adjusted",
#                   )

memodel <- paste0(protIntensity$config$table$getWorkIntensity() , memodel)
modelFunction <- make_custom_model_lm( memodel, model_name = "Model")

unique(protIntensity$data$Condition_)

if (TRUE) {
  resXXmixmodel <- application_run_modelling_V2(
    outpath = outpath,
    data = protIntensity$data,
    config = protIntensity$config,
    modelFunction = modelFunction,
    contrasts = Contrasts,
    modelling_dir = "modelling_results_peptide")

  #names(modelFunction)
  xx <- resXXmixmodel(DEBUG = TRUE)
  dd <- xx$modellingResult_fun()
  dd$modellingResult$modelProtein
  m <- dd$modellingResult$modelProtein$linear_model[[1]]
  #modelFunction$isSingular(m)

  bb <- resXXmixmodel(do = "write_modelling")
  bb <- resXXmixmodel(do = "write_contrasts")
}
reportColumns <- ""
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



