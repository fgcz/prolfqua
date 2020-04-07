rm(list = ls())
library(LFQService)
library(tidyverse)

###### Setting up output dir
# specify model definition

inputMQfile <-  "../samples/timstof/MSstats.zip"
inputFile <- readr::read_csv(unz(inputMQfile, filename = "MSstats.csv"))
inputFile$BioReplicate <- paste("br", inputFile$BioReplicate, sep = "")
inputFile$Condition <- gsub("-",".",inputFile$Condition)
inputFile$pep <- 0

annotation <- inputFile %>%
  dplyr::select(Run, Condition, BioReplicate) %>%
  distinct()

## Tell LFQ Service what column is what.
atable <- AnalysisTableAnnotation$new()
# measurement levels.
atable$hierarchy[["protein_Id"]] <- c("ProteinName")
atable$hierarchy[["peptide_Id"]] <- c("PeptideSequence","PrecursorCharge")
atable$fileName = "Run"
atable$ident_qValue = "pep"
atable$setWorkIntensity("Intensity")
atable$isotopeLabel = "IsotopeLabelType"

anaparam <- AnalysisParameters$new()
anaparam$min_peptides_protein <- 2
config <- AnalysisConfiguration$new(atable, anaparam)

memodel <- "~ Condition_"
# how to model repeated measurements
#memodel <- "~ Condition_ + (1|BioReplicate)"
# how to model techreps
#memodel <- "~ condition_ + (1|BioReplicate) + (1|Run)"
config$table$factors[["Condition_"]] <- "Condition"
config$table$factors[["BioReplicate_"]] <- "BioReplicate"
atable$hierarchyLevel <- 1

Contrasts <- c("RKOvsRKO.R" = "Condition_RKO - Condition_RKO.R")


config$order_Id = "o1"
config$project_Id = "p3061"
config$workunit_Id = "Resource : 1577525 - 20200317_TK.zip"


outpath <- "results_modelling_protAggregate_msstats"
qcdir = "qc_results"




qc_path <- file.path(outpath, qcdir )
if (!dir.exists(qc_path)) {
  dir.create(qc_path, recursive = TRUE)
}

dirlayout <- list()
dirlayout$outdir <- outpath
dirlayout$qc_path = qc_path
#####

DEBUG <- TRUE
RUN_ALL <- TRUE


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
  summarised("plotprot")
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



