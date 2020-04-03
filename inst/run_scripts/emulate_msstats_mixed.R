rm(list = ls())
library(LFQService)
library(tidyverse)

###### Setting up output dir
outpath <- "results_modelling_pepAggregate_msstats"
qcdir = "qc_results"

qc_path <- file.path(outpath, qcdir )
if (!dir.exists(qc_path)) {
  dir.create(qc_path, recursive = TRUE)
}

dirlayout <- list()
dirlayout$outdir <- outpath
dirlayout$qc_path = qc_path
#####


inputMQfile <-  "../samples/timstof/MSstats.zip"

inputFile <- readr::read_csv(unz(inputMQfile, filename = "MSstats.csv"))
inputFile$BioReplicate <- paste("br", inputFile$BioReplicate, sep = "")
inputFile$Condition <- gsub("-",".",inputFile$Condition)

unique(inputFile$Condition)

pepInputFile <- inputFile %>%
  dplyr::group_by(ProteinName , PeptideSequence , IsotopeLabelType, Condition, BioReplicate, Run) %>%
  dplyr::summarize( pepIntensity = sum(Intensity, na.rm = TRUE)) %>% dplyr::ungroup()

dim(inputFile)/dim(pepInputFile)


#mean(is.na(inputFile$FragmentIon))
#mean(is.na(inputFile$ProductCharge))
pepInputFile$pep <- 0
annotation <- pepInputFile %>% dplyr::select(Run, Condition, BioReplicate) %>% distinct()


msstats_config <- function(
  ident_qValue = "pep",
  intensity = "pepIntensity",
  isotopeLabel = "IsotopeLabelType"){
  atable <- AnalysisTableAnnotation$new()
  atable$fileName = "raw.file"
  # measurement levels.
  atable$hierarchy[["protein_Id"]] <- c("ProteinName")
  atable$hierarchy[["peptide_Id"]] <- c("PeptideSequence")

  atable$fileName = "Run"
  atable$hierarchyLevel <- 1
  #
  atable$ident_qValue = ident_qValue
  atable$setWorkIntensity(intensity)
  atable$isotopeLabel = isotopeLabel

  anaparam <- AnalysisParameters$new()
  anaparam$min_peptides_protein <- 2
  configuration <- AnalysisConfiguration$new(atable, anaparam)
  return(configuration)
}

config <- msstats_config()
config$table$factors[["Condition_"]] <- "Condition"
config$table$factors[["BioReplicate_"]] <- "BioReplicate"


config$order_Id = "o5748"
config$project_Id = "p2558"
config$workunit_Id = "20191120_MQ_repack.zip"


# specify model definition
modelName  <- "Model"
memodel <- "~ Condition_  + (1|peptide_Id) + (1|Run)"

#lmmodel <- "~ drug_ * SCI_"

DEBUG <- TRUE
RUN_ALL <- TRUE

Contrasts <- c("RKOvsRKO.R" = "Condition_RKO - Condition_RKO.R")
assign("lfq_write_format", "xlsx", envir = .GlobalEnv)

resDataStart <- setup_analysis(pepInputFile, config)
unique(resDataStart$Condition_)
resDataStart <- remove_small_intensities( resDataStart, config, threshold = 4 ) %>%  complete_cases(config)

res <- list()
res$data <- resDataStart
res$config <- config
res$qc_path <- dirlayout$qc_path

summarised <- data_pep_to_prot(res$data,
                               res$config,
                               res$qc_path)

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

summarisedData <- summarised(DEBUG = TRUE)
message("######################## fit mixed #######################")
reportColumns <- c("p.value",
                   "p.value.adjusted",
                   "statistic")

memodel <- paste0(summarisedData$results$config_pepIntensityNormalized$table$getWorkIntensity() , memodel)
modelFunction <- make_custom_model_lmer( memodel,
                                         model_name = "Model", report_columns =  reportColumns)

unique(summarised$results$pepIntensityNormalized$Condition_)

if (TRUE) {

  resXXmixmodel <- application_run_modelling_V2(
    outpath = outpath,
    data = summarisedData$results$pepIntensityNormalized,
    config = summarisedData$results$config_pepIntensityNormalized,
    modelFunction = modelFunction,
    contrasts = Contrasts,
    modelling_dir = "modelling_results_peptide")

  undebug(resXXmixmodel)

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



