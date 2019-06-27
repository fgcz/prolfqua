rm(list=ls())

library(readr)
library(tidyverse)
library(LFQService)
library(tidyr)
library(dplyr)

HEATMAP <- TRUE
ALLPROTEINPLOTS <- FALSE
MQSUMMARY<- TRUE
path <- "."


#if(!dir.exists(path)){
#  dir.create(path)
#}



#resMetaboDataProgenesis <- readRDS(file="c:/Users/wewol/prog/LFQService/inst/samples/resMetabo.rda")
#usethis::use_data(resMetaboDataProgenesis, overwrite = TRUE)
resMetabo <- LFQService::resMetaboDataProgenesis
resMetabo <- resMetabo  %>% mutate(NRS = gsub( "NRS","NRS_",NRS))
resMetabo %>% rename(Mortality = Outcome, Intervention = Treatment) -> resMetabo

createMetaboCompoundConfiguration <- function(isotopeLabel="isotope",
                                                  qValue="Score"){
  atable <- AnalysisTableAnnotation$new()
  atable$fileName = "rawname"

  # measurement levels.
  atable$hierarchy[["Compound"]] <- c("Compound","Compound_ID")
  #
  atable$ident_qValue = qValue
  atable$workIntensity = "rawIntensity"
  #atable$workIntensity = "FG.MS1PeakArea"
  atable$isotopeLabel = isotopeLabel

  atable$factors[["Mortality"]] = "Mortality"
  atable$factors[["Intervention"]] = "Intervention"
  atable$factors[["NRS"]] = "NRS"
  #atable$factors[["AgeAtDeath.Days"]] = "AgeAtDeath_Days"
  atable$factors[["run_ID"]] = "run_ID"

  atable$factorLevel = 2
  anaparam <- AnalysisParameters$new()
  anaparam$min_peptides_protein <- 2
  configuration <- AnalysisConfiguration$new(atable, anaparam)
  return(configuration)
}

config <- createMetaboCompoundConfiguration()
precursorData <- setup_analysis(resMetabo, config)
precursorData <- remove_small_intensities(precursorData, config,threshold = 4)

# filter qvalues and aggregate peptides -----
hierarchyCounts(precursorData,config)
# This code should be the same for maxquant ----

resDataStart <- LFQService::make_interaction_column_config(precursorData, config)

#readr::write_csv(resDataStart, path = file.path(path, "annotatedTable_Peptide_RAW_Data.csv"))
#saveRDS(config,file.path(path,"config.Rdata"))

if(MQSUMMARY){
  LFQService::render_METABO_Summary_rmd(resDataStart, config , dest_path = ".",  workdir=".")
}


x3 <- summarizeHierarchy(resDataStart, config)
x3 %>% inner_join(resDataStart, by="Compound") -> resDataStart

############################
# Start filtering
config$table$getWorkIntensity()

results <- list()
results$resDataStart <- completeCases(resDataStart, config)
results$config_resDataStart <- config$clone(deep=TRUE)

filteredPep <- filter_factor_levels_by_missing(resDataStart,
                                               config,
                                               percent =50,
                                               factor_level = 1
)

results$filteredPep <- filteredPep
results$config_filteredPep <- config
results$removedFeatures <- anti_join(results$resDataStart, results$filteredPep, by="Compound")
results$removedFeatures %>% dplyr::select(Compound) %>% distinct()

results$config_dataTransformed <- config$clone(deep=TRUE)
filteredPep <- LFQService::transform_work_intensity(filteredPep, results$config_dataTransformed, log2)
results$dataTransformed <- LFQService::applyToIntensityMatrix(filteredPep,
                                                              results$config_dataTransformed ,
                                                              robust_scale)


#                                                              vsn::justvsn)
#results$config_dataTransformed$table$getWorkIntensity()
#head(results$dataTransformed)

results$HEATMAP <- TRUE
results$path <- "."

#results_MetaboData <- results
usethis::use_data(results_MetaboData)
#saveRDS(results, file="allData.rds")

if(MQSUMMARY){
  #rmarkdown::render("METABO_Summarize_Filtering.Rmd",params=results,  envir = new.env())
  #file.copy("METABO_Summarize_Filtering.pdf", file.path(path,"METABO_Summarize_Filtering.pdf" ), overwrite = TRUE)

  LFQService::render_METABO_SummarizeFiltering_rmd(results , dest_path = ".",  workdir=".")
}

