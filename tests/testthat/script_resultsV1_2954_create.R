rm(list=ls())

library(readr)
library(tidyverse)
library(LFQService)
library(tidyr)
library(dplyr)

flevel <- 1
path <- "results_FULL_Phonix_Filter"

# resPepProtAnnot <- read_csv(file = "c:/Users/wewol/Dropbox/DataAnalysis/p2954_MSC_IVD_Christina/data/annotatedPeptide_PhonixDS_1097969.csv")
# resPepProtAnnot %>% dplyr::select(top_protein, protein.group.id) %>%
#   distinct() %>% dplyr::sample_n(size=20) -> proteinSel
# resPepProtAnnot <-  dplyr::inner_join(proteinSel, resPepProtAnnot)
# resPepProtAnnot_p2954 <- resPepProtAnnot
# usethis::use_data(resPepProtAnnot_p2954, overwrite=TRUE)

resPepProtAnnot <- LFQService::resPepProtAnnot_p2954
resPepProtAnnot$isotope <- "light"

createMQProteinPeptideConfiguration <- function(ident_qValue = "pep",
                                                intensity = "peptide.intensity",
                                                isotopeLabel = "isotope"){
  atable <- AnalysisTableAnnotation$new()
  atable$fileName = "raw.file"
  # measurement levels.
  atable$hierarchy[["protein_Id"]] <- c("top_protein","protein.group.id")
  atable$hierarchy[["peptide_Id"]] <- c("sequence","peptide.id")
  #
  atable$ident_qValue = ident_qValue
  atable$setWorkIntensity(intensity)
  atable$isotopeLabel = isotopeLabel
  atable$factors[["Condition"]] = "condition"
  atable$factors[["patient_id"]] = "patient_id"
  atable$factors[["runId"]] = "runId"
  anaparam <- AnalysisParameters$new()
  configuration <- AnalysisConfiguration$new(atable, anaparam)
  return(configuration)
}


config <- createMQProteinPeptideConfiguration()
config$table$factorLevel <- 2
resPepProtAnnot %>% filter(reverse == FALSE) -> resPepProtAnnot
head(resPepProtAnnot)
config$table$hierarchy


resDataStart <- setup_analysis(resPepProtAnnot, config)


resDataStart <- remove_small_intensities(resDataStart, config) %>% completeCases(config)
resDataStart <- LFQService::make_interaction_column_config(resDataStart, config)

#LFQService::render_MQSummary_rmd(resDataStart, config , dest_path = path,  workdir=".")
#rmarkdown::render("MQSummary.Rmd", params = list(data = resDataStart, configuration=config))
#rmarkdown::render("MQSummary.Rmd", params=list(data = resDataStart, configuration=config$clone(deep=TRUE)), envir = new.env())


x3 <- summarizeHierarchy(resDataStart, config)
x3 %>% dplyr::inner_join(resDataStart, by="protein_Id") -> resDataStart

# Start filtering
config$table$factorLevel <-flevel


results <- workflow_MQ_protoV1(resDataStart, config, path ,
                                           peptideFilterFunction = LFQService:::.workflow_MQ_filter_peptides_V2 )

protintensity <- LFQService::workflow_MQ_protein_quants( results )



#readr::write_csv(protintensity$data,
#                 path = file.path(path,"transformed_ProteinIntensities.csv"))

#rmarkdown::render("Summarize_Filtering.Rmd", params=results, envir = new.env())
#rmarkdown::render("Summarize_Filtering.Rmd", params=results)

#render_SummarizeFiltering_rmd(results, dest_path=path, dest_file_name = "SummarizeFiltering.pdf",workdir = getwd())

#params <- results
#rmarkdown::render("RunAnalysis_WithParams.Rmd", params = results)
#file.copy("RunAnalysis_WithParams.pdf", file.path(path, "RunAnalysis_WithParams.pdf"), overwrite = TRUE)
results$HEATMAP <-TRUE
resultsV12954 <- results
#usethis::use_data(resultsV12954, overwrite=TRUE)
#saveRDS(results, file = "allData_PhonixDS_1097969.Rds")
names(results)
results$path

res <- workflow_MQ_protoV1_vis(results)
res$figs_raw[[1]]
res$figs_normalized[[1]]
