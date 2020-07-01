rm(list = ls())

library(readr)
library(tidyverse)
library(LFQService)
library(tidyr)
library(dplyr)

flevel <- 1
path <- "results_FULL_Phonix_Filter"


fp <- .find.package.file("LFQService", "samples/testdata/annotatedPeptide_PhonixDS_1097969.csv")

resPepProtAnnot <- read_csv(file = fp)

resPepProtAnnot$isotope <- "light"

createMQProteinPeptideConfiguration <-
  function(ident_qValue = "pep",
           intensity = "peptide.intensity",
           isotopeLabel = "isotope") {
    atable <- AnalysisTableAnnotation$new()
    atable$fileName = "raw.file"
    # measurement levels.
    atable$hierarchy[["protein_Id"]] <-
      c("top_protein", "protein.group.id")
    atable$hierarchy[["peptide_Id"]] <- c("sequence", "peptide.id")
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
config$table$factorDepth <- 2
resPepProtAnnot %>% dplyr::filter(reverse == FALSE) -> resPepProtAnnot


resDataStart <- setup_analysis(resPepProtAnnot, config)


resDataStart <- remove_small_intensities(resDataStart, config) %>%
  complete_cases(config)
resDataStart <-
  LFQService::make_interaction_column(resDataStart, config$table$fkeysDepth())

# LFQService::render_MQSummary_rmd(resDataStart, config , dest_path = path,  workdir=".")
# rmarkdown::render("MQSummary.Rmd", params = list(data = resDataStart, configuration=config))
# rmarkdown::render("MQSummary.Rmd", params=list(data = resDataStart, configuration=config$clone(deep=TRUE)), envir = new.env())

# Start filtering
config$table$factorDepth <- flevel


resFilt <-
  LFQService::filter_proteins_by_peptide_count(resDataStart, config)


results <- LFQService::normalize_log2_robscale(resFilt$data, config)

protintensity <- LFQService::medpolish_protein_quants(
  results$data,results$config)

protintensity <- protintensity("unnest")
LFQService::toWideConfig(protintensity$data, protintensity$config)


results <- list(
  dataStart = LFQData$new(resDataStart, config),
  dataEnd = LFQData$new(resFilt$data, config)
)

render_SummarizeFiltering_rmd(results, workdir = ".")

# dataStart = LFQData$new(resDataStart, config)
# dataEnd = LFQData$new(resFilt$data, resFilt$config)

rmarkdown::render("Summarize_Filtering.Rmd",
                  params = results,
                  envir = new.env(),
                  output_format = "html_document")


rmarkdown::render("Summarize_Filtering.Rmd", params = results)
#render_SummarizeFiltering_rmd(results, dest_path=path, dest_file_name = "SummarizeFiltering.pdf",workdir = getwd())


#params <- results
#rmarkdown::render("RunAnalysis_WithParams.Rmd", params = results)
#file.copy("RunAnalysis_WithParams.pdf", file.path(path, "RunAnalysis_WithParams.pdf"), overwrite = TRUE)
results$HEATMAP <- TRUE
resultsV12954 <- results
#usethis::use_data(resultsV12954, overwrite=TRUE)
#saveRDS(results, file = "allData_PhonixDS_1097969.Rds")
# names(results)
# results$path
#
# res <-
#   plot_hierarchies_line_df(results$pepIntensityNormalized,
#                            results$config_pepIntensityNormalized)
# res[[1]]
#
# res <-
#   plot_hierarchies_boxplot_df(results$pepIntensityNormalized,
#                               results$config_pepIntensityNormalized)
# res[[1]]
#
