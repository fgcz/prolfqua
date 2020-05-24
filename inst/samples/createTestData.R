library(LFQService)

fp <- file.path( find.package("LFQService"),"inst/samples/testdata/annotatedPeptide_PhonixDS_1097969.csv")
resPepProtAnnot <- read_csv(file = fp)
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
config$table$factorDepth <- 2
resPepProtAnnot %>% dplyr::filter(reverse == FALSE) -> resPepProtAnnot


resDataStart <- setup_analysis(resPepProtAnnot, config)
resDataStart <- remove_small_intensities(resDataStart, config) %>%
  complete_cases(config)

smallData <- resDataStart %>% dplyr::filter(protein_Id %in% sample(resDataStart$protein_Id,20))
testData2954 <- list(resDataStart = smallData, config = config)
usethis::use_data(testData2954, overwrite = TRUE)

resultsV12954 <- LFQService::workflow_MQ_protoV1(resDataStart,
                                                 config,
                                                 path = NULL ,
                                                 peptideFilterFunction = LFQService:::.filter_proteins_by_peptide_count )

usethis::use_data(resultsV12954, overwrite = TRUE)


config <- testDataStart2954$config

patchOldConfig <- function(config){
  atable <- AnalysisTableAnnotation$new()
  atable$fileName = config$table$fileName
  # measurement levels.
  atable$hierarchy <- config$table$hierarchy
  #
  atable$ident_qValue = config$table$ident_qValue
  atable$workIntensity <- config$table$workIntensity
  atable$isotopeLabel = config$table$isotopeLabel
  atable$factors <- config$table$factors
  anaparam <- AnalysisParameters$new()
  configuration <- AnalysisConfiguration$new(atable, anaparam)
  return(configuration)
}

testDataStart2954$config <- patchOldConfig(testDataStart2954$config)
usethis::use_data(testDataStart2954, overwrite = TRUE)
