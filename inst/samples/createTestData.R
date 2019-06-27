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
config$table$factorLevel <- 2
resPepProtAnnot %>% dplyr::filter(reverse == FALSE) -> resPepProtAnnot


resDataStart <- setup_analysis(resPepProtAnnot, config)
resDataStart <- remove_small_intensities(resDataStart, config) %>% completeCases(config)

smallData <- resDataStart %>% dplyr::filter(protein_Id %in% sample(resDataStart$protein_Id,20))
testData2954 <- list(resDataStart = smallData, config = config)
usethis::use_data(testData2954)
