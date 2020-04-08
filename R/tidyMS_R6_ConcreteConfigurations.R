#' This function sets up an example configuration
#' @export
#' @examples
#' skylineconfig <- createSkylineConfiguration()
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#' skylineconfig$table$factorKeys()
#' skylineconfig$table$hierarchyKeys()
createSkylineConfiguration <- function(isotopeLabel="Isotope.Label", ident_qValue="annotation_QValue"){
  atable <- AnalysisTableAnnotation$new()
  atable$fileName = "Replicate.Name"

  # measurement levels.
  atable$hierarchy[["protein_Id"]] <- "Protein.Name"
  atable$hierarchy[["peptide_Id"]] <- "Peptide.Sequence"
  atable$hierarchy[["precursor_Id"]] <-  c("Peptide.Sequence","Precursor.Charge")
  atable$hierarchy[["fragment_Id"]] <- c("Peptide.Sequence","Precursor.Charge","Fragment.Ion", "Product.Charge")

  #
  atable$ident_qValue = ident_qValue
  atable$setWorkIntensity("Area")
  atable$isotopeLabel = isotopeLabel
  anaparam <- AnalysisParameters$new()
  configuration <- AnalysisConfiguration$new(atable, anaparam)
}

#' This function sets up an example spectronaut configuration
#' @export
#' @examples
#' spectronautconfig <- createSpectronautPeptideConfiguration()
#' config <- createSpectronautPeptideConfiguration()
#' config$table$factors[["coding"]] = "coding"
#' config$table$factors[["sex"]] = "sex"
#' config$table$factors[["age"]] = "age"
#' config$table$factors[["Sample_id"]] = "Sample.Name"
createSpectronautPeptideConfiguration <- function(isotopeLabel="Isotope.Label",
                                                  ident_qValue="EG.Qvalue"){
  atable <- AnalysisTableAnnotation$new()
  atable$fileName = "R.FileName"

  # measurement levels.
  atable$hierarchy[["protein_Id"]]    <-  "PG.ProteinAccessions"
  atable$hierarchy[["peptide_Id"]]    <-  "PEP.StrippedSequence"
  atable$hierarchy[["modPeptide_Id"]] <-  "EG.ModifiedSequence"
  atable$hierarchy[["precursor_Id"]]   <-  c("EG.ModifiedSequence", "FG.Charge")

  #
  atable$ident_qValue = ident_qValue
  atable$workIntensity = "FG.Quantity"
  atable$isotopeLabel = isotopeLabel
  anaparam <- AnalysisParameters$new()
  configuration <- AnalysisConfiguration$new(atable, anaparam)
}

#' MQ peptide file configuration - file most be read with tidyMQ_Peptides or tidyMQ_modificationSpecificPeptides
#' @export
create_MQ_peptide_Configuration <- function(ident_qValue = "pep",
                                  intensity = "peptide.intensity",
                                  isotopeLabel = "isotope"){
  atable <- AnalysisTableAnnotation$new()
  atable$fileName = "raw.file"
  # measurement levels.
  atable$hierarchy[["protein_Id"]] <- c("top_protein","protein.group.id")
  atable$hierarchy[["peptide_Id"]] <- c("sequence","peptide.id")
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

#' Create configuration for MSFragger output
#' @export
create_MSFragger_MSstats_config <- function(min_peptides_protein = 2){
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
  return(config)
}

