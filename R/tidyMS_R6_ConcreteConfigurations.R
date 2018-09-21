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
createSpectronautPeptideConfiguration <- function(isotopeLabel="Isotope.Label", ident_qValue="EG.Qvalue"){
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
