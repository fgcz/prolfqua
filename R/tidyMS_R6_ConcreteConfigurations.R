#' Generate instances of AnalysisConfiguration
#'
#' Configuration examples for various signal processing software outputs.
#' @rdname concrete_AnalysisConfiguration
#' @family configuration
#' @name concrete_AnalysisConfiguration
#' @examples
#' # Skyline configuration
#' skylineconfig <- AnalysisConfiguration$new()
#' skylineconfig$fileName <- "Replicate.Name"
#' skylineconfig$hierarchy[["protein_Id"]] <- "Protein.Name"
#' skylineconfig$hierarchy[["peptide_Id"]] <- "Peptide.Sequence"
#' skylineconfig$hierarchy[["precursor_Id"]] <- c("Peptide.Sequence", "Precursor.Charge")
#' skylineconfig$hierarchy[["fragment_Id"]] <- c("Peptide.Sequence", "Precursor.Charge", "Fragment.Ion", "Product.Charge")
#' skylineconfig$ident_qValue <- "annotation_QValue"
#' skylineconfig$set_response("Area")
#' skylineconfig$isotopeLabel <- "Isotope.Label"
#' skylineconfig$factors[["Time"]] = "Sampling.Time.Point"
#' skylineconfig$factor_keys()
#' skylineconfig$hierarchy_keys()
#'
#' # Spectronaut configuration
#' spectronautconfig <- AnalysisConfiguration$new()
#' spectronautconfig$fileName <- "R.FileName"
#' spectronautconfig$hierarchy[["protein_Id"]] <- "PG.ProteinAccessions"
#' spectronautconfig$hierarchy[["peptide_Id"]] <- "PEP.StrippedSequence"
#' spectronautconfig$hierarchy[["modPeptide_Id"]] <- "EG.ModifiedSequence"
#' spectronautconfig$hierarchy[["precursor_Id"]] <- c("EG.ModifiedSequence", "FG.Charge")
#' spectronautconfig$ident_qValue <- "EG.Qvalue"
#' spectronautconfig$workIntensity <- "FG.Quantity"
#' spectronautconfig$isotopeLabel <- "Isotope.Label"
#' spectronautconfig$factors[["coding"]] = "coding"
#' spectronautconfig$factors[["sex"]] = "sex"
#' spectronautconfig$factors[["age"]] = "age"
#' spectronautconfig$factors[["Sample_id"]] = "Sample.Name"
NULL

#' MaxQuant peptide file configuration
#'
#' file must be read with tidyMQ_Peptides, you will still need to add the
#' factors (explanatory variables).
#'
#' @param ident_qValue pep
#' @param intensity peptide.intensity
#' @param isotopeLabel isotope
#' @rdname concrete_AnalysisConfiguration
#' @export
#' @examples
#' tmp <- create_config_MQ_peptide()
#'
create_config_MQ_peptide <- function(ident_qValue = "pep", intensity = "peptide.intensity", isotopeLabel = "isotope") {
  config <- AnalysisConfiguration$new()
  config$fileName <- "raw.file"
  # measurement levels.
  config$hierarchy[["protein_Id"]] <- c("leading.razor.protein")
  config$hierarchy[["peptide_Id"]] <- c("sequence")
  config$hierarchyDepth <- 1

  config$ident_qValue <- ident_qValue
  config$set_response(intensity)
  config$isotopeLabel <- isotopeLabel
  config$min_peptides_protein <- 2

  return(config)
}
