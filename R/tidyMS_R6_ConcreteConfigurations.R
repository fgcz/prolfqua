#' Generate instances of AnalysisConfiguration
#'
#' configurations examples of or various signal processing software outputs
#' @rdname concrete_AnalysisConfiguration
#' @family configuration
#' @name concrete_AnalysisConfiguration
NULL


#' Create a Skyline configuration
#'
#'
#' @param isotopeLabel Isotope.Label
#' @param ident_qValue annotation_QValue
#' @rdname concrete_AnalysisConfiguration
#' @export
#' @examples
#' skylineconfig <- create_config_Skyline()
#' skylineconfig$factors[["Time"]] = "Sampling.Time.Point"
#' skylineconfig$factor_keys()
#' skylineconfig$hierarchy_keys()
#'
#'
#'
create_config_Skyline <- function(isotopeLabel = "Isotope.Label", ident_qValue = "annotation_QValue") {
  config <- AnalysisConfiguration$new()
  config$fileName = "Replicate.Name"

  # measurement levels.
  config$hierarchy[["protein_Id"]] <- "Protein.Name"
  config$hierarchy[["peptide_Id"]] <- "Peptide.Sequence"
  config$hierarchy[["precursor_Id"]] <- c("Peptide.Sequence", "Precursor.Charge")
  config$hierarchy[["fragment_Id"]] <- c("Peptide.Sequence", "Precursor.Charge", "Fragment.Ion", "Product.Charge")

  #
  config$ident_qValue = ident_qValue
  config$set_response("Area")
  config$isotopeLabel = isotopeLabel
  config
}

#' Create Spectronaut configuration
#' @param isotopeLabel Isotope.Label
#' @param ident_qValue EG.Qvalue
#' @export
#' @rdname concrete_AnalysisConfiguration
#' @examples
#' spectronautconfig <- create_config_Spectronaut_Peptide()
#' config <- create_config_Spectronaut_Peptide()
#' config$factors[["coding"]] = "coding"
#' config$factors[["sex"]] = "sex"
#' config$factors[["age"]] = "age"
#' config$factors[["Sample_id"]] = "Sample.Name"
#'
create_config_Spectronaut_Peptide <- function(isotopeLabel = "Isotope.Label", ident_qValue = "EG.Qvalue") {
  config <- AnalysisConfiguration$new()
  config$fileName = "R.FileName"

  # measurement levels.
  config$hierarchy[["protein_Id"]] <- "PG.ProteinAccessions"
  config$hierarchy[["peptide_Id"]] <- "PEP.StrippedSequence"
  config$hierarchy[["modPeptide_Id"]] <- "EG.ModifiedSequence"
  config$hierarchy[["precursor_Id"]] <- c("EG.ModifiedSequence", "FG.Charge")

  #
  config$ident_qValue = ident_qValue
  config$workIntensity = "FG.Quantity"
  config$isotopeLabel = isotopeLabel
  config
}

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
  config$fileName = "raw.file"
  # measurement levels.
  config$hierarchy[["protein_Id"]] <- c("leading.razor.protein")
  config$hierarchy[["peptide_Id"]] <- c("sequence")
  config$hierarchyDepth <- 1
  #
  config$ident_qValue = ident_qValue
  config$set_response(intensity)
  config$isotopeLabel = isotopeLabel
  config$min_peptides_protein <- 2

  return(config)
}

#' Create configuration for MSFragger output
#'
#' @rdname concrete_AnalysisConfiguration
#' @export
#' @examples
#'create_config_MSFragger_MSstats()
#'
create_config_MSFragger_MSstats <- function() {
  ## Tell LFQ Service what column is what.
  config <- AnalysisConfiguration$new()
  # measurement levels.
  config$hierarchy[["protein_Id"]] <- c("ProteinName")
  config$hierarchy[["peptide_Id"]] <- c("PeptideSequence", "PrecursorCharge")
  config$fileName = "Run"
  config$ident_qValue = "pep"
  config$set_response("Intensity")
  config$isotopeLabel = "IsotopeLabelType"
  config$min_peptides_protein <- 2
  return(config)
}
