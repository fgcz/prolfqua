#' Cenerate instances of AnalysisTableAnnotation
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
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#' skylineconfig$table$factor_keys()
#' skylineconfig$table$hierarchyKeys()
create_config_Skyline <- function(isotopeLabel="Isotope.Label",
                                  ident_qValue="annotation_QValue"){
  atable <- AnalysisTableAnnotation$new()
  atable$fileName = "Replicate.Name"

  # measurement levels.
  atable$hierarchy[["protein_Id"]] <- "Protein.Name"
  atable$hierarchy[["peptide_Id"]] <- "Peptide.Sequence"
  atable$hierarchy[["precursor_Id"]] <-  c("Peptide.Sequence","Precursor.Charge")
  atable$hierarchy[["fragment_Id"]] <- c("Peptide.Sequence","Precursor.Charge","Fragment.Ion", "Product.Charge")

  #
  atable$ident_qValue = ident_qValue
  atable$set_response("Area")
  atable$isotopeLabel = isotopeLabel
  anaparam <- AnalysisParameters$new()
  AnalysisConfiguration$new(atable, anaparam)
}

#' Create Spectronaut configuration
#' @param isotopeLabel Isotope.Label
#' @param ident_qValue EG.Qvalue
#' @export
#' @rdname concrete_AnalysisConfiguration
#' @examples
#' spectronautconfig <- create_config_Spectronaut_Peptide()
#' config <- create_config_Spectronaut_Peptide()
#' config$table$factors[["coding"]] = "coding"
#' config$table$factors[["sex"]] = "sex"
#' config$table$factors[["age"]] = "age"
#' config$table$factors[["Sample_id"]] = "Sample.Name"
#'
create_config_Spectronaut_Peptide <- function(isotopeLabel="Isotope.Label",
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
  AnalysisConfiguration$new(atable, anaparam)
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
create_config_MQ_peptide <- function(ident_qValue = "pep",
                                  intensity = "peptide.intensity",
                                  isotopeLabel = "isotope"){
  atable <- AnalysisTableAnnotation$new()
  atable$fileName = "raw.file"
  # measurement levels.
  atable$hierarchy[["protein_Id"]] <- c("leading.razor.protein")
  atable$hierarchy[["peptide_Id"]] <- c("sequence")
  atable$hierarchyDepth <- 1
  #
  atable$ident_qValue = ident_qValue
  atable$set_response(intensity)
  atable$isotopeLabel = isotopeLabel

  anaparam <- AnalysisParameters$new()
  anaparam$min_peptides_protein <- 2
  configuration <- AnalysisConfiguration$new(atable, anaparam)

  return(configuration)
}

#' Create configuration for MSFragger output
#'
#' @rdname concrete_AnalysisConfiguration
#' @export
#' @examples
#'create_config_MSFragger_MSstats()
#'
create_config_MSFragger_MSstats <- function(){
  ## Tell LFQ Service what column is what.
  atable <- AnalysisTableAnnotation$new()
  # measurement levels.
  atable$hierarchy[["protein_Id"]] <- c("ProteinName")
  atable$hierarchy[["peptide_Id"]] <- c("PeptideSequence","PrecursorCharge")
  atable$fileName = "Run"
  atable$ident_qValue = "pep"
  atable$set_response("Intensity")
  atable$isotopeLabel = "IsotopeLabelType"
  anaparam <- AnalysisParameters$new()
  anaparam$min_peptides_protein <- 2
  config <- AnalysisConfiguration$new(atable, anaparam)
  return(config)
}

