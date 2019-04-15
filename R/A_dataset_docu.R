#' List with peptide intensities
#'
#' @format A list of data frames
"correlatedPeptideList"

# skylinePRMSampleData ----
#' Data frame which can be transformed into \link{sample_analysis}
#' by applying the \link{skylineconfig} using the function \link{setup_analysis}.
#'
#' @examples
#' data(skylineconfig)
#' data(skylinePRMSampleData)
#' x<-setup_analysis(skylinePRMSampleData,skylineconfig)
#' all.equal(x,sample_analysis)
#'
#' @format A AnalysisConfiguration R6 class
"skylinePRMSampleData"

#' A configuration which matches the \link{sample_analysis} data.
#'
#'
#' @format A AnalysisConfiguration R6 class
#' @examples
#' skylineconfig <- LFQService::createSkylineConfiguration(isotopeLabel="Isotope.Label.Type",
#'  ident_qValue="Detection.Q.Value")
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#' #skylineconfig$table$factors
#' # usethis::use_data( skylineconfig , overwrite = TRUE )
#' tt <- R6extractValues(skylineconfig)
#' #yaml::write_yaml(tt,file=file.path("skylineconfig.yml"))
#' skylineconfig$table$hkeysLevel()
"skylineconfig"

#' A data frame wich goes along with the \link{skylineconfig}
#' and was generate from \link{skylinePRMSampleData}
#'
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{protein_Id}{protein id}
#'   \item{peptide_Id}{peptide id - stripped sequence}
#'   ...
#' }
#' @source \url{http://www.fgcz.ch/}
#' @examples
#'
#' skylineconfig <- createSkylineConfiguration(isotopeLabel="Isotope.Label.Type",
#'  ident_qValue="Detection.Q.Value")
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#' data(skylinePRMSampleData)
#' sample_analysis <- setup_analysis(skylinePRMSampleData, skylineconfig)
#' sample_analysis <- remove_small_intensities(sample_analysis, skylineconfig, threshold = 30)
#' #usethis::use_data( sample_analysis , overwrite = TRUE )
"sample_analysis"

# skylineSRM_HL_data ----
#' Data frame with a skyline export for an heavy to light experiment
#'
#'
#' @format a data frame
"skylineSRM_HL_data"

#' A configuration which matches the \link{sample_analysis} data.
#'
#'
#' @format A AnalysisConfiguration R6 class
#' @examples
#' skylineconfig_HL <- createSkylineConfiguration(isotopeLabel="Isotope.Label",
#'  ident_qValue="annotation_QValue")
#' skylineconfig_HL$table$factors[["treatment_c"]] <- "Condition2"
#' skylineconfig_HL$table$factors[["time_c"]] <- "time"
#' skylineconfig_HL$parameter$is_intensity_transformed = FALSE
#' #usethis::use_data( skylineconfig_HL , overwrite = TRUE )
"skylineconfig_HL"

#' A data frame wich goes along with the \link{skylineconfig_HL}.
#'
#'
#' @source \url{http://www.fgcz.ch/}
#' @examples
#'
#' skylineconfig_HL <- createSkylineConfiguration(isotopeLabel="Isotope.Label",
#'  ident_qValue="annotation_QValue")
#' skylineconfig_HL$table$factors[["treatment_c"]] <- "Condition2"
#' skylineconfig_HL$table$factors[["time_c"]] <- "time"
#' skylineconfig_HL$parameter$is_intensity_transformed = FALSE
#' data(skylineSRM_HL_data)
#' skylineSRM_HL_data$Area[skylineSRM_HL_data$Area == 0] <- NA
#' sample_analysis_HL <- setup_analysis(skylineSRM_HL_data, skylineconfig_HL)
#' #usethis::use_data( sample_analysis_HL , overwrite = TRUE )
"sample_analysis_HL"

#spectronautDIAData250----

#' Data frame of a Spectronaut export with DIA data
#'
#' @format a data frame
"spectronautDIAData250"

#' A configuration which matches the \link{spectronautDIAData250} and \link{spectronautDIAData250} data.
#'
#'
#' @format A AnalysisConfiguration R6 class
#' @examples
#' spectronautDIAData250_config <- createSpectronautPeptideConfiguration(isotopeLabel="Isotope.Label",
#'  ident_qValue="EG.Qvalue")
#' spectronautDIAData250_config$table$factors[["coding"]] = "coding"
#' spectronautDIAData250_config$table$factors[["sex"]] = "sex"
#' spectronautDIAData250_config$table$factors[["age"]] = "age"
#' spectronautDIAData250_config$table$factors[["Sample_id"]] = "Sample.Name"
#' #usethis::use_data( spectronautDIAData250_config , overwrite = TRUE )
#'
"spectronautDIAData250_config"

#' A data frame wich goes along with the \link{spectronautDIAData250_config}.
#'
#'
#' @source \url{http://www.fgcz.ch/}
#' @examples
#' data <- LFQService::spectronautDIAData250
#' xx <-  LFQService::spectronautDIAData250_config
#'
#' spectronautDIAData250_analysis <- setup_analysis(data,spectronautDIAData250_config)
#' #usethis::use_data( spectronautDIAData250_analysis , overwrite = TRUE )
"spectronautDIAData250_analysis"
