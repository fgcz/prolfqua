rm(list = ls())

library(tidyverse)
library(readr)
library(LFQService)

resPepProtAnnot <- read_csv(file = LFQService:::.find.package.file("LFQService",  "samples/testdata/annotatedData.csv"))

createMQSTYConfiguration <- function(ident_qValue = "pep",
                                     intensity = "intensity",
                                     isotopeLabel = "isotope"){
  atable <- AnalysisTableAnnotation$new()
  atable$fileName = "raw.file"
  # measurement levels.
  atable$hierarchy[["protein_Id"]] <- c("protein","protein.group.ids")
  #atable$hierarchy[["peptide_Id"]] <- c("sequence","peptide.id")
  atable$hierarchy[["peptide_Id"]] <- c("peptide.sequence.prob","charge","multiplicity","site.id")
  atable$hierarchyDepth <- 2
  #
  atable$ident_qValue = ident_qValue
  atable$setWorkIntensity(intensity)
  atable$isotopeLabel = isotopeLabel
  atable$factors[["NASH_"]] = "NASH_"
  atable$factors[["Steatosis_"]] = "Steatosis_"
  atable$factors[["sample_name"]] = "sample_name"
  atable$factorDepth <- 2
  anaparam <- AnalysisParameters$new()
  configuration <- AnalysisConfiguration$new(atable, anaparam)

  return(configuration)
}

config <- createMQSTYConfiguration()

config$parameter$min_peptides_protein <- 1
resDataStart <- setup_analysis(resPepProtAnnot, config)
head(resDataStart)

#resDataStart %>% filter(sample_name != "bbh_1688") -> resDataStart

resDataStart <- remove_small_intensities( resDataStart, config ) %>%
  complete_cases(config)


configuration <- config
config$table$factorDepth <- 2

#rmarkdown::render("MQSummary.Rmd",params = list(configuration = config$clone(deep=TRUE), data = resDataStart))

bb <- interaction_missing_stats(resDataStart, config, factors = NULL)
bb$data %>% dplyr::filter(nrNAs < nrReplicates -4) -> bb
resDataStart <- inner_join(bb,resDataStart)

#rmarkdown::render("MQSummary2.Rmd", params=list(configuration = config$clone(deep=TRUE), data = resDataStart), output_format = bookdown::pdf_document2())

config$table$workIntensity
#LFQService::render_MQSummary_rmd(resDataStart, config,dest_path = qc_path, workdir = ".")
filteredPep <- LFQService::transform_work_intensity(resDataStart, config, log2)
#
impfac_int <- missigness_impute_factors_interactions(filteredPep, config, value = "imputed")



Contrasts <- c("NASH_nafld-NASH_nash" = "NASH_nafld - NASH_nash",
               "Steatosis_s1-Steatosis_s2" = "Steatosis_s1 - Steatosis_s2" ,
               "Steatosis_s1-Steatosis_s3" = "Steatosis_s1 - Steatosis_s3",
               "Steatosis_s2-Steatosis_s3" = "Steatosis_s1 - Steatosis_s3",
               "NASH_nafld: Steatosis_s1 - Steatosis_s2" = "`NASH_nafld:Steatosis_s1` - `NASH_nafld:Steatosis_s2`",
               "NASH_nash: Steatosis_s1 - Steatosis_s2" = "`NASH_nash:Steatosis_s1` - `NASH_nash:Steatosis_s2`",
               "Interaction NASH:Steatosis 1 - 2" = "`NASH_nafld: Steatosis_s1 - Steatosis_s2` - `NASH_nash: Steatosis_s1 - Steatosis_s2`")

tmp <- missigness_impute_contrasts(impfac_int, config, Contrasts)


