# Run mixed models on benchmark dataset.



rm(list = ls())

library(conflicted)
library(LFQService)
library(tidyverse)
library(dplyr)

conflict_prefer("filter", "dplyr")



datadir <- file.path(find.package("prolfquaData") , "quantdata")
inputMQfile <-  file.path(datadir, "MaxQuant_Ionstar2018_PXD003881.zip")
inputAnnotation <- file.path(datadir, "annotation_Ionstar2018_PXD003881.xlsx")

ps <- ProjectStructure$new(outpath = ".",
                     project_Id = 3000,
                     order_Id = "IonStar" ,
                     workunit_Id = "IonStar",
                     inputData = inputMQfile,
                     inputAnnotation = inputAnnotation)

mqdata <- tidyMQ_Peptides_Config(ps$inputData)

annotation <- readxl::read_xlsx(ps$inputAnnotation)

mqdata$config$table$factors[["dilution."]] = "sample"
#mqdata$config$table$factors[["fPairing."]] = "fakePair_Id"
mqdata$config$table$factors[["run_Id"]] = "run_ID"
mqdata$config$table$factorDepth <- 1


# specify model definition

res <- application_add_annotation(
  mqdata$data,
  annotation
)

mqdata$data <- setup_analysis(res, mqdata$config)


IonstarData <- R6::R6Class(
  "IonstarData",
  public = list(
    data = NULL,
    config = NULL,
    initialize = function(data, config){
      self$data = data
      self$config = config
    },
    Pep = function(){
      return(list(data = self$data, config = self$config))
    },
    filtered = function(){
      res <- LFQService::filter_proteins_by_peptide_count(self$data, self$config)
      return(list(data = res$data, config = self$config$clone(deep = TRUE)))
    },
    normalized = function(){
      tmp <- self$filtered()
      conf <- tmp$config$clone(deep = TRUE)
      res <- normalize_log2_robscale(tmp$data,
                                     conf)
      res <- list(data =  res$data,
                  config = res$conf$clone(deep = TRUE))
      return(res)
    },
    subset_normalized = function(){
      tmp <- self$filtered()
      conf <- tmp$config$clone(deep = TRUE)
      data <- transform_work_intensity(tmp$data, conf, log2)
      subset <-
        data %>%
        dplyr::filter(grepl("HUMAN", data$protein_Id))
      scaldata <-
        scale_with_subset(data, subset, conf)
      ndata <- scaldata %>%
        dplyr::rename(transformedIntensity = conf$table$getWorkIntensity())
      conf$table$popWorkIntensity()
      conf$table$setWorkIntensity("transformedIntensity")
      return(list(data = ndata, config = conf))
    }
  )
)

ionstar <- IonstarData$new(mqdata$data,
                           mqdata$config$clone(deep = TRUE))

ionstar$Pep()
ionstar$filtered()
xx <- ionstar$normalized()

xx$config$table$is_intensity_transformed


xx <- ionstar$subset_normalized()
xx$config$table$is_intensity_transformed

usethis::use_data(ionstar, overwrite = TRUE)



hist(xx$data$log2_peptide.intensity)
hist(xx$data$transformedIntensity)


protL <- medpolish_protein_quants(xx$data, xx$config)
#protL("plot")

dataIonstarProtein_subsetNorm <- list(
  data = protL("unnest")$data,
  config = protL("unnest")$config
)
usethis::use_data(dataIonstarProtein_subsetNorm, overwrite = TRUE)
