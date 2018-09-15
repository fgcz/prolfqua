params <-
list(configuration = SRMService::skylineconfig_HL, data = SRMService::sample_analysis_HL)

## ----setup, include=FALSE------------------------------------------------
library(tidyverse)
library(SRMService)

knitr::opts_chunk$set(echo = FALSE, message=FALSE)
data <- eval(params$data)
configuration <- eval(params$configuration)

## ----hierarchyCounts-----------------------------------------------------

x <- hierarchyCounts( data , configuration )
knitr::kable(data.frame(NR = x), caption="summary")


## ----proteinSummaries----------------------------------------------------

x2 <- summarizeProteins(data, configuration)
x3 <- summarizeHierarchy(data, configuration)
tmp <- inner_join(x3, x2)
tmp <- tmp %>% select(-c(starts_with("min"), starts_with("max")))
knitr::kable(tmp, caption="Protein Summaries")


## ----missingFigIntensityHistorgram, fig.cap="histogram of mean condition intensities per transition"----
missignessHistogram(data, configuration)


## ----missingFigBarplot, fig.cap = "nr of missing"------------------------
xx <- missingPerCondition(data, configuration)
xx$figure

## ----missingFigBarplotCumsum, fig.cap="cumulative sum of missing"--------
xx <- missingPerConditionCumsum(data, configuration)
xx$figure

