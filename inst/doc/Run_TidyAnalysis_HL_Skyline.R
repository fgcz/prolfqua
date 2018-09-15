## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
echo = FALSE,
message = FALSE,
warning = FALSE,
fig.width = 10,
fig.height = 10
)

## ------------------------------------------------------------------------
rm(list=ls())
library(conflicted)
library(tidyverse)
library(SRMService)

## ------------------------------------------------------------------------
outdir <- tempdir()

## ------------------------------------------------------------------------
allDataM <- SRMService::skylineSRM_HL_data
head(allDataM)


skylineconfig <- createSkylineConfiguration(isotopeLabel="Isotope.Label", ident_qValue="annotation_QValue")
skylineconfig$table$factors[["treatment_c"]] <- "Condition2"
skylineconfig$table$factors[["time_c"]] <- "time"
skylineconfig$parameter$is_intensity_transformed = FALSE

resData <- setup_analysis(allDataM, skylineconfig)

resData$Area[resData$Area == 0] <- NA

tt <- R6extractValues(skylineconfig)

## ----eval=FALSE----------------------------------------------------------
#  yaml::write_yaml(tt, file=file.path(outdir,"testSkyline.yml"))

## ----fig.cap="transitions in one plot."----------------------------------
proteinIDsymbol <- sym(skylineconfig$table$hierarchyKeys()[1])

xnested <- resData %>% group_by(UQ(proteinIDsymbol)) %>% nest()
linePlotHierarchy_configuration(xnested$data[[2]], xnested$protein_Id[[2]], skylineconfig)

## ----fig.cap="transitions in separate plot."-----------------------------
linePlotHierarchy_configuration(xnested$data[[3]], xnested$protein_Id[[3]], skylineconfig, separate = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  figs <- xnested %>%
#    mutate(plot = map2(data, UQ(proteinIDsymbol) , linePlotHierarchy_configuration, skylineconfig, separate = TRUE))
#  
#  pdf(file.path(outdir,"allProteins.pdf"), width = 10, height = 10)
#  invisible(lapply(figs$plot[1:3], print))
#  dev.off()

## ----eval = FALSE--------------------------------------------------------
#  rmarkdown::render("SRMReport.Rmd",
#                    params = list(data=resData, configuration = skylineconfig),
#                    envir = new.env(),
#                    output_file = "SRMReport.pdf")

## ------------------------------------------------------------------------
HLData <- spreadValueVarsIsotopeLabel(resData, skylineconfig)
HLData <- HLData %>% mutate(log2L_log2H = log2(light_Area)- log2(heavy_Area))
HLData$Isotope.Label <- "L/H"

skylineconfigHL <- skylineconfig
skylineconfigHL$table$workIntensity = "log2L_log2H"
skylineconfigHL$parameter$is_intensity_transformed = TRUE
skylineconfigHL$table$isotopeLabel

xnested <- HLData %>% group_by(UQ(proteinIDsymbol)) %>% nest()
skylineconfigHL$table$getWorkIntensity()

xnested$data[[2]][["log2L_log2H"]]

linePlotHierarchy_configuration(xnested$data[[2]], xnested$protein_Id[[2]], skylineconfigHL)
HLfigs <- xnested %>% mutate(plot = map2(data, UQ(proteinIDsymbol) , linePlotHierarchy_configuration, skylineconfig))


## ----eval=FALSE----------------------------------------------------------
#  
#  pdf(file.path(outdir,"allHLProteins.pdf"), width = 10, height = 10)
#  invisible(lapply(HLfigs$plot[1:3], print))
#  dev.off()

## ------------------------------------------------------------------------

HLfigs3 <- SRMService::applyToHierarchyBySample(HLData, skylineconfigHL, medpolishPly)
HLfigs3 <- inner_join(HLfigs3,HLfigs, by=skylineconfigHL$table$hierarchyKeys()[1])

linePlotHierarchy_configuration(xnested$data[[2]], xnested$protein_Id[[2]], skylineconfigHL) %>%
  linePlotHierarchy_QuantLine(HLfigs3$medpolishPly[[2]],"medpolish", skylineconfigHL)

HLfigs3 <- HLfigs3 %>% mutate(figsMed = map2(plot, medpolishPly, linePlotHierarchy_QuantLine, "medpolish", skylineconfig))


## ----eval=FALSE----------------------------------------------------------
#  pdf(file.path(outdir,"allProteinsWithMed.pdf"), width = 10, height = 10)
#  invisible(lapply(HLfigs3$figsMed[2:5], print))
#  dev.off()

## ----fig.cap="Plot protein intensities"----------------------------------
prots <- HLfigs3 %>% dplyr::select(as.character(proteinIDsymbol),  medpolishPly) %>% unnest()
write_tsv(prots, path = file.path(outdir,"proteins.tsv"))

ggplot(prots, aes(x=time_c, y=medpolish, group=treatment_c, color=treatment_c )) +
  geom_point() +
  geom_line() +
  facet_wrap(~protein_Id) +
  theme_classic()


## ----fig.cap="Plot protein intensities"----------------------------------
ggplot(prots, aes(x=time_c, y=medpolish, group=treatment_c, color=treatment_c )) +
  geom_point() +
  geom_line() +
  facet_wrap(~protein_Id, scales="free_y" ) +
  theme_classic()

