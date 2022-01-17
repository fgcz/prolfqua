# Author Witold Wolski
# wew@fgcz.ethz.ch
library(prolfqua)
library(tidyverse)

outpath <- "results_SamplesRemoved_V2/"

PG2a <- readxl::read_xlsx("data/Export_for_WW/o25914_Proteins..xlsx")
res <- pivot_longer(PG2a, starts_with("Abundance"),names_to = "colname", values_to = "Abundance" )
colnames(res) <- make.names(colnames(res))
resA <- res |> separate(colname, c(NA, "F", "channel", NA,"Condition"))
resA <- rename(resA, nr.peptides = X..Peptides)


annotation <- readxl::read_xlsx("TMT_lableling_o25914.xlsx")
colnames(annotation) <- make.names(colnames(annotation))
annotation <- annotation |> mutate(channel = gsub("^TMT","", TMT.Tag)) |> select(channel,Sample )
annotation$channel[!annotation$channel %in% unique(resA$channel)]
annotation <- annotation |> separate(Sample, c(NA, "ID"), remove = FALSE)

resAa <- inner_join(annotation, resA)
resAa <- resAa |> filter(!ID %in% c("SL1105","SL1665") )



################### CREATE SOME annotations
GRP2 <- list()
GRP2$Contrasts <- c("DiseasedvsControl" = "Condition_Diseased - Condition_Control")


GRP2$projectID <- "25914"
GRP2$projectName <- "Order_25914"
GRP2$workunitID <- "Manual analysis"

GRP2$nrPeptides <- 2
resA <- filter(resAa, nr.peptides >= GRP2$nrPeptides)


GRP2$log2FCthreshold <- 0.26
GRP2$FDRthreshold <- 0.1

# setup configuration

atable <- AnalysisTableAnnotation$new()
atable$fileName = "Sample"
atable$hierarchy[["protein_Id"]] <- c("Accession")
atable$hierarchyDepth <- 1


atable$factors[["Condition_"]] = "Condition"
atable$factors[["ID"]] = "ID"
atable$factorDepth <- 1

atable$setWorkIntensity("Abundance")

prolfqua::make2grpReport(resAa, atable, GRP2,  Description = "Description", outpath)



##################################### ProteinID statistics #######


