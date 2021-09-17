rm(list = ls())
library(prolfqua)
library(tidyverse)

PG2a <- readxl::read_xlsx("data/Export_for_WW/o25914_Proteins..xlsx")
res <- pivot_longer(PG2a, starts_with("Abundance"),names_to = "colname", values_to = "Abundance" )
colnames(res) <- make.names(colnames(res))
resA <- res %>% separate(colname, c(NA, "F", "channel", NA,"Condition"))
resA <- rename(resA, nr.peptides = X..Peptides)



annotation <- readxl::read_xlsx("TMT_lableling_o25914.xlsx")
colnames(annotation) <- make.names(colnames(annotation))
annotation <- annotation %>% mutate(channel = gsub("^TMT","", TMT.Tag)) %>% select(channel,Sample )
annotation$channel[!annotation$channel %in% unique(resA$channel)]
annotation <- annotation %>% separate(Sample, c(NA, "ID"), remove = FALSE)

resAa <- inner_join(annotation, resA)
resAa <- resAa %>% filter(!ID %in% c("SL1105","SL1665") )

################### CREATE SOME annotations
GRP2 <- list()

factorDisplayName <- "Condition_"
factorAnnotationName <- "Condition"
Contrasts <- c("DiseasedvsControl" = "Condition_Diseased - Condition_Control")

GRP2$projectID <- "25914"
GRP2$projectName <- "Order_25914"
GRP2$workunitID <- "Manual analysis"

GRP2$nrPeptides <- 2
GRP2$log2FCthreshold <- 0.26
GRP2$FDRthreshold <- 0.1



##################################### ProteinID statistics #######

startdata <- resAa
startdata <- filter(startdata, nr.peptides >= GRP2$nrPeptides)
startdata <- startdata %>% mutate(proteinAnnot = "FW")


distinctprotid <- startdata %>% select(protein_Id = Accession, fasta.headers = Description) %>% distinct()

distinctprotid <- distinctprotid %>% mutate(proteinAnnot = case_when(grepl("^REV_",protein_Id) ~ "REV",
                                              grepl("^Y-FGCZ",protein_Id) ~ "CON",
                                              TRUE ~ "FW"))

distinctprotid$geneName <- sub(".* GN=(.+) PE.*", "\\1", distinctprotid$fasta.headers)
GRP2$percentOfContaminants <-  table(distinctprotid$proteinAnnot)["CON"]/sum(table(distinctprotid$proteinAnnot)) * 100
GRP2$percentOfFalsePositives <- table(distinctprotid$proteinAnnot)["REV"]/sum(table(distinctprotid$proteinAnnot)) * 100
GRP2$totalNrOfProteins <- sum(table(distinctprotid$proteinAnnot))
GRP2$NrOfProteinsNoDecoys <- table(distinctprotid$proteinAnnot)["FW"]


############################## Create configuration For MQ ####
atable <- AnalysisTableAnnotation$new()
atable$fileName = "Sample"
atable$hierarchy[["protein_Id"]] <- c("Accession")
atable$hierarchyDepth <- 1
atable$setWorkIntensity("Abundance")
config <- AnalysisConfiguration$new(atable)


config$table$factors[[factorDisplayName]] = factorAnnotationName
config$table$factors[["ID"]] = "ID"
config$table$factorDepth <- 1

startdata$cha <- startdata$channel
adata <- setup_analysis(startdata, config)

##################### Preprocess intensities ###################################

lfqdata <- LFQData$new(adata, config)
lfqdata$data$Abundance <- lfqdata$data$Abundance * 1e6
lfqdata$remove_small_intensities()

### Do some type of data normalization (or do not)
lt <- lfqdata$get_Transformer()
transformed <- lt$log2()$robscale()$lfq


GRP2$lfqData <- lfqdata
GRP2$transformedlfqData <- transformed
#GRP2$transformedlfqData$data$transformedIntensity <- GRP2$transformedlfqData$data$transformedIntensity * 3

################## Run Modelling ###############


formula_Condition <-  strategy_lm(paste0(transformed$config$table$getWorkIntensity(), " ~ ",
                                         transformed$config$table$fkeysDepth()))
# specify model definition
modelName  <- "Model"

mod <- prolfqua::build_model(
  transformed$data,
  formula_Condition,
  subject_Id = transformed$config$table$hierarchyKeys() )

GRP2$models <- mod


GRP2$Contrast <- Contrasts

contr <- prolfqua::Contrasts$new(mod, Contrasts)
conrM <- ContrastsModerated$new(contr, modelName = "Linear_Model_Moderated")
mC <- ContrastsSimpleImpute$new(lfqdata = transformed, contrasts = Contrasts)
conMI <- ContrastsModerated$new(mC, modelName = "Imputed_Data")

res <- prolfqua::addContrastResults(conrM, conMI)

GRP2$contrResult <- res$merged$get_contrasts()
GRP2$contrMerged <- res$merged$get_Plotter()
GRP2$contrMerged$fcthresh = GRP2$log2FCthreshold
GRP2$contrMerged$volcano_spec[[1]]$thresh = GRP2$FDRthreshold

GRP2$contrMore <- res$more$get_Plotter()

top20 <- GRP2$contrResult %>% dplyr::select( protein_Id,log2FC= estimate,conf.low,conf.high, FDR ) %>%
  arrange(FDR) %>%
  head(20)
GRP2$top20 <- top20
#knitr::kable(top20, caption = "Top 20 proteins sorted by smallest Q Value (adj.P.Val). The effectSize column is the log2 FC of condition vs reference.")

GRP2$top20confint <- ggplot(top20, aes(x = protein_Id, y = log2FC,
                    ymin = conf.low, ymax = conf.high)) +
  geom_hline( yintercept = 0, color = 'red' ) +
  geom_linerange() + geom_point() + coord_flip() + theme_minimal()


protMore <- GRP2$transformedlfqData$get_copy()
protMore$complete_cases()
protMore$data <- protMore$data %>% filter(.data$protein_Id %in% res$more$contrast_result$protein_Id)
GRP2$imputedProteins <- protMore

# Plot proteins without p-values

xx <- res$more$contrast_result[rowSums(is.na(res$more$contrast_result)) > 0,]
if (nrow(xx) > 0) {
  xx <- xx %>% arrange(estimate)
  GRP2$noPvalEstimate <- ggplot2::ggplot(xx ,aes(x = reorder(protein_Id, estimate), y = estimate)) +
    ggplot2::geom_bar(stat = "identity") + coord_flip()
  missing <- GRP2$transformedlfqData$get_copy()
  missing$complete_cases()
  missing$data <- missing$data %>% dplyr::filter(protein_Id %in% xx$protein_Id)
  missing$get_Plotter()$raster()
}


### -----

names(GRP2)
wr <- GRP2$lfqData$get_Writer()
tmp <- wr$get_wide()
tmp$data <- inner_join(distinctprotid, tmp$data)
names(tmp)
tmp2 <- GRP2$transformedlfqData$get_Writer()$get_wide()

names(tmp2) <- paste0(names(tmp2), ".normalized")
tmp2$data.normalized <- inner_join(distinctprotid, tmp2$data.normalized)
res <- inner_join(distinctprotid, GRP2$contrResult)
writexl::write_xlsx(c(tmp, tmp2,  contrasts = list(res)), path = "AnalysisResults.xlsx")
rmarkdown::render("_GRP2Analysis.Rmd",
                  params = list(grp = GRP2) ,
                  output_format = bookdown::html_document2(toc = TRUE,toc_float = TRUE))

