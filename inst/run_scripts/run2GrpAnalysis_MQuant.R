library(prolfqua)
library(tidyverse)


################### CREATE SOME annotations
datadir <- file.path(find.package("prolfquadata") , "quantdata")

inputMQfile <-  file.path(datadir, "MAXQuant_ComboCourse_p2370_March_2017_WU183008.zip")
inputAnnotation <- file.path(datadir, "annotation_ComboCourse_p2370_March_2017_WU183008.xlsx")

GRP2 <- list()

factorDisplayName <- "Condition_"
factorAnnotationName <- "condition"
Contrasts <- c("GvsE" = "Condition_Glucose - Condition_Ethanol")

GRP2$projectID <- 3000
GRP2$projectName <- "bbbbbbbbbbbbbbbbbbbb"
GRP2$workunitID <- "PDrun"

GRP2$nrPeptides <- 2

GRP2$log2FCthreshold <- 1
GRP2$FDRthreshold <- 0.1


##### Read the data.

startdata <- prolfqua::tidyMQ_ProteinGroups(inputMQfile)
startdata$majProtID <- gsub(";.+","",startdata$majority.protein.ids)
annotation <- readxl::read_xlsx(inputAnnotation)


##################################### ProteinID statistics #######

startdata <- inner_join(annotation, startdata, by = "raw.file")
startdata <- filter(startdata, nr.peptides >= GRP2$nrPeptides)
startdata <- startdata %>% mutate(proteinAnnot = case_when(grepl("^REV_",majority.protein.ids) ~ "REV",
                                              grepl("^zz|^CON",majority.protein.ids) ~ "CON",
                                              TRUE ~ "FW"))

distinctprotid <- startdata %>% select(protein_Id = majProtID, fasta.headers, proteinAnnot) %>% distinct()
desc <- distinctprotid %>% select(-proteinAnnot)

GRP2$percentOfContaminants <-  table(distinctprotid$proteinAnnot)["CON"]/sum(table(distinctprotid$proteinAnnot)) * 100
GRP2$percentOfFalsePositives <- table(distinctprotid$proteinAnnot)["REV"]/sum(table(distinctprotid$proteinAnnot)) * 100
GRP2$totalNrOfProteins <- sum(table(distinctprotid$proteinAnnot))
GRP2$NrOfProteinsNoDecoys <- table(distinctprotid$proteinAnnot)["FW"]


############################## Create configuration For MQ ####
atable <- AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("majProtID")
atable$hierarchyDepth <- 1
atable$setWorkIntensity("mq.protein.intensity")
config <- AnalysisConfiguration$new(atable)


config$table$factors[[factorDisplayName]] = factorAnnotationName
config$table$factors[["run"]] = "Run_ID"
config$table$factorDepth <- 1

adata <- setup_analysis(startdata, config)
adata$run <- as.numeric(adata$run)

##################### Preprocess intensities ###################################

lfqdata <- LFQData$new(adata, config)
lfqdata$remove_small_intensities()
lfqdata$factors()

### Do some type of data normalization (or do not)
lt <- lfqdata$get_Transformer()
transformed <- lt$log2()$robscale()$lfq


GRP2$lfqData <- lfqdata
GRP2$transformedlfqData <- transformed

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

wr <- GRP2$lfqData$get_Writer()
tmp <- wr$get_wide()
tmp2 <- GRP2$transformedlfqData$get_Writer()$get_wide()
names(tmp2) <- paste0(names(tmp2), ".normalized")
res <- inner_join(desc, GRP2$contrResult, by = c("protein_Id" = "protein_Id"))

writexl::write_xlsx(c(tmp, tmp2,  contrasts = list(res)), path = "AnalysisResults.xlsx")

rmarkdown::render("_GRP2Analysis.Rmd", params = list(grp = GRP2), output_format = bookdown::html_document2(toc = TRUE, toc_float = TRUE))

