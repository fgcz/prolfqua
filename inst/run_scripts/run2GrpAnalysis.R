library(prolfqua)
library(tidyverse)


datadir <- file.path(find.package("prolfquaData") , "quantdata")
#inputMQfile <-  file.path(datadir, "MAXQuant_ComboCourse_p2691_March_2018_WU183012.zip")
#inputAnnotation <- file.path(datadir, "annotation_ComboCourse_p2691_March_2018_WU183012.xlsx")

inputMQfile <-  file.path(datadir, "MAXQuant_ComboCourse_p2370_March_2017_WU183008.zip")
inputAnnotation <- file.path(datadir, "annotation_ComboCourse_p2370_March_2017_WU183008.xlsx")
startdata <- prolfqua::tidyMQ_ProteinGroups(inputMQfile)

annotation <- readxl::read_xlsx(inputAnnotation)
annotation$experiment = "p2691"

################### CREATE SOME annotations
grp2 <- list()
grp2$projectID <- 3000
grp2$projectName <- "bbbbbbbbbbbbbbbbbbbb"
grp2$workunitID <- "WU444XDD"
grp2$nrPeptides <- 2

##################################### ProteinID statistics

startdata <- inner_join(annotation, startdata, by = "raw.file")
startdata <- filter(startdata, nr.peptides >= grp2$nrPeptides)

startdata <- startdata %>% mutate(proteinAnnot = case_when(grepl("^REV_",majority.protein.ids) ~ "REV",
                                              grepl("^zz|^CON",majority.protein.ids) ~ "CON",
                                              TRUE ~ "FW"))
distinctprotid <- startdata %>% select(protein_Id = majority.protein.ids, fasta.headers, proteinAnnot) %>% distinct()


grp2$percentOfContaminants <-  table(distinctprotid$proteinAnnot)["CON"]/sum(table(distinctprotid$proteinAnnot)) * 100
grp2$percentOfFalsePositives <- table(distinctprotid$proteinAnnot)["REV"]/sum(table(distinctprotid$proteinAnnot)) * 100
grp2$totalNrOfProteins <- sum(table(distinctprotid$proteinAnnot))
grp2$NrOfProteinsNoDecoys <- table(distinctprotid$proteinAnnot)["FW"]


############################## Create configuration ####
atable <- AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("majority.protein.ids")

atable$hierarchyDepth <- 1
atable$setWorkIntensity("mq.protein.intensity")

anaparam <- AnalysisParameters$new()
config <- AnalysisConfiguration$new(atable, anaparam)


config$table$factors[["Condition_"]] = "condition"
config$table$factors[["run"]] = "Run_ID"
config$table$factorDepth <- 1

adata <- setup_analysis(startdata, config)


########################################################

lfqdata <- LFQData$new(adata, config)
lfqdata$remove_small_intensities()
lfqdata$factors()

### Do some type of data normalization (or do not)
lt <- lfqdata$get_Transformer()

transformed <- lt$log2_robscale()


grp2$lfqData <- lfqdata
grp2$transformedlfqData <- transformed


# Do modelling
transformed$config$table$getWorkIntensity()
transformed$config$table$fkeysDepth()

formula_Condition <-  strategy_lm(paste0(transformed$config$table$getWorkIntensity(), " ~ ", transformed$config$table$fkeysDepth()))

# specify model definition
modelName  <- "Model"
unique(transformed$data$Condition_)

mod <- prolfqua::build_model(
  transformed$data,
  formula_Condition,
  subject_Id = transformed$config$table$hierarchyKeys() )

grp2$models <- mod

Contrasts <- c("Glucose - Ethanol" = "Condition_Glucose - Condition_Ethanol")
Contrasts <- c("GvsE" = "Condition_Glucose - Condition_Ethanol")


mxdf <- max(mod$modelDF$df.residual, na.rm =  TRUE)
m <- mod$modelDF %>% dplyr::filter(.data$df.residual ==  mxdf & isSingular == FALSE )
linfct <- linfct_from_model(m$linear_model[[1]])
prolfqua::linfct_all_possible_contrasts(linfct$linfct_factors)
xx <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_interactions)
Contrasts <- rownames(xx)
names(Contrasts) <- gsub(" - ", "_vs_",(rownames(xx)))

contr <- prolfqua::Contrasts$new(mod, Contrasts)
conrM <- ContrastsModerated$new(contr, modelName = "ModeratedLinear")
modLin <- conrM$get_contrasts()



mC <- ContrastsSimpleImpute$new(lfqdata = transformed, contrasts = Contrasts)
conrMI <- ContrastsModerated$new(mC, modelName = "ModeratedImpute")
modtmp <- conrMI$get_contrasts()


more <- setdiff(modtmp$protein_Id , modLin$protein_Id)
modtmpdistinct <- modtmp %>% filter(protein_Id %in% more)

grp2$ContrastsLM <- modLin
grp2$ContrastsImputed <- modtmpdistinct


pl <- prolfqua::Contrasts_Plotter$new(modtmpdistinct,subject_Id = "protein_Id")
pl$volcano()
pl$volcano_plotly()
pl$ma_plot()
pl$histogram()

# Plot proteins without p-values
xx <- modtmpdistinct[rowSums(is.na(modtmpdistinct)) > 0,]
xx <- xx %>% arrange(estimate)
ggplot2::ggplot(xx ,aes(x = reorder(protein_Id, estimate), y = estimate)) + ggplot2::geom_bar(stat="identity") + coord_flip()


### -----

rmarkdown::render("_Grp2Analysis.Rmd", params = list(grp = grp2), output_format = bookdown::html_document2())

