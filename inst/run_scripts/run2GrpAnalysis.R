library(prolfqua)
library(tidyverse)


datadir <- file.path(find.package("prolfquaData") , "quantdata")
inputMQfile <-  file.path(datadir, "MAXQuant_ComboCourse_p2691_March_2018_WU183012.zip")
inputAnnotation <- file.path(datadir, "annotation_ComboCourse_p2691_March_2018_WU183012.xlsx")

startdata <- prolfqua::tidyMQ_ProteinGroups(inputMQfile)
ncol(startdata)
str(startdata)
annotation <- readxl::read_xlsx(inputAnnotation)
annotation$experiment = "p2691"
#View(annotation)
ncol(annotation)


startdata <- inner_join(annotation, startdata, by = "raw.file")
startdata

startdata <- filter(startdata, nr.peptides > 1)

startdata <- startdata %>% mutate(proteinAnnot = case_when(grepl("^REV_",majority.protein.ids) ~ "REV",
                                              grepl("^zz|^CON",majority.protein.ids) ~ "CON",
                                              TRUE ~ "FW"))

distinctprotid <- startdata %>% select(majority.protein.ids, proteinAnnot) %>% distinct()
table(distinctprotid$proteinAnnot)
relevantProteinHeaders <- c("proteinAnnot","fasta.header")


#### Create configuration ####

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
str(adata)


###
lfqdata <- LFQData$new(adata, config)
lfqdata$remove_small_intensities()
lfqdata$factors()


### Do some type of data normalization
lt <- lfqdata$get_Transformer()
transformed <- lt$log2_robscale()
transformed$config$table$is_intensity_transformed


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

Contrasts <- c("Glucose - Ethanol" = "Condition_Glucose - Condition_Ethanol")
Contrasts <- c("GvsE" = "Condition_Glucose - Condition_Ethanol")

m <- mod$modelDF$linear_model[[1]]
linfct <- linfct_from_model(m)
prolfqua::linfct_all_possible_contrasts(linfct$linfct_factors)
prolfqua::linfct_all_possible_contrasts(linfct$linfct_interactions)

contr <- prolfqua::Contrasts$new(mod, Contrasts)
ContrastsModerated$undebug("get_contrasts")
conrM <- ContrastsModerated$new(contr)
conrM$get_contrasts()

ContrastsSimpleImpute$undebug("get_contrasts")
mC <- ContrastsSimpleImpute$new(lfqdata = transformed, contrasts = Contrasts)
tmp <- mC$get_contrasts()

conrMI <- ContrastsModerated$new(mC)
modtmp <- conrMI$get_contrasts()

plot(tmp$p.value, modtmp$p.value)

### -----



LFQdata <- prolfqua::data_Yeast2Factor
LFQdata$data
var <- summarize_stats(LFQdata$data, LFQdata$config, all = TRUE)
pooled <- var %>% filter(!!sym(LFQdata$config$table$fkeysDepth()[1]) == "pooled")



