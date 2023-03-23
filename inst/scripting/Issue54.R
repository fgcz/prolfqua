inputMQfile <-  file.path("path_to/proteinGroups.txt")
inputAnnotation <- file.path("path_to/annotation.xlsx")

startdata <- prolfqua::tidyMQ_ProteinGroups(inputMQfile)
annotation <- readxl::read_xlsx(inputAnnotation)

annotation$raw.file[1]
unique(startdata$raw.file)[1]

annotation$raw.file <- make.names(tolower(annotation$raw.file))
startdata <- dplyr::inner_join(annotation, startdata, by = "raw.file")


startdata <- dplyr::filter(startdata, nr.peptides > 1)
atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("proteinID")
atable$factors[["Condition"]] <- "Condition"
atable$factors[["PatientID"]] <- "sample"
atable$set_response("mq.protein.intensity")

config <- prolfqua::AnalysisConfiguration$new(atable)

adata <- prolfqua::setup_analysis(startdata, config) #the error comes here


