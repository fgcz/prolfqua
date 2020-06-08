rm(list = ls())
library(LFQService)
library(tidyverse)


outpath <- "results_test"
projectstruct <- createProjectStructure(outpath = outpath)
projectstruct$create()

inputMQfile <-  "../samples/p2558_05748_modspec.zip"
inputMQfile <-  "../samples/p2558_o5748_peptides.zip"


inputAnnotation <- "../samples/p2558_05748_annotation.xlsx"

source("c:/Users/wewol/prog/LFQService/R/tidyMS_MQProteinGroups.R")

mqdata <- tidyMQ_Peptides_Config(inputMQfile)
annotation <- readxl::read_xlsx(inputAnnotation)

# creates default configuration
annotation <- annotation %>% dplyr::filter(annotation$SCI != "un")

mqdata$config$table$factors[["drug_"]] = "genotype"
mqdata$config$table$factors[["SCI_"]] = "SCI"
mqdata$config$table$factorDepth <- 2

mqdata$config$order_Id = "o5748"
mqdata$config$project_Id = "p2558"
mqdata$config$workunit_Id = "20191120_MQ_repack.zip"


# specify model definition
memodel <- "~ drug_ * SCI_  + (1|peptide_Id) + (1|sampleName)"
Contrasts <- c("8wk_vs_1wk" = "SCI_8wk - SCI_1wk",
               "t_vs_v" = "drug_t - drug_v",
               "t_vs_v_given_8wk" = "`drug_t:SCI_8wk` - `drug_v:SCI_8wk`",
               "t_vs_v_given_1wk" = "`drug_t:SCI_1wk` - `drug_v:SCI_1wk`",
               "interaction_construct_with_time" = "t_vs_v_given_8wk - t_vs_v_given_1wk"
)


model <- LFQService:::make_model(formula = memodel, contrasts = Contrasts)
annotated_data <- application_add_annotation( mqdata$data,
                                              annotation,
                                              fileName = mqdata$config$table$fileName )

mqdata$data <- setup_analysis(annotated_data, mqdata$config )

filteredData <- filter_proteins_by_peptide_count(mqdata$data,
                                                 mqdata$config)
normalizedData <- normalize_log2_robscale(filteredData$data,
                                          filteredData$config)



protintensity_fun <- medpolish_protein_quants( normalizedData$data,
                                               normalizedData$config )

protdata <- protintensity_fun("unnest")
xx <- protintensity_fun("plot")

if (FALSE) {
  pdf(file.path(projectstruct$qc_path,"protein_inferences.pdf"))
  lapply(xx$plot, print)
  dev.off()
}


#source("c:/Users/wewol/prog/LFQService/R/LFQData.R")

pep <- LFQData$new(filteredData$data,filteredData$config,is_pep = TRUE, prefix = "peptide_")
prot <- LFQData$new(protdata$data,protdata$config, is_pep = FALSE, prefix = "protein_")

protplotter <- prot$get_Plotter()


protplotter$write(path_qc = projectstruct$qc_path)

#dd <- protplotter$boxplots()
#protplotter$write_boxplots(path_qc = projectstruct$qc_path)


protwriter <- prot$get_Writer()
protwriter$write_long(projectstruct$qc_path)
protwriter$write_wide(projectstruct$qc_path)
protwriter$file_paths
pepwriter <- pep$get_Writer()
pepwriter$write_wide(projectstruct$qc_path)
pepwriter$file_paths


message("######################## fit mixed #######################")
memodel <- paste0(normalizedData$config$table$getWorkIntensity() , memodel)
modelFunction <- make_custom_model_lmer( memodel, model_name = "Model")
reportColumns <- c("p.value",
                   "p.value.adjusted")






#source("c:/Users/wolski/prog/LFQService/R/tidyMS_application.R")
if (TRUE) {
  undebug(application_run_modelling_V2)
  resXXmixmodel <- application_run_modelling_V2(
    data = normalizedData$data,
    config = normalizedData$config,
    modelFunction = modelFunction,
    contrasts = Contrasts,
    modelling_dir = projectstruct$modelling_path )

  resXXmixmodel(do = "write_modelling")

  resXXmixmodel(do = "write_contrasts")
}


relevantParameters <- list(outpath = projectstruct$outpath,
                           inputMQfile = inputMQfile,
                           modelling_dir = "modelling_results_peptide",
                           workunit_Id = mqdata$config$workunit_Id,
                           annotation = annotation,
                           reportColumns = reportColumns,
                           config = mqdata$config,
                           model = memodel,
                           Contrasts = Contrasts,
                           project_Id = mqdata$config$project_Id,
                           order_Id = mqdata$config$order_Id,
                           author = "Witold Wolski <wew@fgcz.ethz.ch>"
)


LFQService::copy_mixed_model_analysis_script()
rmarkdown::render("mixed_model_analysis_script_Report.Rmd",
                  params = list(pars = relevantParameters),
                  output_format = "html_document",
                  output_dir = outpath,
                  output_file = "index.html")



