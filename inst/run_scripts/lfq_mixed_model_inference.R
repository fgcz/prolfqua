rm(list = ls())
library(LFQService)
library(tidyverse)


outpath <- "results_test"
ps <- ProjectStructure$new(outpath = outpath,
                                      order_Id = "o5748",
                                      project_Id = "p2558",
                                      workunit_Id = "20191120_MQ_repack.zip",
                                      inputData = "../samples/p2558_o5748_peptides.zip",
                                      inputAnnotation = "../samples/p2558_05748_annotation.xlsx")

ps$create()


mqdata <- tidyMQ_Peptides_Config(ps$inputData)
annotation <- readxl::read_xlsx(ps$inputAnnotation)

# creates default configuration
annotation <- annotation %>% dplyr::filter(annotation$SCI != "un")

mqdata$config$table$factors[["drug_"]] = "genotype"
mqdata$config$table$factors[["SCI_"]] = "SCI"
mqdata$config$table$factorDepth <- 2




# specify model definition
memodel <- "~ drug_ * SCI_  + (1|peptide_Id) + (1|sampleName)"
Contrasts <- c("8wk_vs_1wk" = "SCI_8wk - SCI_1wk",
               "t_vs_v" = "drug_t - drug_v",
               "t_vs_v_given_8wk" = "`drug_t:SCI_8wk` - `drug_v:SCI_8wk`",
               "t_vs_v_given_1wk" = "`drug_t:SCI_1wk` - `drug_v:SCI_1wk`",
               "interaction_construct_with_time" = "t_vs_v_given_8wk - t_vs_v_given_1wk"
)



annotated_data <- add_annotation( mqdata$data,
                                              annotation,
                                              fileName = mqdata$config$table$fileName )

mqdata$data <- setup_analysis(annotated_data, mqdata$config )
filteredData <- filter_proteins_by_peptide_count(mqdata$data,
                                                 mqdata$config)
normalizedData <- normalize_log2_robscale(filteredData$data,
                                          mqdata$config)
protintensity_fun <- medpolish_protein_quants( normalizedData$data,
                                               normalizedData$config )

protdata <- protintensity_fun("unnest")
xx <- protintensity_fun("plot")
prefix <- "protein_"
if (TRUE) {
  pdf(file.path(ps$qc_path(),paste0(prefix ,"inference_figures.pdf")))
  lapply(xx$plot, print)
  dev.off()
}


pep <- LFQData$new(filteredData$data,mqdata$config,is_pep = TRUE, prefix = "peptide_")
prot <- LFQData$new(protdata$data,protdata$config, is_pep = FALSE, prefix = "protein_")




pep$render(qc_path = ps$qc_path())
prot$render(qc_path = ps$qc_path())

protplotter <- prot$get_Plotter()
protplotter$write(path_qc = ps$qc_path())
protplotter$write_boxplots(path_qc = ps$qc_path())


protwriter <- prot$get_Writer()
protwriter$write_long(ps$qc_path())
protwriter$write_wide(ps$qc_path())
#LFQDataWriter$undebug("write_wide")
pepwriter <- pep$get_Writer()
pepwriter$write_wide(ps$qc_path())


message("######################## fit mixed #######################")
memodel <- paste0(normalizedData$config$table$getWorkIntensity() , memodel)
modelFunction <- make_custom_model_lmer( memodel, report_columns = c("p.value", "p.value.adjusted"))



if (TRUE) {
  debug(application_run_modelling_V2)
  resXXmixmodel <- application_run_modelling_V2(
    data = normalizedData$data,
    config = normalizedData$config,
    modelFunction = modelFunction,
    contrasts = Contrasts,
    modelling_dir = ps$modelling_path() )


  xx <- resXXmixmodel(DEBUG = TRUE)
  resXXmixmodel(do = "write_modelling")

  resXXmixmodel(do = "write_contrasts")
}

relevantParameters <- list(ps = ps,

                           prefix = prefix,
                           annotation = prot$factors(),

                           model = modelFunction,
                           Contrasts = Contrasts,


                           author = "Witold Wolski <wew@fgcz.ethz.ch>"
)


LFQService::copy_mixed_model_analysis_script()
#file.copy("../rmarkdown/mixed_model_analysis_script_Report.Rmd",".",overwrite = TRUE)
rmarkdown::render("mixed_model_analysis_script_Report.Rmd",
                  params = list(pars = relevantParameters),
                  output_format = "html_document",
                  output_dir = outpath,
                  output_file = "index.html")


# Cleanup code
if(grepl("inst/run_scripts", getwd())){
  files <- dir()
  if("lfq_mixed_model_inference.R" %in% files){
    toremeove <- files[!grepl("*.R$",files)]
    unlink(toremeove, recursive = TRUE)
  }
}

