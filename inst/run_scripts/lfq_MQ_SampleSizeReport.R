#!/usr/bin/Rscript

suppressMessages(library(readr))
suppressMessages(library(tidyverse))
suppressMessages(library(LFQService))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

"Sample Size Report from MQ file

Usage:
  lfq_MQ_SampleSizeReport.R <mqzip> [--outdir=<outdir>]

Options:
  -o --outdir=<outdir> output directory [default: samplesize_qc_report]

Arguments:
  mqzip input file
" -> doc

library(docopt)
opt <- docopt(doc)

cat("\nParameters used:\n\t",
    "     mqzip:", mqzip <- opt$mqzip, "\n\t",
    "result_dir:", outputDir <- opt[["--outdir"]], "\n\n\n")


if(!dir.exists(outputDir)){
  if(!dir.create(outputDir)){
    stop("\n could not create : ", outputDir, "\n")
  }
}

assign("lfq_write_format", "xlsx", envir = .GlobalEnv)


summarize_cv_raw_transformed <- function(resDataStart, config){
  config_tmp <- config$clone(deep=TRUE)
  stats_res <- summarize_cv(resDataStart, config_tmp)
  wide <- toWideConfig( resDataStart,config_tmp)
  stats_raw <- inner_join(stats_res, wide$data, by=config_tmp$table$hierarchyKeys())

  resDataStart<- transform_work_intensity(resDataStart, config_tmp, log2)
  data <- LFQService::applyToIntensityMatrix(resDataStart, config_tmp, .func = robust_scale)

  stats_res_transformed <- summarize_cv(data, config_tmp)
  wide <- toWideConfig( data,config_tmp)
  stats_transformed <- inner_join(stats_res_transformed, wide$data, by=config_tmp$table$hierarchyKeys())
  peptideStats <- inner_join(stats_raw, stats_transformed, by=c( config_tmp$table$factorKeys(),config_tmp$table$hierarchyKeys() ), suffix =c(".raw",".transformed") )
  return(peptideStats)
}



resPepProt <- tidyMQ_PeptideProtein(mqzip,.all = TRUE)
resPepProtAnnot <- tidyMQ_from_modSpecific_to_peptide(resPepProt$mq_modSpecPeptides, resPepProt$mq_peptides)
resPepProtAnnot <- tidyMQ_top_protein_name(resPepProtAnnot)

resPepProtAnnot$QC <- "QC"
resPepProtAnnot$isotope <- "light"

createMQProteinPeptideConfiguration <- function(ident_qValue = "pep",
                                                intensity = "peptide.intensity",
                                                isotopeLabel = "isotope"){
  atable <- AnalysisTableAnnotation$new()
  atable$fileName = "raw.file"
  # measurement levels.
  atable$hierarchy[["protein_Id"]] <- c("top_protein","protein.group.id")
  atable$hierarchy[["peptide_Id"]] <- c("sequence","peptide.id")
  atable$ident_qValue = ident_qValue
  atable$setWorkIntensity(intensity)
  atable$isotopeLabel = isotopeLabel
  atable$factors[["FACTOR"]] = "QC"
  atable$factorLevel <- 1
  anaparam <- AnalysisParameters$new()
  configuration <- AnalysisConfiguration$new(atable, anaparam)
  return(configuration)
}

config <- createMQProteinPeptideConfiguration()
resDataStart <- setup_analysis(resPepProtAnnot, config)
resDataStart <- remove_small_intensities(resDataStart, config) %>% completeCases(config)


filename <- basename(tools::file_path_sans_ext(mqzip))
LFQService::render_MQSummary_rmd(resDataStart,
                                 config$clone(deep=TRUE),
                                 pep = TRUE,
                                 dest_path = outputDir,
                                 dest_file_name = paste0("r_",filename,"_Peptide.pdf"),workdir = "." )


peptideStats <- summarize_cv_raw_transformed(resDataStart, config)
lfq_write_table(peptideStats, path=file.path(outputDir,paste0("r_",filename,"_PeptideStats.csv")))


# Perform filtering and Protein Aggregation...

path <- ""
results <- LFQService::workflow_MQ_protoV1(resDataStart,
                                           config,
                                           path,
                                           peptideFilterFunction = LFQService:::.workflow_MQ_filter_peptides_V3 )


#rmarkdown::render("MQSummary2.Rmd", params=list(configuration = config$clone(deep=TRUE), data = resDataStart), output_format = bookdown::pdf_document2())
protintensity <- medpolish_protein_quants( results$pepIntensityNormalized, results$config_pepIntensityNormalized )
xx <- protintensity("unnest")
data <- LFQService::applyToIntensityMatrix(xx$data, xx$config, .func = robust_scale)
stats_res <- summarize_cv(data, xx$config, all=FALSE)
data_wide <- inner_join(stats_res,toWideConfig(data, xx$config)$data, by = xx$config$table$hkeysLevel(), suffix =c(".factor",""))


lfq_write_table(data_wide, path=file.path(outputDir,paste0("r_",filename,"_Protein.csv")))
lfq_write_table(toWideConfig(data, xx$config)$annotation, path=file.path(outputDir,paste0("r_",filename,"_annotation.csv")))

##```{r corrplot, fig.cap="Correlation plot."}
##LFQService::plot_sample_correlation(dataTransformed, config)
##```

LFQService::render_MQSummary_rmd(data, xx$config$clone(deep=TRUE),
                                 project_id = "MaxQuant_p3175_o5772_QC",
                                 workunit_id = "200367",
                                 pep = FALSE,
                                 dest_path = outputDir,
                                 dest_file_name = paste0("r_",filename,"_Protein.pdf"), workdir = ".")

#rmarkdown::render("MQSummary2.Rmd", params=list(data =xx$data, configuration = xx$config), output_file = paste0("Protein_",mqzip,".pdf"))#,envir =new.env())
#break()
#}

