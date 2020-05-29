#!/usr/bin/Rscript

suppressMessages(library(readr))
suppressMessages(library(tidyverse))
suppressMessages(library(LFQService))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
library(docopt)

doc <- "Sample Size Report from MQ file

Usage:
  lfq_MQ_SampleSizeReport.R <mqzip> [--outdir=<outdir>] [--project_Id=<project_Id>] [--order_Id=<order_Id>] [--workunit=<workunit>]

Options:
  -o --outdir=<outdir> output directory [default: samplesize_qc_report]
  -p --project_Id=<project_Id> project identifier [default: NULL]
  -q --order_Id=<order_Id> order identifier [default: NULL]
  -w --workunit=<workunit> parent workunit [default: NULL]

Arguments:
  mqzip input file
"

opt <- docopt(doc)

cat("\nParameters used:\n\t",
    "     mqzip:", mqzip <- opt$mqzip, "\n\t",
    "result_dir:", outputDir <- opt[["--outdir"]], "\n\t",
    "project_Id:", project_Id <- opt[["--project_Id"]], "\n\t",
    "  order_Id:", order_Id <- opt[["--order_Id"]], "\n\t",
    "  workunit:", workunit <- opt[["--workunit"]], "\n\n\n")



if (!dir.exists(outputDir)) {
  if (!dir.create(outputDir)) {
    stop("\n could not create : ", outputDir, "\n")
  }
}

summarize_cv_raw_transformed <- function(resDataStart, config){
  config_tmp <- config$clone(deep = TRUE)
  stats_res <- summarize_cv(resDataStart, config_tmp)
  wide <- toWideConfig( resDataStart,config_tmp)
  stats_raw <- inner_join(stats_res, wide$data, by=config_tmp$table$hierarchyKeys())

  resDataStart <- transform_work_intensity(resDataStart, config_tmp, log2)
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
  atable$factorDepth <- 1
  anaparam <- AnalysisParameters$new()
  configuration <- AnalysisConfiguration$new(atable, anaparam)
  return(configuration)
}

config <- createMQProteinPeptideConfiguration()
config$project_Id <- project_Id
config$order_Id <- order_Id
config$workunit_Id <- workunit


resDataStart <- setup_analysis(resPepProtAnnot, config)
resDataStart <- remove_small_intensities(resDataStart, config) %>%
  complete_cases(config)


filename <- basename(tools::file_path_sans_ext(mqzip))


LFQService::render_MQSummary_rmd(
  resDataStart,
  config$clone(deep = TRUE),
  pep = TRUE,
  dest_path = outputDir,
  dest_file_name = paste0("r_", filename, "_Peptide"),
  workdir = ".",
  format = "html"
)



peptideStats <- summarize_cv_raw_transformed(resDataStart, config)
lfq_write_table(separate_hierarchy(peptideStats, config),
                path = file.path(outputDir, paste0("r_", filename, "_PeptideStats.csv")))


# Perform filtering and Protein Aggregation...


path <- ""
res <- filter_proteins_by_peptide_count(resDataStart, config)
results <- normalize_log2_robscale(res$data, res$config)


#rmarkdown::render("MQSummary2.Rmd", params=list(configuration = config$clone(deep=TRUE), data = resDataStart), output_format = bookdown::pdf_document2())
protintensity <- medpolish_protein_quants( results$data, results$config )
xx <- protintensity("unnest")
data <- LFQService::applyToIntensityMatrix(xx$data, xx$config, .func = robust_scale)
stats_res <- summarize_cv(data, xx$config, all = FALSE)
data_wide <-
  inner_join(
    stats_res,
    toWideConfig(data, xx$config)$data,
    by = xx$config$table$hkeysDepth(),
    suffix = c(".factor", "")
  )


lfq_write_table(
  separate_hierarchy(data_wide, config),
  path = file.path(outputDir, paste0("r_", filename, "_Protein.csv")),
  lfq_write_format = "xlsx"
)
lfq_write_table(
  toWideConfig(data, xx$config)$annotation,
  path = file.path(outputDir, paste0("r_", filename, "_annotation.csv")),
  lfq_write_format = "xlsx"
)

##```{r corrplot, fig.cap="Correlation plot."}
##LFQService::plot_sample_correlation(dataTransformed, config)
##```


LFQService::render_MQSummary_rmd(
  data,
  xx$config$clone(deep = TRUE),
  pep = FALSE,
  dest_path = outputDir,
  dest_file_name = paste0("r_", filename, "_Protein"),
  workdir = ".",
  format = "html"
)

#rmarkdown::render("MQSummary2.Rmd", params=list(data =xx$data, configuration = xx$config), output_file = paste0("Protein_",mqzip,".pdf"))#,envir =new.env())
#break()
#}

