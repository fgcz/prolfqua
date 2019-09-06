# used in p3147 for Marcus Masek

####################### End of configuration #################################

assign("lfq_write_format", "xlsx", envir = .GlobalEnv)

# create result structure
qc_path <- file.path(outpath, "qc_results")
modelling_path <- file.path(outpath, "modelling_results")


if(!dir.exists(outpath)){
  dir.create(outpath)
}

if(!dir.exists(qc_path)){
  dir.create(qc_path)
}

if(!dir.exists(modelling_path)){
  dir.create(modelling_path)
}



## read the data

resPepProtAnnot <- tidyMQ_modificationSpecificPeptides(inputMQfile)

pdf(file.path(qc_path, "retention_time_plot.pdf"), height = 15)
resPepProtVis <- resPepProtAnnot %>% filter(mod.peptide.intensity > 4)
ggplot(resPepProtVis, aes(x = retention.time, y= log2(mod.peptide.intensity))) + geom_point(alpha=1/20, size=0.3) + facet_wrap(~raw.file, ncol=2 )
dev.off()


peptidestxt <- tidyMQ_Peptides(inputMQfile)
resPepProtAnnot <- tidyMQ_from_modSpecific_to_peptide(resPepProtAnnot, peptidestxt)
resPepProtAnnot <- tidyMQ_top_protein_name(resPepProtAnnot)



# add annotation
annotation <- readxl::read_xlsx(inputAnntation)
noData <- annotation[!annotation$raw.file %in% resPepProtAnnot$raw.file,]


if(nrow(noData)){
  message("some files in annotation have no measurements")
  message(paste(noData,collapse = " "))
}

measSamples <- unique(resPepProtAnnot$raw.file)
noAnnot <- measSamples[! measSamples%in% annotation$raw.file]

if(length(noAnnot) > 0 ){
  message("some measured samples have no annotation!")
  message(paste(noAnnot,collapse = " "))
}

resPepProtAnnot <- inner_join(annotation, resPepProtAnnot, by= "raw.file")

####




###  Setup analysis ####
#resPepProtAnnot <- resPepProtAnnot %>% filter(raw.file !=  "intensity")
resPepProtAnnot <- resPepProtAnnot %>% filter(reverse !=  TRUE)
resPepProtAnnot$isotope = "light"

resDataStart <- setup_analysis(resPepProtAnnot, config)

resDataStart <- remove_small_intensities( resDataStart, config, threshold = 4 ) %>%
  completeCases(config)


# fitler the data

results <- LFQService::workflow_MQ_protoV1(
  resDataStart,
  config,
  outpath,
  peptideFilterFunction = LFQService:::.workflow_MQ_filter_peptides_V3 )

wideFRAME <- LFQService::toWideConfig(results$pepIntensityNormalized,
                                      results$config_pepIntensityNormalized)


lfq_write_table(separate_hierarchy(wideFRAME$data,
                                   results$config_pepIntensityNormalized),
                file.path(qc_path, "peptide_intensities.csv"))


if(!DEBUG){
  LFQService::render_MQSummary_rmd(resDataStart,
                                   config$clone(deep=TRUE),
                                   pep=TRUE,
                                   workdir = ".",
                                   dest_path = qc_path,
                                   dest_file_name="peptide_intensities_qc",
                                   format = "html")
}


### PROTEIN QUANTIFICATION ####
config$table$hierarchyLevel <- 1
protintensity <- medpolish_protein_quants( results$pepIntensityNormalized,
                                           results$config_pepIntensityNormalized )

prot_summary <- summarizeHierarchy(results$pepIntensityNormalized,
                                   results$config_pepIntensityNormalized)


protein_quants_write(protintensity,
                     qc_path)

figs <- protintensity("plot")

pdf(file.path(qc_path, "protein_intensities_inference_figures.pdf"))
lapply(figs$plot[1:10],print)
dev.off()



# render protein quantification reprot
if(!DEBUG){
  LFQService::render_MQSummary_rmd(protintensity("unnest")$data,
                                   protintensity("unnest")$config$clone(deep=TRUE),
                                   pep=FALSE,
                                   workdir = ".",
                                   dest_path = qc_path,
                                   dest_file_name="protein_intensities_qc",
                                   format="html"
  )
}

#################################################
### Do missing value imputation

filteredPep <- results$pepIntensityNormalized
config <- results$config_pepIntensityNormalized

res_contrasts_imputed <- workflow_missigness_impute_contrasts(filteredPep,
                                                              config,
                                                              Contrasts)


## linear models ----

pepConfig <- config
# first model ----

model <- paste0(pepConfig$table$getWorkIntensity() ,model)
modelFunction <- make_custom_model_lmer( model, model_name = "Model")


### make contrasts -----

modellingResult_fun <- workflow_model_analyse(filteredPep,
                                              modelFunction,
                                              subject_Id = pepConfig$table$hkeysLevel())

modellingResult <- modellingResult_fun()
modellingResult_fun(modelling_path)


m <- get_complete_model_fit(modellingResult$modellingResult$modelProtein)

#factor_contrasts <- linfct_factors_contrasts(m$linear_model[[1]])

linfct <- linfct_from_model(m$linear_model[[1]], as_list = FALSE)
linfct_A <- linfct_matrix_contrasts(linfct, Contrasts)


if(DEBUG){
  pdf(file.path(modelling_path,"Linear_functions.pdf"), width=18, height=10)
  quantable::imageWithLabels(t(linfct_A),
                             col = quantable::getBlueWhiteRed(),
                             marLeft = c(8,10,4.1,2.1))
  dev.off()
}

#rownames(factor_contrasts$linfct_interactions)

modelProteinF <- modellingResult$modellingResult$modelProtein %>% dplyr::filter(exists_lmer == TRUE)

res_contrasts <- workflow_contrasts_linfct(modelProteinF ,
                                           linfct_A,
                                           pepConfig,
                                           modelName = modellingResult$modellingResult$modelName,
                                           prefix =  "Contrasts",
                                           contrastfun = modelFunction$contrast_fun)

#names(xx)
#View(xx$contrast_result)

xx <- res_contrasts(modelling_path,columns = reportColumns)
xx_imputed <- res_contrasts_imputed("long",what = "contrasts")

merge_contrasts_results <- function(contrast_minimal,
                                    xx_imputed, subject_Id){
  res <- right_join(contrast_minimal, xx_imputed, by=c(subject_Id,"lhs" = "contrast"))
  res <- res %>% rename(contrast = lhs)
  res <- res %>% mutate(pseudo_estimate = case_when(is.na(estimate) ~ imputed, TRUE ~ estimate))
  res <- res %>% mutate(is_pseudo_estimate = case_when(is.na(estimate) ~ TRUE, TRUE ~ FALSE))

  for(column in reportColumns){
    print(column)
    res <- res %>% mutate(!!sym(paste0("pseudo_",column)) := case_when(is.na(estimate) ~ 0, TRUE ~ !!sym(column)))
  }
  res <- res %>% select(-imputed, -meanArea)
  return(res)
}

contrast_results <- merge_contrasts_results(xx$contrast_minimal, xx_imputed,
                                            subject_Id = pepConfig$table$hkeysLevel())

lfq_write_table(separate_hierarchy(contrast_results, config), path = file.path(modelling_path, "foldchange_estimates.csv"))
