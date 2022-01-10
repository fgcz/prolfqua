library(readr)
library(tidyverse)
library(prolfqua)


yml <- yaml::read_yaml("config.yaml")

BFABRIC <- list()
BFABRIC$workunitID = yml$job_configuration$workunit_id
BFABRIC$workunitURL = paste0("https://fgcz-bfabric.uzh.ch/bfabric/workunit/show.html?id=",BFABRIC$workunitID,"&tab=details")
BFABRIC$projectID = yml$job_configuration$project_id
BFABRIC$orderID = yml$job_configuration$order_id
BFABRIC$inputID = yml$job_configuration$input[[1]][[1]]$resource_id
BFABRIC$inputURL = yml$job_configuration$input[[1]][[1]]$resource_url

BFABRIC$datasetID <- yml$application$parameters$datasetId


spc <- if ( yml$application$parameters$SpcInt == "spc") { TRUE } else {FALSE}
FCthreshold <- as.numeric( yml$application$parameters$FCthreshold )
SSthreshold <- as.numeric( yml$application$parameters$SAINTscore )


ZIPDIR = paste0("C",BFABRIC$projectID,"WU",BFABRIC$workunitID)
dir.create(ZIPDIR)

FDRthreshold <- 0.25
FCthreshold <- 2
SSthreshold <- 0.75
treat <- "FRAGPIPE_"


annotation <- readr::read_csv("dataset.csv")
pp <- prolfqua::tidy_MSFragger_combined_protein_V16("combined_protein.tsv")



annotation$raw.file <- basename(annotation$`Relative Path`)
annotation <- mutate(annotation, raw.file = paste0("x", gsub(".raw", "", tolower(raw.file))))
annotation$raw.file %in% pp$raw.file
annotation$`Relative Path` <- NULL
pdata <- inner_join(annotation, pp )
pdata <- pdata %>% filter(combined.total.peptides > 1)
ata <- prolfqua::AnalysisTableAnnotation$new()
ata$fileName = "raw.file"

ata$factors[["CorT"]] = "CONTROL"
ata$factors[["bait"]] = "BAIT"
ata$factorDepth <- 2

ata$hierarchy[["protein_Id"]] = "protein"

if (spc) {
    ata$workIntensity = "spectral.count"
} else {
    ata$workIntensity = "intensity"
}

config <- prolfqua::AnalysisConfiguration$new(ata)
sdata <- setup_analysis(pdata, config)
lfqdata <- LFQData$new(sdata, config)
lfqdata$remove_small_intensities()


# remove rev sequences
lfqdata$data <- lfqdata$data %>% filter(!grepl("^REV__|^CON__", protein_Id))
RESULTS <- list()
RESULTS$annotation <- lfqdata$factors()


intdata <- lfqdata$data
intdata <- inner_join(intdata , distinct( select(pdata, protein, protein.length)), by = c(protein_Id = "protein"))

bb <- protein_2localSaint(intdata, quantcolumn = lfqdata$config$table$getWorkIntensity())
RESULTS <- c(RESULTS, bb)

res <- runSaint(bb, spc = spc)
RESULTS <- c(RESULTS, res)

###################
ttt <- res$list
cse <- ContrastsSaintExpress$new(ttt)
pcse <- cse$get_Plotter(fcthreshold = log2(FCthreshold), bfdrthreshold = FDRthreshold)
####################

writexl::write_xlsx(RESULTS, path = file.path(ZIPDIR,paste0(treat, "_data.xlsx")))

if (FALSE) {
    sig <- cse$get_contrasts() %>%
        filter(.data$BFDR  < FDRthreshold & .data$log2FC  >  log2(FCthreshold))
} else {
    sig <- cse$get_contrasts() %>%
        filter(.data$SaintScore  >  SSthreshold & .data$log2FC  >  log2(FCthreshold))
}

tt <- lfqdata$get_Transformer()$log2()
tt <- tt$robscale()
lfqdata <- tt$lfq

ReportData <- list()
ReportData$BFABRIC <- BFABRIC
ReportData$SSthreshold <- SSthreshold
ReportData$FCthreshold <- FCthreshold
ReportData$spc <- spc
ReportData$lfqdata <- lfqdata
ReportData$sig <- sig
ReportData$cse <- cse
ReportData$pcse <- pcse


rm(list = setdiff(ls(), c("ReportData","ZIPDIR","treat")))

rmarkdown::render("SaintExpressReportMsFragger.Rmd",
                  params = list(sep = ReportData),
                  output_format = bookdown::html_document2())

file.copy("SaintExpressReportMsFragger.html",
 file.path(ZIPDIR, paste0(treat, "SaintExpressReportMsFragger.html")),
 overwrite = TRUE)

