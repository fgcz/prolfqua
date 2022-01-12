# Author : Witold Wolski <wew@fgcz.ethz.ch>
# compatible with prolfqua 2.9.0 release available from https://github.com/wolski/prolfqua/releases/tag/v0.2.9


# Read b-fabric related information
yml <- yaml::read_yaml("config.yaml")

BFABRIC <- list()
BFABRIC$workunitID = yml$job_configuration$workunit_id
BFABRIC$workunitURL = paste0("https://fgcz-bfabric.uzh.ch/bfabric/workunit/show.html?id=",BFABRIC$workunitID,"&tab=details")
BFABRIC$projectID = yml$job_configuration$project_id
BFABRIC$orderID = yml$job_configuration$order_id
BFABRIC$inputID = yml$job_configuration$input[[1]][[1]]$resource_id
BFABRIC$inputURL = yml$job_configuration$input[[1]][[1]]$resource_url
BFABRIC$datasetID <- yml$application$parameters$datasetId

ZIPDIR = paste0("C",BFABRIC$projectID,"WU",BFABRIC$workunitID)
dir.create(ZIPDIR)


# list with data used with the markdown report
REPORTDATA <- list()

# Applciation parameters
REPORTDATA$spc <- if ( yml$application$parameters$SpcInt == "spc") { TRUE } else {FALSE}
REPORTDATA$FCthreshold <- as.numeric( yml$application$parameters$FCthreshold )
REPORTDATA$SSthreshold <- as.numeric( yml$application$parameters$SAINTscore )
REPORTDATA$FDRthreshold <- 0.25

# Prefix for exported files
treat <- "FRAGPIPE_"

# load data
annotation <- readr::read_csv("dataset.csv")
pp <- prolfqua::tidy_MSFragger_combined_protein_V16("combined_protein.tsv")

# attach annotation to combined_protein data
annotation$raw.file <- basename(annotation$`Relative Path`)
annotation <- dplyr::mutate(annotation, raw.file = paste0("x", gsub(".raw", "", tolower(raw.file))))
annotation$`Relative Path` <- NULL
stopifnot(sum(annotation$raw.file %in% pp$raw.file) > 0) # check that some files are annotated, if not exit script.
pdata <- dplyr::inner_join(annotation, pp )

# filter for more than 2 peptides per protein
pdata <- pdata |> dplyr::filter(combined.total.peptides > 1)

# configure prolfqua
ata <- prolfqua::AnalysisTableAnnotation$new()
ata$fileName = "raw.file"
ata$factors[["CorT"]] = "CONTROL"
ata$factors[["bait"]] = "BAIT"
ata$factorDepth <- 2

ata$hierarchy[["protein_Id"]] = "protein"

if (REPORTDATA$spc) {
    ata$workIntensity = "razor.spectral.count"
} else {
    ata$workIntensity = "razor.intensity"
}


config <- prolfqua::AnalysisConfiguration$new(ata)
sdata <- prolfqua::setup_analysis(pdata, config)
lfqdata <- prolfqua::LFQData$new(sdata, config)
lfqdata$remove_small_intensities()


# remove rev and contaminant sequences
lfqdata$data <- lfqdata$data |> dplyr::filter(!grepl("^REV__|^CON__", protein_Id))


RESULTS <- list() # RESULT is stored in excel table
RESULTS$annotation <- lfqdata$factors()

# Run Saint Analysis
intdata <- lfqdata$data
intdata <- dplyr::inner_join(intdata ,
                      dplyr::distinct( dplyr::select(pdata, protein, protein.length)),
                      by = c(protein_Id = "protein"))

bb <- prolfqua::protein_2localSaint(
  intdata,
  quantcolumn = lfqdata$config$table$getWorkIntensity())
RESULTS <- c(RESULTS, bb)
res <- prolfqua::runSaint(bb, spc = REPORTDATA$spc)
RESULTS <- c(RESULTS, res)

# write analysis results
writexl::write_xlsx(RESULTS, path = file.path(ZIPDIR,paste0(treat, "_data.xlsx")))

# Prepare result visualization and render report
cse <- prolfqua::ContrastsSaintExpress$new(res$list)
pcse <- cse$get_Plotter(fcthreshold = log2(REPORTDATA$FCthreshold),
                        bfdrthreshold = REPORTDATA$FDRthreshold)

sig <- cse$get_contrasts() |>
  dplyr::filter(.data$SaintScore  >  REPORTDATA$SSthreshold & .data$log2FC  >  log2(REPORTDATA$FCthreshold))

tt <- lfqdata$get_Transformer()$log2()
lfqdata <- tt$robscale()$lfq

REPORTDATA$BFABRIC <- BFABRIC
REPORTDATA$lfqdata <- lfqdata
REPORTDATA$sig <- sig
REPORTDATA$cse <- cse
REPORTDATA$pcse <- pcse


rm(list = setdiff(ls(), c("REPORTDATA","ZIPDIR","treat"))) # delete all variables not needed for rendering

rmarkdown::render("SaintExpressReportMsFragger.Rmd",
                  params = list(sep = REPORTDATA),
                  output_format = bookdown::html_document2(),
                  envir = new.env())

file.copy("SaintExpressReportMsFragger.html",
 file.path(ZIPDIR, paste0(treat, "SaintExpressReportMsFragger.html")),
 overwrite = TRUE)


