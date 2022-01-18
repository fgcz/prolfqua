library(prolfqua)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  ymlfile <- "config_WU273293.yaml"
} else {
  ymlfile <- args[1]

}
yml = yaml::read_yaml(ymlfile)



WORKUNITID = yml$job_configuration$workunit_id
PROJECTID = yml$job_configuration$project_id
ORDERID = yml$job_configuration$order_id
INPUT_ID = yml$job_configuration$input[[1]][[1]]$resource_id
INPUT_URL = yml$job_configuration$input[[1]][[1]]$resource_url

peptidef <- "peptides.txt"
proteinf <- "proteinGroups.txt"
dsf <- "dataset.csv"
REPEATED <- TRUE


stopifnot(file.exists(peptidef), file.exists(proteinf), file.exists(dsf))

protein <- prolfqua::tidyMQ_ProteinGroups(proteinf)
peptide <- prolfqua::tidyMQ_Peptides(peptidef)
annot <- read.csv(dsf)

annot <- annot |> dplyr::mutate(
  raw.file = gsub("^x|.d.zip$|.raw$","",
                  tolower(make.names(basename(annot$Relative.Path)))
  )
)


annot$Relative.Path <- NULL

proteinAnnot <- dplyr::select(protein, proteinID, fasta.headers ) |> distinct()
peptide <- dplyr::inner_join(annot, peptide)
peptide <- dplyr::inner_join(proteinAnnot, peptide, by = c(proteinID = "leading.razor.protein"))

################### annotations
GRP2 <- list()

GRP2$projectID <- PROJECTID
GRP2$orderID <- ORDERID
GRP2$workunitID <- WORKUNITID

GRP2$Software <- "MaxQuant"

GRP2$inputID <- INPUT_ID
GRP2$inputURL <- INPUT_URL
GRP2$nrPeptides <- 2
GRP2$log2FCthreshold <- 1
GRP2$FDRthreshold <- 0.1

# Setup configuration

atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("proteinID")
atable$hierarchy[["peptide_Id"]] <- c("sequence")

#
atable$hierarchyDepth <- 1
atable$factors[["Experiment_"]] = "Experiment"

if (!is.null(annot$Subject) & REPEATED) {
  atable$factors[["Subject"]] = "Subject"
}
atable$factorDepth <- 1
atable$setWorkIntensity("peptide.intensity")


if (FALSE) {
  ps <- prolfqua::ProjectStructure$new(outpath = ".",
                                       project_Id = "",
                                       workunit_Id = basename(getwd()),
                                       order_Id = "",
                                       inputAnnotation = NULL,
                                       inputData = NULL)

  prolfqua::render_MQSummary_rmd(lfqdata$data,
                                 config$clone(deep = TRUE),
                                 ps, format = "html")
}



# Compute all possible 2 Grps to avoid specifying reference.
levels <- annot$Experiment |> unique()
outdir <- "xyz"
dir.create(outdir)


for (i in 1:length(levels)) {
  for (j in 1:length(levels)) {
    if (i != j) {
      i <- 1
      j <- 2
      cat(levels[i], levels[j], "\n")
      GRP2$Contrasts <- paste0("Experiment_",levels[i], " - ", "Experiment_",levels[j])
      names(GRP2$Contrasts) <- paste0("Experiment" , levels[i], "_vs_", levels[j])
      message(GRP2$Contrasts)
      outpath <- file.path( outdir, paste0("Experiment_" , levels[i], "_vs_", levels[j]))
      proteinF <- peptide |> dplyr::filter( .data$Experiment == levels[i] | .data$Experiment == levels[j])

      debug(prolfqua::make2grpReport)
      grp2 <- prolfqua::make2grpReport(proteinF,
                                       atable,
                                       GRP2,
                                       protein_annot = "fasta.headers",
                                       remove = TRUE)


      prolfqua::write_2GRP(grp2, outpath = outpath)
      prolfqua::render_2GRP(grp2, outpath = outpath, htmlname = outpath)

    }
  }
}

