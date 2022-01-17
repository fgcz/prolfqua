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


cpf <- "combined_protein.tsv"
dsf <- "dataset.csv"

message("Checking input files for two-group prolfqua script ...")
stopifnot(file.exists(cpf), file.exists(dsf))


protein <- as_tibble(
  read.csv(cpf,
           header = TRUE,
           sep = "\t",
           stringsAsFactors = FALSE))

annot <- read.csv(dsf)
protein <- prolfqua::tidy_MSFragger_combined_protein_V16(protein)

annot <- annot |> dplyr::mutate(raw.file =
                                   make.names(paste0("x",tolower(gsub(".d.zip$|.raw$","",basename(Relative.Path))))))


annot$Relative.Path <- NULL
protein <- dplyr::inner_join(annot, protein)
protein <- protein |> dplyr::filter(unique.spectral.count > 1)


################### annotations
GRP2 <- list()


GRP2$projectID <- PROJECTID
GRP2$orderID <- ORDERID
GRP2$workunitID <- WORKUNITID

GRP2$Software <- "FragPipe IonQuant"

GRP2$inputID <- INPUT_ID
GRP2$inputURL <- INPUT_URL
GRP2$nrPeptides <- 2
GRP2$log2FCthreshold <- 1
GRP2$FDRthreshold <- 0.1

# Setup configuration

atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("entry.name")
##
atable$hierarchyDepth <- 1
atable$factors[["Experiment_"]] = "Experiment"
## only factor 1 ''group'' is important
atable$factorDepth <- 1
atable$setWorkIntensity("intensity")

# Compute all possible 2 grps to avoid specifying reference.

levels <- annot$Experiment |> unique()

for (i in 1:length(levels)) {
  for (j in 1:length(levels)) {
    if (i != j) {

      cat(levels[i], levels[j], "\n")
      GRP2$Contrasts <- paste0("Experiment_",levels[i], " - ", "Experiment_",levels[j])
      names(GRP2$Contrasts) <- paste0("Experiment" , levels[i], "_vs_", levels[j])
      message(GRP2$Contrasts)
      outpath <- paste0("Experiment_" , levels[i], "_vs_", levels[j])
      proteinF <- peptide |> dplyr::filter( .data$Experiment == levels[i] | .data$Experiment == levels[j])

      grp2 <- prolfqua::make2grpReport(proteinF, atable, GRP2, protein_annot = "fasta.headers",
                                       remove = TRUE)

      prolfqua::write_2GRP(grp2, outpath = outpath)
      prolfqua::render_2GRP(grp2, outpath = outpath, htmlname = outpath)

    }
  }
}

