library(prolfqua)
library(tidyverse)

peptidef <- "peptides.txt"
proteinf <- "proteinGroups.txt"
dsf <- "dataset.csv"
REPEATED <- TRUE

unzip("../data/2062292.zip",peptidef)
unzip("../data/2062292.zip",proteinf)

stopifnot(file.exists(peptidef), file.exists(proteinf), file.exists(dsf))

protein <- prolfqua::tidyMQ_ProteinGroups(proteinf)
peptide <- prolfqua::tidyMQ_Peptides(peptidef)

annot <- read.csv("../data/datasetAVSA.csv")
annot <- read.csv("../data/datasetCompareControlRepeated.csv")

annot <- annot %>% dplyr::mutate(
  raw.file = gsub("^x|.d.zip$|.raw$","",
                  tolower(make.names(basename(annot$Relative.Path)))
  )
)


annot$Relative.Path <- NULL

proteinAnnot <- select(protein, proteinID, fasta.headers ) %>% distinct()
peptide <- dplyr::inner_join(annot, peptide)
peptide <- dplyr::inner_join(proteinAnnot, peptide, by = c(proteinID = "leading.razor.protein"))

################### annotations
GRP2 <- list()


GRP2$projectID <- NA
GRP2$projectName <- NA
GRP2$workunitID <- "MaxQuant"
GRP2$nrPeptides <- 2
GRP2$log2FCthreshold <- 1
GRP2$FDRthreshold <- 0.1

# Setup configuration

atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("proteinID")
atable$hierarchy[["peptide_Id"]] <- c("sequence")

##
atable$hierarchyDepth <- 1
atable$factors[["Experiment_"]] = "Experiment"
if (!is.null(annot$Subject) & REPEATED) {
  atable$factors[["Subject"]] = "Subject"

}
atable$factorDepth <- 1
## only factor 1 ''Experiment_'' is important

atable$setWorkIntensity("peptide.intensity")


# compute all possible 2 grps to avoid specifying reference.

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

