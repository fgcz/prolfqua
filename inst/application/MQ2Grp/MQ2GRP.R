library(prolfqua)
library(tidyverse)

peptide <- "peptide.txt"
dsf <- "dataset.csv"

message("Checking input files for two-group prolfqua script ...")
stopifnot(file.exists(cpf), file.exists(dsf))


protein <- as_tibble(
  read.csv(cpf,
           header = TRUE, sep = "\t",
           stringsAsFactors = FALSE))

annot <- read.csv(dsf)
protein <- prolfqua::tidy_MSFragger_combined_protein_V16(protein)
peptide <- prolfqua::tidyMQ_Peptides(peptide)

annot <- annot %>% dplyr::mutate(raw.file =
                                   make.names(paste0("x",tolower(gsub(".d.zip$|.raw$","",basename(Relative.Path))))))


annot$Relative.Path <- NULL

protein <- dplyr::inner_join(annot, peptide)

################### annotations
GRP2 <- list()


GRP2$projectID <- NA
GRP2$projectName <- NA
GRP2$workunitID <- "MSFragger IonQuant"
GRP2$nrPeptides <- 2
GRP2$log2FCthreshold <- 1
GRP2$FDRthreshold <- 0.1

# setup configuration

atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("entry.name")
##
atable$hierarchyDepth <- 1
atable$factors[["Experiment_"]] = "Experiment"
## only factor 1 ''Experiment_'' is important
atable$factorDepth <- 1
atable$setWorkIntensity("intensity")


# compute all possible 2 grps to avoid specifying reference.

levels <- annot$Experiment |> unique()

for (i in 1:length(levels)) {
  for (j in 1:length(levels)) {
    if (i != j) {
      cat(levels[i], levels[j], "\n")
      GRP2$Contrasts <- paste0("Experiment_",levels[i], " - ", "Experiment_",levels[j])
      names(GRP2$Contrasts) <- paste("Experiment" , levels[i], "vs", levels[j], sep = " ")
      message(GRP2$Contrasts)
      outpath <- paste("Experiment" , levels[i], "vs", levels[j], sep = "")
      proteinF <- protein |> dplyr::filter(.data$Experiment == levels[i] | .data$Experiment == levels[j])
      grp2 <- prolfqua::make2grpReport(proteinF, atable, GRP2, protein_annot = "description",
                                       remove = TRUE)

      prolfqua::write_2GRP(grp2, outpath = outpath)
      prolfqua::render_2GRP(grp2, outpath = outpath, htmlname = outpath)

    }
  }
}

