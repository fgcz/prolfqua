#' extract intensities and annotations from MQ proteinGroups.txt
#' @export
#' @param MQProteinGroups data.frame generated with read.csv("peptide.txt",sep="\\t", stringsAsFactors=FALSE)
#' @examples
#' protein_txt <- system.file("samples/maxquant_txt/MSQC1/proteinGroups.txt",package = "LFQService")
#' protein_txt <- read.csv(protein_txt, header=TRUE, stringsAsFactors = FALSE, sep="\t")
#' mq_proteins <-tidyMQ_ProteinGroups(protein_txt)
#' head(mq_proteins)
tidyMQ_ProteinGroups <- function(MQProteinGroups){
  if(is.character(MQProteinGroups)){
    MQProteinGroups <- read.csv(MQProteinGroups, header=TRUE, stringsAsFactors = FALSE, sep="\t")
  }
  colnames(MQProteinGroups) <- tolower(colnames(MQProteinGroups))

  pint <- dplyr::select(MQProteinGroups, "protein.group.id" = "id", starts_with("intensity."))
  pintLFQ <- dplyr::select(MQProteinGroups, "protein.group.id" = "id", starts_with("lfq.intensity."))
  meta <- dplyr::select(MQProteinGroups,
                        "protein.ids" = "protein.ids",
                        "majority.protein.ids" = "majority.protein.ids",
                        "nr.peptides" = "peptides",
                        "fasta.headers",
                        "protein.group.id" = "id",
                        "protein.score" = "score"
  )

  pint <- pint %>%
    gather(key="raw.file", value="mq.protein.intensity", starts_with("intensity.")) %>%
    mutate(raw.file = gsub("intensity.","",raw.file))

  pintLFQ <- pintLFQ %>%
    gather(key="raw.file", value="mq.protein.lfq.intensity", starts_with("lfq.intensity.")) %>%
    mutate(raw.file = gsub("lfq.intensity.","",raw.file))

  pint <- inner_join(pint, pintLFQ , by=c("protein.group.id","raw.file"))
  res <- inner_join(meta, pint , by="protein.group.id")
  return(res)
}


#' read evidence file
#' @export
#' @examples
#' library(tidyverse)
#' evidence_txt <- system.file("samples/maxquant_txt/MSQC1/evidence.txt",package = "LFQService")
#' evidence_txt <- read.csv(evidence_txt, header=TRUE, stringsAsFactors = FALSE, sep="\t")
#' mq_evidence <- tidyMQ_Evidence(evidence_txt)
tidyMQ_Evidence <- function(Evidence){
  if(is.character(Evidence)){
    Evidence <- read.csv(Evidence, header=TRUE, stringsAsFactors = FALSE, sep="\t")
  }
  colnames(Evidence) <- tolower(colnames(Evidence))
  res <- dplyr::select(Evidence,
                       "evidence.id" = "id",
                       "peptide.id",
                       "raw.file",
                       "protein.group.id"="protein.group.ids",
                       "evidence.score" = "score",
                       "delta.score",
                       "calibrated.retention.time",
                       "charge",
                       "mass",
                       "ms.ms.count",
                       "ms.ms.scan.number",
                       "evidence.intensity" = "intensity")
  res %>% mutate(raw.file = tolower(raw.file)) -> res
  res$proteotypic <-!grepl(";",res$protein.group.id)
  res <- res %>% separate_rows(protein.group.id, sep=";",convert =TRUE)
  return(res)
}
#' Generating mq all level file.
#'
#' @export
#' @examples
#'
#' txt_directory <- system.file("samples/maxquant_txt/MSQC1", package = "LFQService")
#' allData <- tidyMQ_All(txt_directory)
#' zip_archive <- "inst/samples/maxquant_txt/twoGroup3Reps.zip"
#' res <- tidyMQ_All(zip_archive)
tidyMQ_All <- function(txt_directory){
  if(grepl("\\.zip$",txt_directory)){
    proteins_txt <- read.csv(unz(txt_directory,"proteinGroups.txt"),
                             header=TRUE, sep="\t", stringsAsFactors = FALSE)
    peptides_txt <- read.csv(unz(txt_directory,"peptides.txt"),
                             header=TRUE, sep="\t", stringsAsFactors = FALSE)
    evidence_txt <- read.csv(unz(txt_directory,"evidence.txt"),
                             header=TRUE, sep="\t", stringsAsFactors = FALSE)
  }else{
    proteins_txt <- file.path(txt_directory, "proteinGroups.txt")
    peptides_txt <- file.path(txt_directory, "peptides.txt")
    evidence_txt <- file.path(txt_directory , "evidence.txt")
  }
  mq_proteins <- tidyMQ_ProteinGroups(proteins_txt)
  mq_peptides <- tidyMQ_Peptides(peptides_txt)
  mq_evidence <- tidyMQ_Evidence(evidence_txt)

  resProt_Pep_Evidence <- inner_join(resProt_Pep, mq_evidence, by = c("protein.group.id", "raw.file", "peptide.id"))
  return(resProt_Pep_Evidence)
}

#' Generating mq all level file.
#'
#' @export
#' @examples
#'
#' txt_directory <- system.file("samples/maxquant_txt/MSQC1", package = "LFQService")
#' allData <- tidyMQ_All(txt_directory)
#' zip_archive <- "inst/samples/maxquant_txt/twoGroup3Reps.zip"
#' res <- tidyMQ_PeptideProtein(zip_archive)
tidyMQ_PeptideProtein <- function(txt_directory, .all = FALSE){
  if(grepl("\\.zip$",txt_directory)){
    proteins_txt <- read.csv(unz(txt_directory,"proteinGroups.txt"),
                             header=TRUE, sep="\t", stringsAsFactors = FALSE)
    peptides_txt <- read.csv(unz(txt_directory,"peptides.txt"),
                             header=TRUE, sep="\t", stringsAsFactors = FALSE)
    mod_spec_peptides_txt <- read.csv(unz(txt_directory,"modificationSpecificPeptides.txt"),
                                      header=TRUE, sep="\t", stringsAsFactors = FALSE)

  }else{
    proteins_txt <- file.path(txt_directory, "proteinGroups.txt")
    peptides_txt <- file.path(txt_directory, "peptides.txt")
    mod_spec_peptides_txt <- file.path(txt_directory, "modificationSpecificPeptides.txt")
  }
  mq_proteins <- tidyMQ_ProteinGroups(proteins_txt)
  mq_peptides <- tidyMQ_Peptides(peptides_txt)
  mq_modSpecPeptides <- tidyMQ_modificationSpecificPeptides(mod_spec_peptides_txt)

  resProt_Pep <- inner_join(mq_proteins,mq_peptides, by = c("protein.group.id", "raw.file"))

  if(.all){
    return(list(resProt_Pep = resProt_Pep,
                mq_proteins = mq_proteins,
                mq_peptides = mq_peptides,
                mq_modSpecPeptides = mq_modSpecPeptides))
  }else{
    return(resProt_Pep)
  }
}

#' parse MQ modificationSpecificPeptides.txt
#' @export
#' @param MQPeptides data.frame generated with read.csv("peptide.txt",sep="\\t", stringsAsFactors=FALSE)
#' @examples
#' library(tidyverse)
#' peptides_txt <- "c:/Users/wewol/Dropbox/DataAnalysis/p2621_HumanAgeInteraction/data/721705/modificationSpecificPeptides.txt"
#' peptides_txt <- read.csv(peptides_txt, header=TRUE, stringsAsFactors = FALSE, sep="\t")
#' MQPeptides <- peptides_txt
#' View(MQPeptides)
#' mq_peptides <- tidyMQ_modificationSpecificPeptides(peptides_txt)
#'
#' head(mq_peptides)
#'
tidyMQ_modificationSpecificPeptides <- function(MQPeptides){
  if(is.character(MQPeptides)){
    MQPeptides <- read.csv(MQPeptides, header=TRUE, stringsAsFactors = FALSE, sep="\t")
  }
  colnames(MQPeptides) <- tolower(colnames(MQPeptides))
  #return(MQPeptides)
  sc <- sym("potential.contaminant")
  colnames(MQPeptides)
  meta <- dplyr::select(MQPeptides,
                        "mod.spec.peptide.id" = "id",
                        "peptide.id",
                        "sequence",
                        "modifications",
                        "proteins",
                        "protein.group.id"="protein.group.ids",
                        "mass",
                        "retention.time",
                        "peptide.score" ="score",
                        "delta.score",
                        "pep",
                        "missed.cleavages",
                        "unique.groups" = "unique..groups.",
                        "unique.proteins" = "unique..proteins.",
                        "potential.contaminant" = ends_with("contaminant")) %>%
    mutate(!!"potential.contaminant" := case_when( !!sc == "" ~ FALSE, !!sc == "+" ~ TRUE)) %>%
    mutate(!!"unique.groups" := case_when( !!sym("unique.groups") == "yes" ~ TRUE,
                                           !!sym("unique.groups") == "no" ~ FALSE))


  pint <- dplyr::select(MQPeptides,"mod.spec.peptide.id"= "id", starts_with("intensity."))
  PepIntensities <- pint %>%
    gather(key="raw.file", value="mod.spec.peptide.intensity", starts_with("intensity.")) %>%
    mutate(raw.file = gsub("intensity.","",raw.file))

  idtype <- dplyr::select(MQPeptides, "mod.spec.peptide.id"="id", starts_with("identification.type."))
  if(ncol(idtype) > 1){ # if only one file no id type is provided
    PepIDType <- idtype %>%
      gather(key="raw.file", value="mod.spec.id.type", starts_with("identification.type.")) %>%
      mutate(raw.file = gsub("identification.type.","",raw.file))
    PepIntensities <-inner_join(PepIntensities,PepIDType, by=c("mod.spec.peptide.id", "raw.file" ))
  }else{
    PepIntensities$id.type <- "By MS/MS"
  }
  xx <- inner_join(meta , PepIntensities, by="mod.spec.peptide.id")
  xx$proteotypic <-!grepl(";",xx$protein.group.id)
  xx <- xx %>% separate_rows(protein.group.id, sep=";",convert =TRUE)
  return(xx)
}

#' parse MQ peptides.txt
#' @export
#' @param MQPeptides data.frame generated with read.csv("peptide.txt",sep="\\t", stringsAsFactors=FALSE)
#' @examples
#' library(tidyverse)
#' peptide_txt <- "c:/Users/wewol/Dropbox/DataAnalysis/p2621_HumanAgeInteraction/data/721705/peptides.txt"
#' #peptide_txt <- system.file("samples/maxquant_txt/MSQC1/peptides.txt",package = "LFQService")
#' peptides_txt <- read.csv(peptide_txt, header=TRUE, stringsAsFactors = FALSE, sep="\t")
#' mq_peptides <-tidyMQ_Peptides(peptides_txt)
#' peptides_txt <- system.file("samples/maxquant_txt/tiny/peptides.txt",package = "LFQService")
#' peptides_txt <- read.csv(peptides_txt, header=TRUE, stringsAsFactors = FALSE, sep="\t")
#' tmp <-paste(peptides_txt$Evidence.IDs, collapse = ";")
#' tmp <- strsplit(tmp, ";")
#' length(unique(tmp[[1]]))
#'
#' mq_peptides <-tidyMQ_Peptides(peptides_txt)
#' head(mq_peptides)
tidyMQ_Peptides <- function(MQPeptides){
  if(is.character(MQPeptides)){
    MQPeptides <- read.csv(MQPeptides, header=TRUE, stringsAsFactors = FALSE, sep="\t")
  }
  colnames(MQPeptides) <- tolower(colnames(MQPeptides))
  sc <- sym("potential.contaminant")
  meta <- dplyr::select(MQPeptides,
                        "peptide.id" = "id",
                        "sequence",
                        "proteins",
                        "leading.razor.protein",
                        "protein.group.id"="protein.group.ids",
                        "peptide.score" ="score",
                        "pep",
                        "missed.cleavages",
                        "unique.groups" = "unique..groups.",
                        "potential.contaminant" = ends_with("contaminant")) %>%
    mutate(!!"potential.contaminant" := case_when( !!sc == "" ~ FALSE, !!sc == "+" ~ TRUE)) %>%
    mutate(!!"unique.groups" := case_when( !!sym("unique.groups") == "yes" ~ TRUE,
                                           !!sym("unique.groups") == "no" ~ FALSE))

  pint <- dplyr::select(MQPeptides,"peptide.id"= "id", starts_with("intensity."))

  PepIntensities <- pint %>%
    gather(key="raw.file", value="peptide.intensity", starts_with("intensity.")) %>%
    mutate(raw.file = gsub("intensity.","",raw.file))

  idtype <- dplyr::select(MQPeptides, "peptide.id"="id", starts_with("identification.type."))
  if(ncol(idtype) > 1){ # if only one file no id type is provided
    PepIDType <- idtype %>%
      gather(key="raw.file", value="id.type", starts_with("identification.type.")) %>%
      mutate(raw.file = gsub("identification.type.","",raw.file))
    PepIntensities <-inner_join(PepIntensities,PepIDType, by=c("peptide.id", "raw.file" ))
  }else{
    PepIntensities$id.type <- "By MS/MS"
  }
  xx<-inner_join(meta , PepIntensities, by="peptide.id")

  xx$proteotypic <-!grepl(";",xx$protein.group.id)
  xx <- xx %>% separate_rows(protein.group.id, sep=";",convert =TRUE)
  return(xx)
}

#' parse MQ allPeptides.txt
#' @export
#' @param MQPeptides data.frame generated with read.csv("peptide.txt",sep="\\t", stringsAsFactors=FALSE)
#' @examples
#' peptides_txt <- "c:/Users/wewol/Dropbox/DataAnalysis/p2621_HumanAgeInteraction/data/721705/allPeptides.txt"
#' peptides_txt <- read.csv(peptides_txt, header=TRUE, stringsAsFactors = FALSE, sep="\t")
#' MQPeptides <- peptides_txt
#' head(MQPeptides)
#' mq_peptides <- tidyMQ_allPeptides(peptides_txt)
#' dim(mq_peptides)
#' mq_peptides %>% dplyr::filter(ms.ms.count != 0) -> idid
#' dim(idid)
#' mq_peptides %>% dplyr::filter(sequence != "") -> idid
#' head(idid)
#' unique(idid$modified.sequence)
#' head(mq_peptides)
#'
tidyMQ_allPeptides <- function(MQPeptides){
  if(is.character(MQPeptides)){
    MQPeptides <- read.csv(MQPeptides, header=TRUE, stringsAsFactors = FALSE, sep="\t")
  }
  colnames(MQPeptides) <- tolower(colnames(MQPeptides))
  colnames(MQPeptides)
  xx <- dplyr::select(MQPeptides,
                        "raw.file",
                        "type",
                        "charge",
                        "m.z",
                        "mass",
                        "retention.time",
                        "sequence",
                        "modified.sequence",
                        "proteins",
                        "peptide.score" ="score",
                        "intensity",
                        "ms.ms.count") %>%
    mutate(sequence = str_trim(sequence), modified.sequence = str_trim(modified.sequence))
  return(xx)
}


