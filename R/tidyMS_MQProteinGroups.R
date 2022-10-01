#' Methods for reading  MaxQuant outputs
#'
#' Convert MaxQuant outputs into tidy tables. For more details see functions listed in the see also section.
#'
#' @family MaxQuant
#' @name MaxQuant
NULL

#' extract intensities and annotations from MQ proteinGroups.txt
#' @export
#' @keywords internal
#' @param MQProteinGroups data.frame generated with read.csv("peptide.txt",sep="\\t", stringsAsFactors=FALSE)
#' @family MaxQuant
#'
#' @examples
#' protein_txt <- system_file("samples/maxquant_txt/tiny2.zip",package = "prolfqua")
#' protein_txt <- read.csv(unz(protein_txt,"proteinGroups.txt"), header=TRUE, stringsAsFactors = FALSE, sep="\t")
#' mq_proteins <-tidyMQ_ProteinGroups(protein_txt)
#'
tidyMQ_ProteinGroups <- function(MQProteinGroups) {
  if (is.character(MQProteinGroups)) {
    if (grepl("\\.zip$",tolower(MQProteinGroups))) {
      message(MQProteinGroups)
      MQProteinGroups <- read.csv(unz(MQProteinGroups,"proteinGroups.txt"),
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    }else{
      MQProteinGroups <- read.csv(MQProteinGroups, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    }
  }
  colnames(MQProteinGroups) <- tolower(colnames(MQProteinGroups))

  pint <- dplyr::select(MQProteinGroups, "protein.group.id" = "id", starts_with("intensity."))
  pintLFQ <- dplyr::select(MQProteinGroups, "protein.group.id" = "id", starts_with("lfq.intensity."))
  pintMSCount <- dplyr::select(MQProteinGroups, "protein.group.id" = "id", starts_with("ms.ms.count."))


  meta <- dplyr::select(MQProteinGroups,
                        "protein.ids" = "protein.ids",
                        "majority.protein.ids" = "majority.protein.ids",
                        "nr.peptides" = "peptides",
                        "fasta.headers",
                        "protein.group.id" = "id",
                        "protein.score" = one_of("score")
  )
  meta <- meta |> mutate( proteinID = gsub(";.*$","", .data$majority.protein.ids) )
  pint <- pint |>
    tidyr::gather(key = "raw.file", value = "mq.protein.intensity", starts_with("intensity.")) |>
    dplyr::mutate(raw.file = gsub("intensity.","",.data$raw.file))

  pintLFQ <- pintLFQ |>
    tidyr::gather(key = "raw.file", value = "mq.protein.lfq.intensity", starts_with("lfq.intensity.")) |>
    dplyr::mutate(raw.file = gsub("lfq.intensity.","",.data$raw.file))

  pintCount <- pintMSCount |>
    tidyr::gather(key = "raw.file", value = "mq.protein.ms.ms.count", starts_with("ms.ms.count.")) |>
    dplyr::mutate(raw.file = gsub("ms.ms.count.","",.data$raw.file))


  pint <- dplyr::inner_join(pint, pintLFQ , by = c("protein.group.id","raw.file"))
  pint <- dplyr::inner_join(pint, pintCount , by = c("protein.group.id","raw.file"))


  res <- dplyr::inner_join(meta, pint , by = "protein.group.id")
  return(res)
}


#' read evidence file
#' @param Evidence MQ evidence file or zip archive with evidence file
#' @export
#' @keywords internal
#' @family MaxQuant
#' @examples
#'
#' evidence_txt <- system_file("samples/maxquant_txt/tiny2.zip",package = "prolfqua")
#' evidence_txt <- read.csv(unz(evidence_txt,"evidence.txt"), header=TRUE, stringsAsFactors = FALSE, sep="\t")
#' mq_evidence <- tidyMQ_Evidence(evidence_txt)
tidyMQ_Evidence <- function(Evidence){
  if (is.character(Evidence)) {
    if (grepl("\\.zip$",Evidence)) {
      message(Evidence)
      Evidence <- read.csv(unz(Evidence,"evidence.txt"),
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    }else{
      Evidence <- read.csv(Evidence, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    }
  }
  colnames(Evidence) <- tolower(colnames(Evidence))
  res <- dplyr::select(Evidence,
                       "evidence.id" = "id",
                       "peptide.id",
                       "raw.file",
                       "protein.group.id" = "protein.group.ids",
                       "mod.peptide.id" = "mod..peptide.id",
                       "leading.razor.protein",
                       "evidence.score" = "score",
                       "delta.score",
                       "pep",
                       "calibrated.retention.time",
                       "retention.time",
                       "retention.length",
                       "charge",
                       "mass",
                       "ms.ms.count",
                       "ms.ms.scan.number",
                       "evidence.intensity" = "intensity",
                       "modifications",
                       "modified.sequence",
                       "missed.cleavages",
                       "reverse"
                       )
  res <- res |> dplyr::mutate(reverse = dplyr::case_when(reverse == "+" ~ TRUE, TRUE ~ FALSE))
  res |> dplyr::mutate(raw.file = tolower(.data$raw.file)) -> res
  res$proteotypic <- !grepl(";",res$protein.group.id)
  res <- res |> separate_rows(.data$protein.group.id, sep = ";",convert  = TRUE)
  return(res)
}


#' parse MQ peptides.txt
#' @param MQPeptides data.frame generated with read.csv("peptide.txt",sep = "\\t", stringsAsFactors = FALSE)
#' @family MaxQuant
#' @export
#' @keywords internal
#' @examples
#'
#'
#' peptide_txt <- system_file("samples/maxquant_txt/tiny2.zip",
#'  package = "prolfqua")
#'
#' peptides_txt <- read.csv(unz(peptide_txt, "peptides.txt"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#' mq_peptides <- tidyMQ_Peptides(peptides_txt)
#'
tidyMQ_Peptides <- function(MQPeptides, proteotypic_only = TRUE){
  if (is.character(MQPeptides)) {
    if (grepl("\\.zip$",tolower(MQPeptides))) {
      MQPeptides <- read.csv(unz(MQPeptides,"peptides.txt"),
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    }else{
      MQPeptides <- read.csv(MQPeptides, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    }
  }
  colnames(MQPeptides) <- tolower(colnames(MQPeptides))
  sc <- sym("potential.contaminant")
  meta <- dplyr::select(MQPeptides,
                        "peptide.id" = "id",
                        "sequence",
                        "proteins",
                        "leading.razor.protein",
                        "protein.group.id" = "protein.group.ids",
                        "peptide.score"  = "score",
                        "pep",
                        "ms.ms.count" = "ms.ms.count",
                        dplyr::one_of("missed.cleavages"),
                        "unique.groups" = "unique..groups.",
                        "potential.contaminant" = ends_with("contaminant"),
                        "reverse" = "reverse"
                        ) |>
    dplyr::mutate(!!"potential.contaminant" := case_when( !!sc == "" ~ FALSE, !!sc == "+" ~ TRUE)) |>
    dplyr::mutate(!!"unique.groups" := case_when( !!sym("unique.groups") == "yes" ~ TRUE,
                                                  !!sym("unique.groups") == "no" ~ FALSE)) |>
    dplyr::mutate(!!"reverse" := case_when( !!sym("reverse") == "+" ~ TRUE,
                                            !!sym("reverse") == "" ~ FALSE))

  pint <- dplyr::select(MQPeptides,"peptide.id" =  "id", starts_with("intensity."))

  PepIntensities <- pint |>
    tidyr::gather(key = "raw.file", value = "peptide.intensity", starts_with("intensity.")) |>
    dplyr::mutate(raw.file = gsub("intensity.","",.data$raw.file))

  # review what happens here
  idtype <- dplyr::select(MQPeptides, "peptide.id" = "id", starts_with("identification.type."))
  if (ncol(idtype) > 1) { # if only one file no id type is provided
    PepIDType <- idtype |>
      tidyr::gather(key = "raw.file", value = "id.type", starts_with("identification.type.")) |>
      dplyr::mutate(raw.file = gsub("identification.type.","",.data$raw.file))
    PepIntensities <- dplyr::inner_join(PepIntensities,PepIDType, by = c("peptide.id", "raw.file" ))
  }else{
    PepIntensities$id.type <- "By MS/MS"
  }

  xx <- dplyr::inner_join(meta , PepIntensities, by = "peptide.id")
  xx$proteotypic <- !grepl(";", xx$protein.group.id)
  xx <- xx |> separate_rows( .data$protein.group.id, sep = ";", convert  = TRUE)

  xx <- xx |> mutate(proteins = case_when(proteins == "" ~ leading.razor.protein, TRUE ~ proteins))
  xx$isotope <- "light"
  if (proteotypic_only) {
    xx <- xx |> dplyr::filter(.data$proteotypic == TRUE)
  }
  return(xx)
}






