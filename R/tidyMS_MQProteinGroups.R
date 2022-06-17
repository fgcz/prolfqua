#' methods for reading  MaxQuant outputs
#'
#' Convert outputs into tidy tables
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
#' library(prolfqua)
#' protein_txt <- system.file("samples/maxquant_txt/tiny2.zip",package = "prolfqua")
#' protein_txt <- read.csv(unz(protein_txt,"proteinGroups.txt"), header=TRUE, stringsAsFactors = FALSE, sep="\t")
#' mq_proteins <-tidyMQ_ProteinGroups(protein_txt)
#' head(mq_proteins)
#' plot(mq_proteins$mq.protein.intensity, mq_proteins$mq.protein.lfq.intensity, log="xy")
#' plot(mq_proteins$mq.protein.intensity, mq_proteins$mq.protein.ms.ms.count, log="xy")
#' abline(0,1, col=2)
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
#' evidence_txt <- system.file("samples/maxquant_txt/tiny2.zip",package = "prolfqua")
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
#' Generating mq file from proteinGroups.txt, peptide.txt and evidence.txt all level file incuding evidence
#' @param txt_directory or zip
#' @family MaxQuant
#' @export
#' @keywords internal
#' @examples
#'
#' txt_directory <- system.file("samples/maxquant_txt/tiny2.ZIP", package = "prolfqua")
#' allData <- tidyMQ_merged(txt_directory)
#'
tidyMQ_merged <- function(txt_directory){
  if (grepl("\\.zip$",tolower(txt_directory))) {
    proteins_txt <- read.csv(unz(txt_directory,"proteinGroups.txt"),
                             header = TRUE, sep = "\t", stringsAsFactors  =  FALSE)
    peptides_txt <- read.csv(unz(txt_directory,"peptides.txt"),
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    evidence_txt <- read.csv(unz(txt_directory,"evidence.txt"),
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }else{
    proteins_txt <- file.path(txt_directory, "proteinGroups.txt")
    peptides_txt <- file.path(txt_directory, "peptides.txt")
    evidence_txt <- file.path(txt_directory, "evidence.txt")
  }
  mq_proteins <- tidyMQ_ProteinGroups(proteins_txt)
  mq_peptides <- tidyMQ_Peptides(peptides_txt)
  mq_evidence <- tidyMQ_Evidence(evidence_txt)

  resProt_Pep <- dplyr::inner_join(mq_proteins,mq_peptides, by = c("protein.group.id", "raw.file"))
  resProt_Pep_Evidence <- dplyr::inner_join(resProt_Pep, mq_evidence, by = c("protein.group.id", "raw.file", "peptide.id"))
  return(resProt_Pep_Evidence)
}

#' Generating mq from proteinGroups, peptide.txt and modifiedPeptideSequence
#' @param txt_directory or zip
#'
#' @family MaxQuant
#' @export
#' @keywords internal
#' @examples
#'
#' #txt_directory <- system.file("samples/maxquant_txt/tiny2.ZIP", package = "prolfqua")
#' #allData <- tidyMQ_PeptideProtein(txt_directory)
#'
tidyMQ_PeptideProtein <- function(txt_directory, .all = FALSE){
  if (grepl("\\.zip$",tolower(txt_directory))) {
    proteins_txt <- read.csv(unz(txt_directory,"proteinGroups.txt"),
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    peptides_txt <- read.csv(unz(txt_directory,"peptides.txt"),
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    mod_spec_peptides_txt <- read.csv(unz(txt_directory,"modificationSpecificPeptides.txt"),
                                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  }else{
    proteins_txt <- file.path(txt_directory, "proteinGroups.txt")
    peptides_txt <- file.path(txt_directory, "peptides.txt")
    mod_spec_peptides_txt <- file.path(txt_directory, "modificationSpecificPeptides.txt")
  }
  mq_proteins <- tidyMQ_ProteinGroups(proteins_txt)
  mq_peptides <- tidyMQ_Peptides(peptides_txt)
  mq_modSpecPeptides <- tidyMQ_modificationSpecificPeptides(mod_spec_peptides_txt)

  resProt_Pep <- inner_join(mq_proteins,mq_peptides, by = c("protein.group.id", "raw.file"))

  if (.all) {
    return(list(resProt_Pep = resProt_Pep,
                mq_proteins = mq_proteins,
                mq_peptides = mq_peptides,
                mq_modSpecPeptides = mq_modSpecPeptides))
  }else{
    return(resProt_Pep)
  }
}

#' parse MQ modificationSpecificPeptides.txt
#' @param MQPeptides data.frame generated with read.csv("peptide.txt",sep = "\\t", stringsAsFactors = FALSE)
#' @family MaxQuant
#' @export
#' @keywords internal
#' @examples
#'
#' if(FALSE){
#' peptides_txt <- "d:/Dropbox/DataAnalysis/p2621_HumanAgeInteraction/data/721705/modificationSpecificPeptides.txt"
#' peptides_txt <- read.csv(peptides_txt,
#'  header = TRUE,
#'   stringsAsFactors = FALSE,
#'    sep = "\t")
#' head(peptides_txt)
#' MQPeptides <- peptides_txt
#' #View(MQPeptides)
#' mq_peptides <- tidyMQ_modificationSpecificPeptides(peptides_txt)
#'
#' head(mq_peptides)
#' }
tidyMQ_modificationSpecificPeptides <- function(MQPeptides){
  if (is.character(MQPeptides)) {
    if (grepl("\\.zip$",tolower(MQPeptides))) {
      message(MQPeptides)
      MQPeptides <- read.csv(unz(MQPeptides,"modificationSpecificPeptides.txt"),
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    }else{
      MQPeptides <- read.csv(MQPeptides, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    }
  }
  colnames(MQPeptides) <- tolower(colnames(MQPeptides))

  ### get mod columns
  meta <- dplyr::select(MQPeptides,
                        "mod.peptide.id" = "id",
                        "peptide.id",
                        "sequence",
                        "modifications",
                        "proteins",
                        "protein.group.id" = "protein.group.ids",
                        "mass",
                        "retention.time",
                        "mod.peptide.score" = "score",
                        "delta.score",
                        "pep",
                        dplyr::one_of("missed.cleavages"),
                        "unique.groups" = "unique..groups.",
                        "unique.proteins" = "unique..proteins.",
                        "reverse" = "reverse",
                        "potential.contaminant" = ends_with("contaminant")
                        )

  # add columns with modification to data.
  if ("missed.cleavages" %in% colnames(meta)) {
    stMODcol <- grep(pattern = "unique..proteins." ,colnames(MQPeptides)) + 1
    endMODcol <- grep(pattern = "missed.cleavages" ,colnames(MQPeptides)) - 1
    if (endMODcol > stMODcol) {
      mod_cols <- dplyr::select(MQPeptides, "id", stMODcol:endMODcol)
      mod_cols <- setNames(mod_cols , paste0("modification.", colnames(mod_cols)))
      #return(list(mod_cols = mod_cols, meta = meta))
      meta <- dplyr::inner_join(meta, mod_cols, by = c("mod.peptide.id" = "modification.id"))
    }
  }
  if (sum(grepl("site.ids",colnames(MQPeptides)))) {
    mod_cols2 <- dplyr::select(MQPeptides, "id", dplyr::ends_with("site.ids"))
    mod_cols2 <- setNames(mod_cols2 , paste0("site.ids.", colnames(mod_cols2)))
    meta <- dplyr::inner_join(meta, mod_cols2, by = c("mod.peptide.id" = "site.ids.id"))
  }
  sc <- sym("potential.contaminant")
  meta <- meta |>  dplyr::mutate(!!"potential.contaminant" := case_when( !!sc == "" ~ FALSE, !!sc == "+" ~ TRUE)) |>
    dplyr::mutate(!!"unique.groups" := case_when( !!sym("unique.groups") == "yes" ~ TRUE,
                                                  !!sym("unique.groups") == "no" ~ FALSE)) |>
    dplyr::mutate(!!"unique.proteins" := case_when( !!sym("unique.proteins") == "yes" ~ TRUE,
                                                    !!sym("unique.proteins") == "no" ~ FALSE)) |>
    dplyr::mutate(!!"reverse" := case_when( !!sym("reverse") == "+" ~ TRUE,
                                            !!sym("reverse") == "" ~ FALSE))


  pint <- dplyr::select(MQPeptides,"mod.peptide.id" = "id", starts_with("intensity."))
  PepIntensities <- pint |>
    tidyr::gather(key = "raw.file", value = "mod.peptide.intensity", starts_with("intensity.")) |>
    dplyr::mutate(raw.file = gsub("intensity.","",.data$raw.file))

  idtype <- dplyr::select(MQPeptides, "mod.peptide.id" = "id", starts_with("identification.type."))
  if (ncol(idtype) > 1) { # if only one file no id type is provided
    PepIDType <- idtype |>
      tidyr::gather(key = "raw.file", value = "mod.id.type", starts_with("identification.type.")) |>
      dplyr::mutate(raw.file = gsub("identification.type.","",.data$raw.file))
    PepIntensities <- dplyr::inner_join(PepIntensities,PepIDType, by = c("mod.peptide.id", "raw.file" ))
  }else{
    PepIntensities$id.type <- "By MS/MS"
  }
  xx <- dplyr::inner_join(meta , PepIntensities, by = "mod.peptide.id")
  xx$proteotypic <- !grepl(";",xx$protein.group.id)
  xx <- xx |> separate_rows(.data$protein.group.id, sep = ";",convert  = TRUE)
  return(xx)
}

#' parse MQ peptides.txt
#' @param MQPeptides data.frame generated with read.csv("peptide.txt",sep = "\\t", stringsAsFactors = FALSE)
#' @family MaxQuant
#' @export
#' @keywords internal
#' @examples
#'
#'
#' peptide_txt <- system.file("samples/maxquant_txt/tiny2.ZIP",package = "prolfqua")
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

#' parse MQ allPeptides.txt
#' @param MQPeptides data.frame generated with read.csv("peptide.txt",sep = "\\t", stringsAsFactors = FALSE)
#' @family MaxQuant
#' @export
#' @keywords internal
#' @examples
#' if(FALSE){
#' peptides_txt <- "c:/Users/wewol/Dropbox/DataAnalysis/p2621_HumanAgeInteraction/data/721705/allPeptides.txt"
#' peptides_txt <- read.csv(peptides_txt, header = TRUE, stringsAsFactors = FALSE, sep="\t")
#' MQPeptides <- peptides_txt
#' head(MQPeptides)
#' mq_peptides <- tidyMQ_allPeptides(peptides_txt)
#' dim(mq_peptides)
#' mq_peptides |> dplyr::filter(ms.ms.count != 0) -> idid
#' dim(idid)
#' mq_peptides |> dplyr::filter(sequence != "") -> idid
#' head(idid)
#' unique(idid$modified.sequence)
#' head(mq_peptides)
#'}
tidyMQ_allPeptides <- function(MQPeptides){
  if (is.character(MQPeptides)) {
    MQPeptides <- read.csv(MQPeptides, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
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
                      "peptide.score" = "score",
                      "intensity",
                      "ms.ms.count") |>
    dplyr::mutate(sequence = str_trim(sequence), modified.sequence = str_trim(.data$modified.sequence))
  return(xx)
}

#' convert modification specific to peptide level
#' aggregates mod.peptide.intensity, takes min of pep and max of peptide.score
#' use if you want to disable MQ default precursor filter
#' (which is not to use modified peptide sequences for peptide quantification)
#'
#' @family MaxQuant
#' @param mq_modSpecPeptides ouput of tidyMQ_modificationSpecificPeptides
#' @param mq_peptides  output of tidyMQ_Peptides
#' @export
#' @keywords internal
#'
tidyMQ_from_modSpecific_to_peptide <- function(mq_modSpecPeptides, mq_peptides) {
  warning("use only if you want to disable MQ default precursor filter")
  mq_modSpecPeptides <- mq_modSpecPeptides |> dplyr::filter(.data$unique.groups)
  relevantColumns <- setdiff(colnames(mq_peptides) , c("leading.razor.protein","id.type"))

  xx <- mq_modSpecPeptides |>
    dplyr::group_by(.data$peptide.id, .data$raw.file ) |>
    dplyr::mutate(peptide.intensity = sum(.data$mod.peptide.intensity, na.rm = TRUE),
                  pep = min(.data$pep, na.rm = TRUE),
                  peptide.score = max(.data$mod.peptide.score, na.rm = TRUE)) |>  dplyr::ungroup()

  dimcheck <- mq_modSpecPeptides |> dplyr::select(.data$peptide.id, .data$raw.file ) |>
    dplyr::distinct() |> nrow()

  peptides <- xx |> dplyr::select( one_of(relevantColumns) ) |>
    dplyr::distinct()
  stopifnot( dimcheck == nrow(peptides) )
  return(peptides)
}


#' Read Sites.txt and convert into longish format
#' @param dPat sitespecific file read by read.csv()
#'
#' @family MaxQuant
#' @export
#' @keywords internal
tidyMQ_from_Sites <- function(pDat){
  colnames(pDat) <- tolower(colnames(pDat) )
  # Cut PhosphoSTY down to useful columns

  annotation <-  pDat |> dplyr::select("localization.prob",
                                        "protein",
                                        "leading.proteins",
                                        "positions.within.proteins",
                                        "protein.group.ids" = "protein.group.ids",
                                        "pep",
                                        "score",
                                        "delta.score",
                                        "position.in.peptide",
                                        "charge",
                                        "amino.acid",
                                        "peptide.window.coverage",
                                        "sequence.window",
                                        "modification.window",
                                        "reverse",
                                        "potential.contaminant",
                                        "mod.peptide.ids" = "mod..peptide.ids",
                                        "site.id" = "id",
                                        number.of.sites = starts_with("number.of."),
                                        peptide.sequence.prob = ends_with("..probabilities"),
                                        peptide.sequence.scorediff = ends_with(".score.diffs"))
  #contains("___"))

  annotation <- mutate(annotation, potential.contaminant = case_when(potential.contaminant == "+" ~ TRUE, TRUE ~ FALSE))
  annotation <- mutate(annotation, reverse = case_when(reverse == "+" ~ TRUE, TRUE ~ FALSE))
  annotation <- mutate(annotation, stripped.sequence = gsub("\\(.+\\)","", .data$peptide.sequence.prob))

  intensities <- pDat |> dplyr::select("site.id" = "id", contains("___") )
  longish <- intensities |> gather(key = "raw.file", value = "intensity", contains("___"))
  longish <- longish |> separate(.data$raw.file, into = c("raw.file", "multiplicity"), sep="___")
  longish <- longish |> mutate(raw.file = gsub("intensity\\.", "", .data$raw.file))

  res <- inner_join(longish, annotation, by = "site.id")
  return(list(data = res, annotation = annotation))
}




#' Create retention time plot for MQ
#' @param MQtxtfolder txt folder
#' @param qc_path output folder
#' @family MaxQuant
#' @keywords internal
#' @export
plot_MQ_intensity_vs_retention_time <- function(MQtxtfolder, qc_path ) {

  mod_peptides_available <- "modificationSpecificPeptides.txt" %in% unzip(MQtxtfolder, list = TRUE)$Name

  res <- if (mod_peptides_available) {
    mqData <- tidyMQ_modificationSpecificPeptides(MQtxtfolder)

    {# create visualization for modified peptide sequence
      height <- length(unique(mqData$raw.file))/2 * 300
      png(file.path(qc_path, "retention_time_plot.png"), height = height, width = 1200)
      {
        resPepProtVis <- mqData |> dplyr::filter(.data$mod.peptide.intensity > 4)
        tmp <- ggplot(resPepProtVis, aes(x = .data$retention.time, y = log2(.data$mod.peptide.intensity))) +
          geom_point(alpha = 1/20, size = 0.3) +
          facet_wrap(~raw.file, ncol = 2)

        print(tmp)
      }
      dev.off()
    }
    TRUE
  }else{
    FALSE
  }

  return(res)
}





