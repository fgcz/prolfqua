# Helper functions -----

#' Keep only those proteins with 2 IDENTIFIED peptides
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @return list with data.frame (data) and name of new column (name)
#' @export
#' @keywords internal
#' @family preprocessing
#' @examples
#'
#'
#'
#' istar <- prolfqua_data('data_ionstar')$Pep()
#' istar_data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' filterPep <- prolfqua::filter_proteins_by_peptide_count( istar_data ,  istar$config )
#'  x <- prolfqua::summarize_hierarchy(filterPep$data , istar$config)
#' stopifnot(x$peptide_Id_n >= istar$config$parameter$min_peptides_protein)
#'
#'
filter_proteins_by_peptide_count <-
  function(pdata,
           config){

    # remove single hit wonders
    tmp <- nr_B_in_A(pdata,config)
    if (!is.null(tmp)) {
      res <- dplyr::filter(tmp$data, !!sym(tmp$name) >= config$parameter$min_peptides_protein )
      name <- tmp$name
    }else{
      res <- pdata
      name <- NULL
    }
    return(list(data = res, name = name))
  }



#' get the difference of two dataset where one is a subset of the other.
#' @param x data.frame
#' @param y data.frame
#' @param config AnlysisConfiguration
#' @return data.frame
#' @export
#' @keywords internal
#'
#' @examples
#'
#'
#'
#' istar <-prolfqua_data('data_ionstar')$Pep()
#' istar_data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' filterPep <- prolfqua:::filter_proteins_by_peptide_count( istar_data ,  istar$config )
#' tmp <- filter_difference(istar_data, filterPep$data, istar$config)
#' stopifnot(nrow(istar_data )  - nrow(filterPep$data) == nrow(tmp))
#' tmp <- filter_difference(filterPep$data, istar_data , istar$config)
#' stopifnot(nrow(istar_data )  - nrow(filterPep$data) == nrow(tmp))
#'
filter_difference <- function(x, y, config){
  if (nrow(y) > nrow(x)) {
    dplyr::anti_join(y, x, by = config$table$idVars())
  }else{
    dplyr::anti_join(x, y, by = config$table$idVars())
  }
}


#' Create 2 grp report in html and write data to xlsx table
#'
#' For use examples see run_scripts directory
#' @rdname make2grpReport
#' @param startdata table in long format
#' @param atable AnalysisTableAnnotation annotate startdata table
#' @param GRP2 list with named arguments i.e. Contrasts, projectID, projectName, workunitID, nrPeptides, Diffthreshold, FDRthreshold
#' @param protein_annot column with protein description e.g. (fasta header)
#' @param revpattern default "REV_"
#' @param contpattern default "^zz|^CON__"
#' @param remove do you want to remove contaminants default (TRUE)
#' @param transform which transformation to use to normalize the data, default robscale
#' @param aggregate method to use for peptide aggregation
#' @export
#' @family workflow
#' @examples
#'
#' istar <- prolfqua_data('data_ionstar')$filtered()
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' data$Description <-"AAAAA"
#'
#' GRP2 <- list()
#' GRP2$Bfabric <- list()
#' GRP2$Bfabric$projectID <- "3765"
#' GRP2$Bfabric$projectName <- "Order_26863"
#' GRP2$Bfabric$orderID <- "3765"
#'
#' GRP2$Bfabric$workunitID <- "2057368.zip"
#' GRP2$Bfabric$inputID <- "2057368.zip"
#' GRP2$Bfabric$inputURL <- "https://www.fgcz.ch"
#'
#' #at least 2 peptides per protein
#' GRP2$pop <- list()
#' GRP2$pop$transform <- "vsn"
#' GRP2$pop$aggregate <- "medpolish"
#' GRP2$pop$Diffthreshold <- 0.5
#' GRP2$pop$FDRthreshold <- 0.25
#' GRP2$pop$Contrasts <- c(b_vs_a = "dilution.b - dilution.a")
#'
#' GRP2$Software <- "MaxQuant"
#'
#' data <- dplyr::filter(data, dilution. == "a" |  dilution. == "b")
#' atab <- AnalysisTableAnnotation$new()
#'
#' atab$fileName = "raw.file"
#' atab$hierarchy["protein_Id"] = "protein_Id"
#' atab$hierarchy["peptide_Id"] = "peptide_Id"
#' atab$factors["dilution."] = "dilution."
#' atab$setWorkIntensity("peptide.intensity")
#' atab$isotopeLabel = "isotope"
#' config <- prolfqua::AnalysisConfiguration$new(atab)
#'
#' protein_annot = "Description"
#' grp <- make2grpReport(data, atab, GRP2)
#'
#' \dontrun{
#'
#' render_2GRP(grp, ".")
#' render_2GRP(grp, "." ,word = TRUE)
#' write_2GRP(grp,".")
#' }
make2grpReport <- function(startdata,
                           atable,
                           GRP2,
                           protein_annot = "Description",
                           revpattern = "^REV_",
                           contpattern = "^zz|^CON__",
                           remove = FALSE
                           ) {
  GRP2$pop$nrPeptides <- 2

  # Preprocess Data
  config <- prolfqua::AnalysisConfiguration$new(atable)
  adata <- setup_analysis(startdata, config)
  proteinID <- atable$hkeysDepth()

  prot_annot <- select(startdata , c( atable$hierarchy[[proteinID]], protein_annot)) |> distinct()
  prot_annot <- rename(prot_annot, !!proteinID := (!!atable$hierarchy[[proteinID]]))

  lfqdata <- LFQData$new(adata, config)
  lfqdata$remove_small_intensities()


  ### Do some type of data normalization (or do not)
  lt <- lfqdata$get_Transformer()
  if (GRP2$pop$transform == "robscale") {
    transformed <- lt$log2()$robscale()$lfq
  } else if (GRP2$pop$transform == "vsn") {
    transformed <- lt$intensity_matrix( .func = vsn::justvsn)$lfq
  } else if (GRP2$pop$transform == "none") {
    transformed <- lt$log2()$lfq
  } else {
    logger::log_warn("no such transformaton : {GRP2$pop$transform}")
  }
  logger::log_info("data transformed : {GRP2$pop$transform}.")


  ### Aggregate peptides to proteins
  if ( length(transformed$config$table$hierarchyKeys()) > transformed$config$table$hierarchyDepth ) {
    message("AGGREGATING PEPTIDE DATA!")
    transformed$filter_proteins_by_peptide_count()
    aggregator <- transformed$get_Aggregator()
    if (GRP2$pop$aggregate == "medpolish") {
      aggregator$medpolish()
    } else if (GRP2$pop$aggregate == "topN") {
      aggregator$sum_topN()
    } else if (GRP2$pop$aggregate == "lmrob") {
      aggregator$lmrob()
    } else {
      logger::log_warn("no such aggregator {GRP2$pop$aggregate}.")
    }
    logger::log_info("data aggregated: {GRP2$pop$aggregate}.")

    transformed <- aggregator$lfq_agg
  }

  ## count contaminants.
  protAnnot <- RowAnnotProtein$new(
    transformed,
    row_annot = prot_annot)

  allProt <- nrow( protAnnot$row_annot )
  GRP2$RES <- list()
  GRP2$RES$Summary <- data.frame(
    totalNrOfProteins = allProt,
    percentOfContaminants = round(protAnnot$annotateCON(contpattern)/allProt * 100 , digits = 2),
    percentOfFalsePositives  = round(protAnnot$annotateREV(revpattern)/allProt * 100 , digits = 2),
    NrOfProteinsNoDecoys = protAnnot$nr_clean()
  )
  GRP2$RES$rowAnnot <- protAnnot

  if (remove) {
    lfqdata <- lfqdata$get_subset(protAnnot$clean())
    transformed <- transformed$get_subset(protAnnot$clean())
    logger::log_info("removing contaminants and reverse sequences with patterns: {contpattern} {revpattern}")
  }

  GRP2$RES$lfqData <- lfqdata
  GRP2$RES$transformedlfqData <- transformed

  ################## Run Modelling ###############

  formula <- paste0(transformed$config$table$getWorkIntensity(), " ~ ",
         paste(transformed$config$table$factorKeys(), collapse = " + "))
  message("FORMULA :",  formula)
  GRP2$RES$formula <- formula
  formula_Condition <-  strategy_lm(formula)
  # specify model definition
  modelName  <- "Model"

  mod <- prolfqua::build_model(
    transformed,
    formula_Condition,
    subject_Id = transformed$config$table$hierarchyKeys() )

  logger::log_info("fitted model with formula : {formula}")
  GRP2$RES$models <- mod

  contr <- prolfqua::Contrasts$new(mod, GRP2$pop$Contrasts)
  conrM <- ContrastsModerated$new(contr, modelName = "Linear_Model_Moderated")
  mC <- ContrastsSimpleImpute$new(lfqdata = transformed, contrasts = GRP2$pop$Contrasts)
  conMI <- ContrastsModerated$new(mC, modelName = "Imputed_Mean")

  res <- prolfqua::addContrastResults(conrM, conMI)
  GRP2$RES$contrMerged <- res$merged
  GRP2$RES$contrMore <- res$more
  return(GRP2)
}


#' write 2 grp results
#' @rdname make2grpReport
#' @param GRP2 return value of \code{\link{make2grpReport}}
#' @param outpath path to place output
#' @param xlsxname file name for xlsx
#' @export
#' @family workflow
#'
write_2GRP <- function(GRP2, outpath, xlsxname = "AnalysisResults"){
  dir.create(outpath)
  rd <- GRP2$RES$lfqData
  tr <- GRP2$RES$transformedlfqData
  ra <- GRP2$RES$rowAnnot
  formula <- data.frame(formula = GRP2$RES$formula, contrast_name = names(GRP2$pop$Contrasts), contrast = GRP2$pop$Contrasts)
  wideraw <- inner_join(ra$row_annot, rd$to_wide()$data)
  widetr <- inner_join(ra$row_annot , tr$to_wide()$data )
  ctr <- inner_join(ra$row_annot , GRP2$RES$contrMerged$get_contrasts())
  resultList <- list()
  resultList$annotation = tr$to_wide()$annot
  resultList$raw_data = wideraw
  resultList$transformed_data = widetr
  resultList$contrasts = ctr
  resultList$formula = formula
  resultList$summary = GRP2$RES$Summary
  writexl::write_xlsx(resultList, path = file.path(outpath, paste0(xlsxname, ".xlsx")))
}

#' render 2GRP analysis report
#' @rdname make2grpReport
#' @param GRP2 return value of \code{\link{make2grpReport}}
#' @param outpath path to place output
#' @param htmlname name for html file
#' @param word default FALSE, if true create word document.
#' @export
#' @family workflow
render_2GRP <- function(GRP2, outpath, htmlname="Result2Grp", word = FALSE){
  prolfqua::copy_2grp_markdown()
  dir.create(outpath)

  rmarkdown::render(
    "_Grp2Analysis.Rmd",
    params = list(grp = GRP2) ,
    output_format = if(word){
      bookdown::word_document2(toc = TRUE, toc_float = TRUE) } else {
        bookdown::html_document2(toc = TRUE, toc_float = TRUE)
      }
  )
  fname <- paste0("_Grp2Analysis", if(word) {".docx"} else {".html"})
  if (file.copy(fname, file.path(outpath, paste0(htmlname,if(word) {".docx"} else {".html"})), overwrite = TRUE)) {
    file.remove(fname)
  }
}

