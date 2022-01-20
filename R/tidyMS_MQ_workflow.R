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
#' @param GRP2 list with named arguments i.e. Contrasts, projectID, projectName, workunitID, nrPeptides, log2FCthreshold, FDRthreshold
#' @param protein_annot column with portein desciription e.g. (fasta header)
#' @param revpattern default "REV_"
#' @param contpattern default "^zz|^CON__"
#' @param remove do you want to remove contaminants.
#' @export
#' @family workflow
#' @examples
#'
#'
#'
#' istar <- prolfqua_data('data_ionstar')$filtered()
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#'
#' GRP2 <- list()
#' GRP2$projectID <- "3765"
#' GRP2$projectName <- "Order_26863"
#' GRP2$workunitID <- "2057368.zip"
#'
#' #at least 2 peptides per protein
#' GRP2$nrPeptides <- 2
#'
#' # Set FC to >= |2| and FRD to 0.1
#' GRP2$log2FCthreshold <- 1
#' GRP2$FDRthreshold <- 0.1
#' GRP2$Contrasts <- c(avsb = "dilution.b - dilution.a")
#' data$Description <-"AAAAA"
#'
#'
#'
#'  atab <- AnalysisTableAnnotation$new()
#'
#' atab$fileName = "raw.file"
#' atab$hierarchy["protein_Id"] = "protein_Id"
#' atab$hierarchy["peptide_Id"] = "peptide_Id"
#'
#' atab$factors["dilution."] = "dilution."
#' atab$setWorkIntensity("peptide.intensity")
#' atab$isotopeLabel = "isotope"
#' config <- prolfqua::AnalysisConfiguration$new(atab)
#'
#' protein_annot = "Description"
#' #debug(make2grpReport)
#' grp <- make2grpReport(data,atab, GRP2, NULL)
#'
#' \dontrun{
#' write_2GRP(grp,tempdir())
#' render_2GRP(grp, tempdir())
#' }
make2grpReport <- function(startdata,
                           atable,
                           GRP2,
                           protein_annot = "Description",
                           revpattern = "^REV_",
                           contpattern = "^zz|^CON__",
                           remove = FALSE) {

  proteinID <- atable$hkeysDepth()
  config <- prolfqua::AnalysisConfiguration$new(atable)


  adata <- setup_analysis(startdata, config)
  prot_annot <- select( startdata , c( atable$hierarchy[[atable$hkeysDepth()]], protein_annot)) |> distinct()
  prot_annot <- rename(prot_annot, !!atable$hkeysDepth() := (!!atable$hierarchy[[atable$hkeysDepth()]]))

  lfqdata <- LFQData$new(adata, config)
  lfqdata$remove_small_intensities()


  ### Do some type of data normalization (or do not)
  lt <- lfqdata$get_Transformer()
  transformed <- lt$log2()$robscale()$lfq
  if ( length(transformed$config$table$hierarchyKeys()) > transformed$config$table$hierarchyDepth ) {
    message("AGGREGATING PEPTIDE DATA!")
    transformed$filter_proteins_by_peptide_count()
    aggregator <- transformed$get_Aggregator()

    aggregator$medpolish()
    transformed <- aggregator$lfq_agg
  }


  protAnnot <- RowAnnotProtein$new(
    transformed,
    row_annot = prot_annot)

  allProt <- nrow( protAnnot$row_annot )
  GRP2$totalNrOfProteins <- allProt
  GRP2$percentOfContaminants <- round(protAnnot$annotateREV(revpattern)/allProt * 100 , digits = 2)
  GRP2$percentOfFalsePositives  <- round(protAnnot$annotateCON(contpattern)/allProt * 100 , digits = 2)
  GRP2$NrOfProteinsNoDecoys <- protAnnot$nr_clean()
  GRP2$rowAnnot <- protAnnot

  if (remove) {
    message("REMOVING: contaminants and reverse sequences")
    lfqdata <- lfqdata$get_subset(protAnnot$clean())
    transformed <- transformed$get_subset(protAnnot$clean())
  }


  GRP2$lfqData <- lfqdata
  GRP2$transformedlfqData <- transformed

  ################## Run Modelling ###############

  formula <- paste0(transformed$config$table$getWorkIntensity(), " ~ ",
         paste(transformed$config$table$factorKeys(), collapse = " + "))
  message("FORMULA :",  formula)
  GRP2$formula <- formula
  formula_Condition <-  strategy_lm(formula)
  # specify model definition
  modelName  <- "Model"

  mod <- prolfqua::build_model(
    transformed,
    formula_Condition,
    subject_Id = transformed$config$table$hierarchyKeys() )

  GRP2$models <- mod

  contr <- prolfqua::Contrasts$new(mod, GRP2$Contrasts)
  conrM <- ContrastsModerated$new(contr, modelName = "Linear_Model_Moderated")
  mC <- ContrastsSimpleImpute$new(lfqdata = transformed, contrasts = GRP2$Contrasts)
  conMI <- ContrastsModerated$new(mC, modelName = "Imputed_Condition_Mean")

  res <- prolfqua::addContrastResults(conrM, conMI)


  GRP2$contrResult <- res$merged$get_contrasts()
  GRP2$contrMerged <- res$merged$get_Plotter()
  GRP2$contrMerged$fcthresh = GRP2$log2FCthreshold
  GRP2$contrMerged$volcano_spec[[1]]$thresh = GRP2$FDRthreshold

  GRP2$contrMore <- res$more

  top20 <- GRP2$contrResult |>
    dplyr::select( !!sym(proteinID ),
                   diff = .data$diff,
                   .data$conf.low,
                   .data$conf.high,
                   .data$FDR ) |>
    arrange(.data$FDR) |>
    head(20)
  GRP2$top20 <- top20
  GRP2$top20confint <- ggplot(top20, aes(x = !!sym(proteinID), y = .data$diff,
                                         ymin = .data$conf.low, ymax = .data$conf.high)) +
    geom_hline( yintercept = 0, color = 'red' ) +
    geom_linerange() + geom_point() + coord_flip() + theme_minimal()


  protMore <- GRP2$transformedlfqData$get_copy()
  protMore$complete_cases()
  moreID <-  res$more$get_contrasts()[[proteinID]]
  protMore$data <-  dplyr::filter(protMore$data ,!!sym(proteinID ) %in% moreID)
  GRP2$imputedProteins <- protMore

  # Plot proteins without p-values

  xx <- res$more$contrast_result[rowSums(is.na(res$more$get_contrasts())) > 0,]
  if (nrow(xx) > 1) {
    xx <- xx |> arrange(.data$diff)
    GRP2$noPvalEstimate <- ggplot2::ggplot(
      xx ,
      aes(x = stats::reorder(!!sym(proteinID),
                             .data$diff),
          y = .data$diff)) +
      ggplot2::geom_bar(stat = "identity") + coord_flip()


    missing <- GRP2$transformedlfqData$get_copy()
    missing$complete_cases()
    missingID <- xx[[ proteinID ]]
    missing$data <- missing$data |> dplyr::filter(!!sym(proteinID ) %in% missingID)
    missing$get_Plotter()$raster()
  }

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
  rd <- GRP2$lfqData
  tr <- GRP2$transformedlfqData
  ra <- GRP2$rowAnnot
  wideraw <- inner_join(ra$row_annot, rd$to_wide()$data)
  widetr <- inner_join(ra$row_annot , tr$to_wide()$data )
  ctr <- inner_join(ra$row_annot , GRP2$contrResult )
  resultList <- list()
  resultList$annotation = tr$to_wide()$annot
  resultList$raw_data = wideraw
  resultList$transformed_data = widetr
  resultList$contrasts = ctr
  writexl::write_xlsx(resultList, path = file.path(outpath, paste0(xlsxname, ".xlsx")))
}

#' render 2GRP analysis report
#' @rdname make2grpReport
#' @param GRP2 return value of \code{\link{make2grpReport}}
#' @param outpath path to place output
#' @param htmlname name for html file
#' @export
#' @family workflow
render_2GRP <- function(GRP2, outpath, htmlname="Result2Grp"){
  prolfqua::copy_2grp_markdown()
  rmarkdown::render("_Grp2Analysis.Rmd",
                    params = list(grp = GRP2) ,
                    output_format = bookdown::html_document2(toc = TRUE,toc_float = TRUE))
  dir.create(outpath)
  file.copy("_Grp2Analysis.html", file.path(outpath, paste0(htmlname,".html")), overwrite = TRUE)
}
