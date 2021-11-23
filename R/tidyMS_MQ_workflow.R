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
#' library(prolfqua)
#' library(tidyverse)
#'
#' istar <- prolfqua_data('data_ionstar')$Pep()
#' istar_data <- istar$data %>% dplyr::filter(protein_Id %in% sample(protein_Id, 100))
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
#' library(prolfqua)
#' library(tidyverse)
#'
#' istar <-prolfqua_data('data_ionstar')$Pep()
#' istar_data <- istar$data %>% dplyr::filter(protein_Id %in% sample(protein_Id, 100))
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
#' @param startdata table in long format
#' @param atable AnalysisTableAnnotation annotate startdata table
#' @param GRP2 list with named arguments i.e. Contrasts, projectID, projectName, workunitID, nrPeptides, log2FCthreshold, FDRthreshold
#' @param protein_annot column with portein desciription e.g. (fasta header)
#' @param outpath directory to write results too.
#' @param revpattern default "REV_"
#' @param contpattern default "^zz|^CON__"
#' @param remove do you want to remove contaminants.
#' @export
#'
make2grpReport <- function(startdata,
                             atable,
                             GRP2,
                             protein_annot = "Description",
                             outpath = ".",
                             revpattern = "^REV_",
                             contpattern = "^zz|^CON__",
                             remove = FALSE){


  proteinID <- atable$hkeysDepth()
  config <- prolfqua::AnalysisConfiguration$new(atable)

  annotProtein <- function(startdata , Accession, Description, revpattern = "^REV_", contpattern = "zzY-FGCZ"){
    GRP2 <- list()
    distinctprotid <- startdata %>% select(pID = !!sym(Accession), fasta.headers = {{Description}}) %>% distinct()
    distinctprotid <- distinctprotid %>% mutate(proteinAnnot = case_when(grepl(revpattern,pID) ~ "REV",
                                                                         grepl(contpattern,pID) ~ "CON",
                                                                         TRUE ~ "FW"))
    GRP2$percentOfContaminants <-  round(mean(distinctprotid$proteinAnnot == "CON") * 100, digits = 2)
    GRP2$percentOfFalsePositives <- round(mean(distinctprotid$proteinAnnot == "REV") * 100, digits = 2)
    GRP2$totalNrOfProteins <- sum(table(distinctprotid$proteinAnnot))
    GRP2$NrOfProteinsNoDecoys <- sum(distinctprotid$proteinAnnot == "FW")
    return(list(stats = GRP2, distinctprotid = distinctprotid))
  }

  res <- annotProtein(startdata, Accession = atable$hierarchy[[1]], !!sym(protein_annot), revpattern = revpattern, contpattern = contpattern)
  GRP2 <- c(GRP2, res$stats)

  if (remove) {
    distinctprotid <- dplyr::filter(res$distinctprotid, .data$proteinAnnot == "FW")
    startdata <- startdata %>% filter(!!sym(atable$hierarchy[[1]]) %in% distinctprotid$pID)
  } else {
    distinctprotid <- res$distinctprotid
  }


  ############################## Create configuration For MQ ####
  adata <- setup_analysis(startdata, config)

  ##################### Preprocess intensities ###################################

  lfqdata <- LFQData$new(adata, config)
  lfqdata$remove_small_intensities()



  ### Do some type of data normalization (or do not)
  lt <- lfqdata$get_Transformer()
  transformed <- lt$log2()$robscale()$lfq
  if (length(transformed$config$table$hierarchyKeys()) > transformed$config$table$hierarchyDepth) {
    transformed$filter_proteins_by_peptide_count()
    aggregator <- transformed$get_Aggregator()
    aggregator$medpolish()
    transformed <- aggregator$lfq_agg
  }

  GRP2$lfqData <- lfqdata
  GRP2$transformedlfqData <- transformed
  #GRP2$transformedlfqData$data$transformedIntensity <- GRP2$transformedlfqData$data$transformedIntensity * 3

  ################## Run Modelling ###############


  formula_Condition <-  strategy_lm(paste0(transformed$config$table$getWorkIntensity(), " ~ ",
                                           transformed$config$table$fkeysDepth()))
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
  conMI <- ContrastsModerated$new(mC, modelName = "Imputed_Data")

  res <- prolfqua::addContrastResults(conrM, conMI)


  GRP2$contrResult <- res$merged$get_contrasts()
  GRP2$contrMerged <- res$merged$get_Plotter()
  GRP2$contrMerged$fcthresh = GRP2$log2FCthreshold
  GRP2$contrMerged$volcano_spec[[1]]$thresh = GRP2$FDRthreshold

  GRP2$contrMore <- res$more$get_Plotter()

  top20 <- GRP2$contrResult %>% dplyr::select( !!sym(proteinID ),log2FC = .data$estimate,.data$conf.low,.data$conf.high, .data$FDR ) %>%
    arrange(.data$FDR) %>%
    head(20)
  GRP2$top20 <- top20
  #knitr::kable(top20, caption = "Top 20 proteins sorted by smallest Q Value (adj.P.Val). The effectSize column is the log2 FC of condition vs reference.")

  GRP2$top20confint <- ggplot(top20, aes(x = !!sym(proteinID), y = .data$log2FC,
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
  if (nrow(xx) > 0) {
    xx <- xx %>% arrange(.data$estimate)
    GRP2$noPvalEstimate <- ggplot2::ggplot(xx ,aes(x = stats::reorder(!!sym(proteinID), .data$estimate), y = .data$estimate)) +
      ggplot2::geom_bar(stat = "identity") + coord_flip()
    missing <- GRP2$transformedlfqData$get_copy()
    missing$complete_cases()
    missingID <- xx[[ proteinID ]]
    missing$data <- missing$data %>% dplyr::filter(!!sym(proteinID ) %in% missingID)
    missing$get_Plotter()$raster()
  }


  ### -----

  wr <- GRP2$lfqData$get_Writer()
  tmp <- wr$get_wide()
  tmp$data <- inner_join(distinctprotid, tmp$data, by = c(pID = proteinID) )
  tmp2 <- GRP2$transformedlfqData$get_Writer()$get_wide()

  names(tmp2) <- paste0(names(tmp2), ".normalized")
  tmp2$data.normalized <- inner_join(distinctprotid, tmp2$data.normalized, by = c(pID = proteinID))
  res <- inner_join(distinctprotid, GRP2$contrResult, by = c(pID = proteinID))

  dir.create(outpath)

  writexl::write_xlsx(c(tmp, tmp2,  contrasts = list(res)),
                      path = file.path(outpath,"AnalysisResults.xlsx"))


  prolfqua::copy_2grp_markdown()
  rmarkdown::render("_GRP2Analysis.Rmd",
                    params = list(grp = GRP2) ,
                    output_format = bookdown::html_document2(toc = TRUE,toc_float = TRUE))
  file.copy("_GRP2Analysis.html", file.path(outpath, "_GRP2Analysis.html"), overwrite = TRUE)
}
