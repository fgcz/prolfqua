#' Extracts uniprot ID
#'
#' @export
#' @keywords internal
#' @family utilities
#' @return data.frame
#' @param obj data.frame
#' @param idcol columns to extract uniprot id's from
#'
#' @examples
#'
#' bb <- old2new(prolfqua_data('data_ionstar')$filtered())
#' colnames(bb$data)
#' tmp <- prolfqua::separate_hierarchy(bb$data, bb$config)
#' tmp$UniprotID <- NULL
#' tmp <- get_UniprotID_from_fasta_header(tmp, idcolumn = "top_protein")
#' stopifnot("UniprotID" %in%  colnames(tmp))
get_UniprotID_from_fasta_header <- function(df, idcolumn = "protein_Id")
{
  map <- df |> dplyr::select(idcolumn) |> distinct() |>
    dplyr::filter(grepl(pattern = "^sp|^tr", !!sym(idcolumn))) |>
    tidyr::separate(col = !!sym(idcolumn),
                    sep = "_",
                    into = c("begin",
                             "end"), remove = FALSE) |>
    tidyr::separate(col = !!sym("begin"),
                    sep = "\\|",
                    into = c("prefix", "UniprotID", "Symbol")) |>
    dplyr::select(-!!sym("prefix"), -!!sym("Symbol"), -!!sym("end"))
  res <- dplyr::right_join(map, df, by = idcolumn)
  return(res)
}


#' Removes rows with more than thresh NA's from matrix
#' @export
#' @keywords internal
#' @family utilities
#' @return matrix
#' @param obj matrix or dataframe
#' @param thresh - maximum number of NA's / row - if more the row will be removed
#' @examples
#'
#' obj = matrix(rnorm(10*10),ncol=10)
#' dim(obj)
#' obj[3,3] = NA
#' x1 = remove_NA_rows(obj, thresh=0)
#' stopifnot(all(c(9,10)==dim(x1)))
#' x2 = remove_NA_rows(obj, thresh=1)
#' stopifnot(all(c(10,10)==dim(x2)))
remove_NA_rows <- function(obj, thresh=0 )
{
  x <- apply(obj,1,function(x){sum(is.na(x))})
  obj <- obj[!(x > thresh),]
}

#' splits names and creates a matrix
#' @export
#' @keywords internal
#' @family utilities
#' @param names vector with names
#' @param split patter to use to split
#' @return matrix
#'
#' @examples
#' dat = c("bla_ra0/2_run0","bla_ra1/2_run0","bla_ra2/2_run0")
#' names_to_matrix(dat,split="\\_|\\/")
names_to_matrix <- function(names,split="\\||\\_")
{
  cnamessplit <- stringr::str_split(as.character(names),pattern = split)
  protnam <- do.call("rbind",cnamessplit)
  return(protnam)
}




#' plot volcano given multiple contrasts
#' @param misspX data in long format
#' @param effect column containing effect sizes
#' @param p.value column containing p-values, q.values etc
#' @param contrast column with contrast
#' @param colour colouring of points
#' @param xintercept fc thresholds
#' @param sigthreshold sigthreshold threshold
#' @param label column containing labels
#' @param size controls size of text
#' @param segment.size controls size of lines
#' @param segment.alpha controls visibility of lines
#' @param ablines adds ablines horizontal and vertical
#' @param scales parameter to ggplot2::facet_wrap
#' @export
#' @keywords internal
#' @family utilities
#' @examples
#'
#'
#' show <- prolfqua_data('data_multigroupFC') |>
#'    dplyr::filter(Condition  %in% c("TTB7_38h_test - TTB7_16h_test","TTB7_96h_test - TTB7_96h_sys") )
#' prolfqua::multigroup_volcano( show,
#' effect="logFC",
#' significance = "adj.P.Val",
#' contrast="Condition",
#' colour=NULL,
#' label="Name",
#' maxNrOfSignificantText = 300)
multigroup_volcano <- function(.data,
                              effect = "fc",
                              significance = "p.adjust",
                              contrast = "condition",
                              colour = "colour",
                              xintercept = c(-2,2),
                              yintercept = 0.05,
                              label = NULL,
                              size = 1,
                              segment.size = 0.3,
                              segment.alpha = 0.3,
                              scales = "fixed",
                              maxNrOfSignificantText = 20)
{
  misspX <- tidyr::unite(.data, "label", label)

  p <- .multigroup_volcano(misspX,
                          effect = effect,
                          significance = significance,
                          contrast = contrast,
                          colour = colour,
                          xintercept = xintercept,
                          yintercept = yintercept,
                          text = "label",
                          scales = scales)
  colname = paste("-log10(", significance , ")" , sep = "")

  if (!is.null(label)) {
    effectX <- misspX[,effect]
    typeX <- misspX[,significance]
    subsetData <- subset(misspX, (effectX < xintercept[1] | xintercept[2] < effectX) & typeX < yintercept ) |>
      head(n = maxNrOfSignificantText)
    if (nrow(subsetData) > 0) {
      p <- p + ggrepel::geom_text_repel(data = subsetData,
                                        aes_string(effect , colname , label = "label"),
                                        size = size,
                                        segment.size = segment.size,
                                        segment.alpha = segment.alpha)
    }
  }
  return(p)
}

.multigroup_volcano <- function(data,
                               effect = "fc",
                               significance = "p.adjust",
                               contrast = "contrast",
                               colour = "colour",
                               xintercept = c(-2,2),
                               yintercept = 0.05,
                               text = NULL,

                               scales = "fixed")
{
  colname = paste("-log10(", significance, ")", sep = "")
  p <- ggplot(
    data,
    aes_string(x = effect, y = colname, color = colour, text = text)) +
    geom_point(alpha = 0.5)
  p <- p + scale_colour_manual(
    values = c("black", "green",
               "blue", "red") )
  p <- p + facet_wrap(
    as.formula(paste("~", contrast)),
    scales = scales) + labs(y = colname)
  log2FC <- effect
  if ( !is.null(xintercept) ) {
    p <- p + geom_vline(
      xintercept = xintercept,
      linetype = "dashed",
      colour = "red")
  }
  if ( !is.null(yintercept) ) {
    p_value <- paste0("-log10(",yintercept,")")
    p <- p + geom_hline(
      yintercept = -log10(yintercept),
      linetype = "dashed",
      color = "blue",
      show.legend = FALSE)
  }

  return(p)
}



#' matrix or data.frame to tibble (taken from tidyquant)
#'
#' @param x a matrix
#' @param preserve_row_names give name to rownames column, if NULL discard rownames
#' @param ... further parameters passed to as_tibble
#' @export
#' @family utilities
#' @keywords internal
#' @examples
#' x <- matrix(rnorm(20), ncol=4)
#' rownames(x) <- LETTERS[seq_len(nrow(x))]
#' matrix_to_tibble(x)
#' !(is.matrix(x) || is.data.frame(x))
#'
matrix_to_tibble <- function(x, preserve_row_names = "row.names", ... )
{
  if (!(is.matrix(x) || is.data.frame(x))) stop("Error: `x` is not a matrix or data.frame object.")
  if (!is.null(preserve_row_names)) {
    row.names <- rownames(x)
    if (!is.null(row.names)  ) {
      dplyr::bind_cols(
        tibble::tibble(!! preserve_row_names := row.names),
        tibble::as_tibble(x, ...)
      )

    } else {

      warning("Warning: No row names to preserve. ",
              "Object otherwise converted to tibble successfully.")
      tibble::as_tibble(x, ...)
    }

  } else {

    tibble::as_tibble(x, ...)

  }
}


# @TODO think of making it public.
#' Compute jackknive resampling on matrix
#' @param xdata matrix
#' @param .method method i.e. cor, parameters
#' @param ... further parameters to .method
#' @return list with all jackknife matrices
#' @export
#' @return list of matrices
#' @family transitioncorrlation
#' @keywords internal
#' @examples
#' xx <- matrix(rnorm(20), ncol=4)
#' cortest <- function(x){print(dim(x));cor(x)}
#' jackknife(xx, cortest)
#' tmp <- jackknife(xx, cor, use="pairwise.complete.obs", method="pearson")
#' stopifnot(length(tmp) == 3)
#' stopifnot(length(tmp[[2]]) == nrow(xx))
jackknife <- function(xdata, .method, ... ) {
  x <- seq_len(nrow(xdata))
  call <- match.call()
  n <- length(x)
  u <- vector( "list", length = n )
  for (i in seq_len(n)) {
    tmp <- xdata[x[-i],]
    u[[i]] <- .method(tmp, ...)
  }
  names(u) <- seq_len(n)
  thetahat <- .method(xdata, ...)
  invisible(list(thetahat = thetahat, jack.values = u, call = call ))
}

#' Compute correlation matrix with jack
#' @seealso jackknife
#' @param dataX data.frame with transition intensities per peptide
#' @param distmethod dist or correlation method working with matrix i.e. cor
#' @param ... further parameters to method
#' @export
#' @family transitioncorrlation
#' @keywords internal
#' @return summarizes results producced with jackknife
#' @examples
#' dataX <- matrix(rnorm(20), ncol=4)
#' rownames(dataX)<- paste("R",seq_len(nrow(dataX)),sep="")
#' colnames(dataX)<- paste("C",seq_len(ncol(dataX)),sep="")
#' tmp <- jackknife(dataX, cor, use="pairwise.complete.obs", method="pearson")
#' res <- jackknife_matrix(dataX, cor)
#' res
#' stopifnot(dim(res) == c(4,4))
#' res <- jackknife_matrix(dataX, cor, method="spearman")
#' stopifnot(dim(res) == c(4,4))

jackknife_matrix <- function(dataX, distmethod , ... ){
  if (is.null(colnames(dataX))) {
    colnames(dataX) <- paste("C", seq_len(ncol(dataX)), sep = "")
  }
  if (is.null(rownames(dataX))) {
    rownames(dataX) <- paste("R", seq_len(nrow(dataX)), sep = "")
  }

  if (nrow(dataX) > 1 & ncol(dataX) > 1) {
    tmp <- jackknife( dataX, distmethod, ... )
    x <- purrr::map_df(tmp$jack.values, prolfqua::matrix_to_tibble)
    dd <- tidyr::gather(x, "col.names" , "correlation" , 2:ncol(x))
    ddd <- dd |>
      group_by(UQ(sym("row.names")), UQ(sym("col.names"))) |>
      summarize_at(c("jcor" = "correlation"), function(x){max(x, na.rm = TRUE)})

    dddd <- tidyr::spread(ddd, UQ(sym("col.names")), UQ(sym("jcor"))  )
    dddd <- as.data.frame(dddd)
    rownames(dddd) <- dddd$row.names
    dddd <- dddd[,-1]
    return(dddd)
  }else{
    message("Could not compute correlation, nr rows : " , nrow(dataX) )
    return(NULL)
  }
}

#' Compute correlation matrix
#' @param dataX data.frame with transition intensities per peptide
#' @export
#' @keywords internal
#' @family transitioncorrlation
#' @examples
#' list <- prolfqua_data('data_correlatedPeptideList')
#' cor_order(data_correlatedPeptideList[[1]])
cor_order <- function(dataX) {
  if (nrow(dataX) > 1) {
    ordt <- (dataX)[order(apply(dataX, 1, mean)), ]
    dd <- stats::cor(t(ordt), use = "pairwise.complete.obs", method = "spearman")
    return(dd)
  }else{
    message("Could not compute correlation, nr rows : " , nrow(dataX) )
  }
}

#' Compute correlation matrix with jackknife resampling
#
#' @param dataX data.frame with e.g. transition intensities per peptide or pepitde intensities per protein
#' @export
#' @keywords internal
#' @family transitioncorrlation
#' @examples
#' dd <- prolfqua_data("data_correlatedPeptideList")
#'
#' class(dd[[1]])
#' dd[[1]][1,2] <- NA
#' cor_jackknife_matrix(dd[[1]])
#'
cor_jackknife_matrix <- function(dataX,
                                       distmethod =
                                         function(x){cor(x, use = "pairwise.complete.obs", method = "pearson")}) {
  if (nrow(dataX) > 1) {
    ordt <- (dataX)[order(apply(dataX, 1, mean)), ]
    xpep <- t(ordt)
    res <- prolfqua::jackknife_matrix(xpep, distmethod)
  }else{
    message("Could not compute correlation, nr rows : ", nrow(dataX))
    res <- NULL
  }
  return(res)
}


#' normal pairs plot with different pch and plus abline
#' @param dataframe data matrix or data.frame as normally passed to pairs
#' @param ... params usually passed to pairs
#' @param legend  add legend to plots
#' @param pch point type default "."
#' @export
#' @family utilities
#' @keywords internal
#' @examples
#' tmp = matrix(rep((1:100),times = 4) + rnorm(100*4,0,3),ncol=4)
#' pairs_w_abline(tmp,log="xy",main="small data")
#' pairs_w_abline(tmp,log="xy",main="small data", legend=TRUE)
#' @seealso also \code{\link{pairs}}
pairs_w_abline <- function(dataframe,
                           legend = FALSE,
                           pch = ".",
                           ...) {
  pairs(
    dataframe,
    panel = function(x, y) {
      graphics::points(x, y, pch = pch)
      graphics::abline(
        a = 0,
        b = 1,
        v = 0,
        h = 0,
        col = 2
      )
      if (legend) {
        cR2 <- stats::cor(x, y, use = "pairwise.complete.obs") ^ 2
        graphics::legend("topleft",
                         legend = paste("R^2=", round(cR2, digits = 2) , sep = ""),
                         text.col = 3)
      }
    }
    ,
    lower.panel = NULL,
    ...
  )
}
#' histogram panel for pairs function (used as default in pairs_smooth)
#' @export
#' @family utilities
#' @keywords internal
#' @param x numeric data
#' @param ... additional parameters passed to rect
panel_hist <- function(x, ...)
{
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[seq_len(2)], 0, 1.5))
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts
  y <- y / max(y)
  rect(breaks[-nB], 0, breaks[-1], y,  ...)
}
#' correlation panel for pairs plot function (used as default in pairs_smooth)
#' @export
#' @family utilities
#' @keywords internal
#' @param x numeric data
#' @param y numeric data
#' @param ... not used
#' @param digits number of digits to display
panel_cor <- function(x, y, digits = 2, ...)
{
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y, use = "pairwise.complete.obs")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.7, txt)

  txt <- format(c(r ^ 2, 0.123456789), digits = digits)[1]
  txt <- bquote(R ^ {
    2
  } ~ "=" ~ .(txt))
  text(0.5, 0.5, txt)

  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if (p < 0.01)
    txt2 <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.3, txt2)
}

#' smoothScatter pairs
#' @param dataframe data matrix or data.frame as normally passed to pairs
#' @param legend  add legend to plots
#' @param ... params usually passed to pairs
#' @export
#' @family utilities
#' @keywords internal
#' @examples
#' tmp = matrix(rep((1:100),times = 4) + rnorm(100*4,0,3),ncol=4)
#' pairs_smooth(tmp,main="small data", legend=TRUE)
#' pairs_smooth(tmp,main="small data", diag.panel=panel_hist)
#' pairs_smooth(tmp,log="xy",main="small data", legend=TRUE)
#' @seealso also \code{\link{pairs}}
pairs_smooth <- function(dataframe, legend = FALSE, ...) {
  pairs(
    dataframe,
    upper.panel = function(x, y) {
      graphics::smoothScatter(x, y, add = TRUE)
      graphics::abline(
        a = 0,
        b = 1,
        v = 0,
        h = 0,
        col = 2
      )
      if (legend) {
        cR2 <- stats::cor(x , y, use = "pairwise.complete.obs") ^ 2
        graphics::legend("topleft",
                         legend = paste("R^2=", round(cR2, digits = 2) ,
                                        sep = ""),
                         text.col = 3)
      }
    }
    , lower.panel = panel_cor,
    ...
  )
}



#' table facade to easily switch implementations
#' @export
#' @family utilities
#' @keywords internal
table_facade <- function(df, caption, digits =  getOption("digits"), kable=TRUE){
  if (kable) {
    knitr::kable(df, digits = digits, caption = caption )
  }
}

#' table facade to easily switch implementations
#' @export
#' @family utilities
#' @keywords internal
table_facade.list <- function(parlist, kable=TRUE){
  table_facade(parlist$content, digits = parlist$digits, caption = parlist$caption )
}



.string.to.colors <- function(string, colors = NULL){
  if (is.factor(string)) {
    string = as.character(string)
  }
  if (!is.null(colors)) {
    if (length(colors) != length(unique(string))) {
      stop("The number of colors must be equal to the number of unique elements.")
    }
    else {
      conv = cbind(unique(string), colors)
    }
  } else {
    conv = cbind(unique(string), rainbow(length(unique(string))))
  }
  unlist(lapply(string, FUN = function(x){conv[which(conv[,1] == x),2]}))
}

.number.to.colors = function(value, colors = c("red", "blue"), num = 100){
  cols = colorRampPalette(colors)(num)
  cols = 	cols[findInterval(value, vec = seq(from = min(value), to = max(value), length.out = num))]
  cols
}

# @examples
#
# dataA <- data.frame(fc = c(-1,0,1,2,8), BFDR = c(0.01,1, 0.01, 0.005,0),
# condition = rep("A",5), Prey = LETTERS[1:5], modelName = c("A","A","B","A","A"))
#
# bc <- .volcano(dataA, contrast = "BLAB", effect = "fc", xintercept = 1, yintercept= 0.01, palette = c(A = "black" , B = "red"))
#
# bc |> plotly::subplot()
.volcano <- function(data,
                     contrast = NULL,
                     effect = "log2_EFCs",
                     significance = "BFDR",
                     x_annot = min(data[[effect]], na.rm = TRUE),
                     y_annot = min(data[[significance]], na.rm = TRUE),
                     proteinID = "Prey",
                     color = "modelName",
                     palette = NULL,
                     xintercept = c(-1,1),
                     yintercept = 0.1,
                     title_size = 25
){
  vline <- function(x = 0, color = "green") {
    list(
      type = "line",
      y0 = 0,
      y1 = 1,
      yref = "paper",
      x0 = x,
      x1 = x,
      line = list(color = color, dash = "dot")
    )
  }

  hline <- function(y = 0, color = "red") {
    list(
      type = "line",
      x0 = 0,
      x1 = 1,
      xref = "paper",
      y0 = y,
      y1 = y,
      line = list(color = color, dash = "dot")
    )
  }

  p <- plotly::plot_ly(data,
                       x = as.formula(paste0("~",effect)),
                       y = as.formula(paste0("~ I(-log10(", significance, "))")),
                       type = "scatter",
                       mode = "markers" ,
                       color = as.formula(paste0("~", color)),
                       colors = palette ,
                       text = as.formula(paste0("~", proteinID)) ,
                       showlegend = FALSE)

  sh2 <- c( lapply(xintercept, vline) , list(hline(-log10(yintercept))))
  p2 <- p |> plotly::layout(shapes = sh2) |>
    plotly::add_annotations(contrast,
                            x = x_annot,
                            y = -log10(y_annot),
                            showarrow = FALSE,
                            xanchor = "left",
                            font = list(size = title_size))
  return(p2)
}

#' volcano plotly
#' @param .data data frame
#' @param effect effect size x-axis
#' @param significance significance
#' @param proteinID column with protein ids
#' @param xintercept vertical abline at x
#' @param yintercept horizontal abline at y
#' @export
#' @examples
#'
#' data <- data.frame(fc = c(-1,0,1,2,8), BFDR = c(0.01,1, 0.01, 0.005,0),
#' condition = rep("A",5), Prey = LETTERS[1:5], modelName = c("A","A","B","A","A"))
#'
#' dataB <- data.frame(fc = c(-1,0,1,2,8), BFDR = c(0.01,1, 0.01, 0.005,0),
#' condition = rep("B",5), Prey = LETTERS[1:5],modelName = c("A","A","B","B","B"))
#' data <- dplyr::bind_rows(data, dataB)
#' bc <- volcano_plotly(data, xintercept = 1, yintercept= 0.01, palette = c(A = "black" , B = "red"))
#' bc |> plotly::subplot()
#'
volcano_plotly <- function(.data,
                           effect = "fc",
                           significance = "BFDR",
                           contrast = "condition",
                           proteinID = "Prey",
                           color = "modelName",
                           palette = NULL,
                           xintercept = c(-2,2),
                           yintercept = 0.1,
                           minsignificance = 1e-4,
                           title_size = 25

)
{
  .data[[significance]] <- pmax(.data[[significance]], minsignificance)

  xx <- .data |>
    dplyr::group_by(!!dplyr::sym(contrast)) |> tidyr::nest()

  makeshared <- function(x, proteinID = "Prey") {
    crosstalk::SharedData$new(x, as.formula( paste0("~", proteinID) ) ,
                              group = "Choose protein")
  }
  xx <- dplyr::mutate( xx,
                       shared_data = purrr::map( data,
                                                 makeshared, proteinID = proteinID  ))


  x_annot = min(.data[[effect]], na.rm = TRUE)
  y_annot = max(minsignificance, min(.data[[significance]], na.rm = TRUE))

  xd <- purrr::map2( xx$shared_data,
                     xx[[contrast]],
                     .volcano, effect = effect,
                     significance = significance,
                     x_annot = x_annot,
                     y_annot = y_annot,
                     proteinID = proteinID,
                     color = color,
                     palette = palette,
                     xintercept = xintercept,
                     yintercept = yintercept,
                     title_size = title_size
  )
  return(xd)
}



#' load data from prolfqua
#' @param datastr name of dataset
#' @param package default prolfqua
#' @export
#' @family data
prolfqua_data <- function(datastr, package="prolfqua"){
  data(list = datastr, package = package)
  return(eval(parse(text = datastr)))
}



