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
#' x1 = removeNArows(obj, thresh=0)
#' stopifnot(all(c(9,10)==dim(x1)))
#' x2 = removeNArows(obj, thresh=1)
#' stopifnot(all(c(10,10)==dim(x2)))
removeNArows <- function(obj, thresh=0 )
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
#' split2table(dat,split="\\_|\\/")
split2table <- function(names,split="\\||\\_")
{
  cnamessplit <- stringr::str_split(as.character(names),pattern = split)
  protnam <- do.call("rbind",cnamessplit)
  return(protnam)
}

#' plot volcano given multiple conditions
#' @param misspX data in long format
#' @param effect column containing effect sizes
#' @param p.value column containing p-values, q.values etc
#' @param condition column with condition
#' @param colour colouring of points
#' @param xintercept fc thresholds
#' @param pvalue pvalue threshold
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
#' library(ggplot2)
#' library(tidyverse)
#' library(ggrepel)
#' multigroupVolcano(LFQServiceData::multigroupFCDATA,effect="logFC",p.value="adj.P.Val",condition="Condition",colour="colour",label="Name" )
multigroupVolcano <- function(misspX,
                              effect = "fc",
                              p.value = "p.adjust",
                              condition = "condition",
                              colour = "colour",
                              xintercept = c(-2,2),
                              pvalue = 0.05,
                              label = NULL,
                              size = 1,
                              segment.size = 0.3,
                              segment.alpha = 0.3,
                              ablines = data.frame(
                                fc = c(0,0),
                                p = c(0.01,0.05),
                                Area = c('p=0.01','p=0.05')
                              ), scales = "fixed",
                              maxNrOfSignificantText = 20)
{
  misspX <- tidyr::unite(misspX, "label", label)
  colname = paste("-log10(", p.value , ")" , sep = "")
  p <- ggplot( misspX, aes_string(x = effect , y = colname, color = colour  )  )  +
    geom_point(alpha = 0.5)
  p <- p + scale_colour_manual(values = c("black", "green", "blue","red"))
  p <- p + facet_wrap(as.formula(paste("~",condition)),scales = scales) + labs(y = colname)

  ablines$neg_log10p <- -log10(ablines$p)
  p <- p + geom_abline(data = ablines, aes_string(slope = "fc", intercept = "neg_log10p",colour = "Area")) +
    geom_vline(xintercept = xintercept,linetype = "dashed", colour = "red")

  if (!is.null(label)) {
    effectX <- misspX[,effect]
    typeX <- misspX[,p.value]
    subsetData <- subset(misspX, (effectX < xintercept[1] | xintercept[2] < effectX) & typeX < pvalue ) %>% head(maxNrOfSignificantText)
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
#' rownames(x) <- LETTERS[1:nrow(x)]
#' matrix_to_tibble(x)
#' !(is.matrix(x) || is.data.frame(x))
#'
matrix_to_tibble <- function(x, preserve_row_names = "row.names", ... )
{
  if (!(is.matrix(x) || is.data.frame(x))) stop("Error: `x` is not a matrix or data.frame object.")
  if (!is.null(preserve_row_names)) {
    row.names <- rownames(x)
    if (!is.null(row.names)  ) {
      #&& !identical(row.names, 1:nrow(x) %>% as.character())
      dplyr::bind_cols(
        tibble::tibble(!! preserve_row_names := row.names),
        tibble::as_tibble(x, ...)
      )

    } else {

      warning(paste0("Warning: No row names to preserve. ",
                     "Object otherwise converted to tibble successfully."))
      tibble::as_tibble(x, ...)
    }

  } else {

    tibble::as_tibble(x, ...)

  }
}


# @TODO think of making it public.
#' copute jack knive
#' @param xdata matrix
#' @param .method method i.e. cor, parameters
#' @param ... further parameters to .method
#' @return list with all jackknife matrices
#' @export
#' @family utilities
#' @keywords internal
#' @examples
#' xx <- matrix(rnorm(20), ncol=4)
#' cortest <- function(x){print(dim(x));cor(x)}
#' my_jackknife(xx, cortest)
#' my_jackknife(xx, cor, use="pairwise.complete.obs", method="pearson")
my_jackknife <- function(xdata, .method, ... ) {
  x <- 1:nrow(xdata)
  call <- match.call()
  n <- length(x)
  u <- vector( "list", length = n )
  for (i in 1:n) {
    tmp <- xdata[x[-i],]
    u[[i]] <- .method(tmp, ...)
  }
  names(u) <- 1:n
  thetahat <- .method(xdata, ...)
  invisible(list(thetahat = thetahat, jack.values = u, call = call ))
}

#' Compute correlation matrix with jack
#' @param dataX data.frame with transition intensities per peptide
#' @param distmethod dist or correlation method working with matrix i.e. cor
#' @param ... further parameters to method
#' @export
#' @family utilities
#' @keywords internal
#' @return summarizes results producced with my_jackknife
#' @examples
#' dataX <- matrix(rnorm(20), ncol=4)
#' rownames(dataX)<- paste("R",1:nrow(dataX),sep="")
#' colnames(dataX)<- paste("C",1:ncol(dataX),sep="")
#' tmp <- my_jackknife(dataX, cor, use="pairwise.complete.obs", method="pearson")
#'
#' jackknifeMatrix(dataX, cor)
#' jackknifeMatrix(dataX, cor, method="spearman")
jackknifeMatrix <- function(dataX, distmethod , ... ){
  if (is.null(colnames(dataX))) {
    colnames(dataX) <- paste("C", 1:ncol(dataX), sep = "")
  }
  if (is.null(rownames(dataX))) {
    rownames(dataX) <- paste("R", 1:nrow(dataX), sep = "")
  }

  if (nrow(dataX) > 1 & ncol(dataX) > 1) {
    tmp <- my_jackknife( dataX, distmethod, ... )
    x <- plyr::ldply(tmp$jack.values, LFQService::matrix_to_tibble)
    dd <- tidyr::gather(x, "col.names" , "correlation" , 3:ncol(x))
    ddd <- dd %>%
      group_by(UQ(sym("row.names")), UQ(sym("col.names"))) %>%
      summarize_at(c("jcor" = "correlation"), function(x){max(x, na.rm = TRUE)})

    dddd <- tidyr::spread(ddd, UQ(sym("col.names")), UQ(sym("jcor"))  )
    dddd <- as.data.frame(dddd)
    rownames(dddd) <- dddd$row.names
    dddd <- dddd[,-1]
    return(dddd)
  }else{
    message("Could not compute correlation, nr rows : " , nrow(dataX) )
  }
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
panel.hist <- function(x, ...)
{
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
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
panel.cor <- function(x, y, digits = 2, ...)
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
#' pairs_smooth(tmp,main="small data", diag.panel=panel.hist)
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
    , lower.panel = panel.cor,
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
