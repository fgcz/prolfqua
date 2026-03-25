# Internal error handler for tryCatch in model fitting and plotting
.error_handler <- function(e) {
  warning("WARN :", e)
  as.character(e)
}

#' Extracts uniprot ID
#'
#' @export
#' @family utilities
#' @return data.frame
#' @param df data.frame
#' @param idcolumn character, column name to extract uniprot IDs from
#'
#' @examples
#'
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' tmp <- prolfqua::separate_hierarchy(bb$data,old2new( bb$config))
#' tmp$UniprotID <- NULL
#' tmp <- get_UniprotID_from_fasta_header(tmp, idcolumn = "top_protein")
#' stopifnot("UniprotID" %in%  colnames(tmp))
#'
get_UniprotID_from_fasta_header <- function(df, idcolumn = "protein_Id") {
  map <- df |>
    dplyr::select(dplyr::all_of(idcolumn)) |>
    distinct() |>
    dplyr::filter(grepl(pattern = "^[a-z_]{2,7}\\|", !!sym(idcolumn))) |>
    tidyr::separate(col = !!sym(idcolumn), sep = "_", into = c("begin", "end"), remove = FALSE) |>
    tidyr::separate(col = !!sym("begin"), sep = "\\|", into = c("prefix", "UniprotID", "Symbol")) |>
    dplyr::select(-!!sym("prefix"), -!!sym("Symbol"), -!!sym("end"))
  if (nrow(map) > 0) {
    res <- dplyr::right_join(map, df, by = idcolumn)
  } else {
    warning("No SwissProt (sp) or Trembl (tr) id's were extracted.")
    res <- df
  }
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
#'
remove_NA_rows <- function(obj, thresh = 0) {
  x <- apply(obj, 1, function(x) {
    sum(is.na(x))
  })
  obj <- obj[!(x > thresh), ]
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
names_to_matrix <- function(names, split = "\\||\\_") {
  cnamessplit <- stringi::stri_split_regex(as.character(names), pattern = split)
  protnam <- do.call("rbind", cnamessplit)
  return(protnam)
}


#' plot volcano given multiple contrasts
#' @param .data data in long format
#' @param effect column containing effect sizes
#' @param significance column containing p-values, q.values etc
#' @param contrast column with contrast
#' @param colour colouring of points
#' @param xintercept fc thresholds
#' @param yintercept significance threshold
#' @param label column containing labels
#' @param size controls size of text
#' @param segment.size controls size of lines
#' @param segment.alpha controls visibility of lines
#' @param scales parameter to ggplot2::facet_wrap
#' @param maxNrOfSignificantText maximum number of significant labels to display
#' @export
#' @keywords internal
#' @family utilities
#' @examples
#'
#'
#' show <- data.frame(logFC = rnorm(100, 0, 1), adj.P.Val = runif(100, 0, 1),
#'   Condition = sample(c("a","b")), colour = "forward", Name = paste0("n", 1:100))
#' prolfqua::multigroup_volcano( show,
#' effect="logFC",
#' significance = "adj.P.Val",
#' contrast="Condition",
#' colour=NULL,
#' label="Name",
#' maxNrOfSignificantText = 300)
multigroup_volcano <- function(
  .data,
  effect = "fc",
  significance = "p.adjust",
  contrast = "condition",
  colour = "colour",
  xintercept = c(-2, 2),
  yintercept = 0.05,
  label = NULL,
  size = 1,
  segment.size = 0.3,
  segment.alpha = 0.3,
  scales = "fixed",
  maxNrOfSignificantText = 20
) {
  misspX <- tidyr::unite(.data, "label", dplyr::all_of(label))

  p <- .multigroup_volcano(
    misspX,
    effect = effect,
    significance = significance,
    contrast = contrast,
    colour = colour,
    xintercept = xintercept,
    yintercept = yintercept,
    text = "label",
    scales = scales
  )
  colname <- paste("-log10(", significance, ")", sep = "")

  if (!is.null(label)) {
    effectX <- misspX[, effect]
    typeX <- misspX[, significance]
    subsetData <- subset(misspX, (effectX < xintercept[1] | xintercept[2] < effectX) & typeX < yintercept) |>
      head(n = maxNrOfSignificantText)
    if (nrow(subsetData) > 0) {
      p <- p +
        ggrepel::geom_text_repel(
          data = subsetData,
          aes(x = .data[[effect]], y = !!rlang::parse_expr(colname), label = .data$label),
          size = size,
          segment.size = segment.size,
          segment.alpha = segment.alpha
        )
    }
  }
  return(p)
}

.multigroup_volcano <- function(
  data,
  effect = "fc",
  significance = "p.adjust",
  contrast = "contrast",
  colour = "colour",
  xintercept = c(-2, 2),
  yintercept = 0.05,
  text = NULL,

  scales = "fixed"
) {
  colname <- paste("-log10(", significance, ")", sep = "")
  p <- ggplot(
    data,
    aes(x = .data[[effect]], y = !!rlang::parse_expr(colname))
  ) +
    geom_point(alpha = 0.5)
  if (!is.null(colour)) {
    p <- p + aes(color = .data[[colour]])
  }
  if (!is.null(text)) {
    p <- p + aes(text = .data[[text]])
  }
  p <- p +
    scale_colour_manual(
      values = c("black", "green", "blue", "red")
    )
  p <- p +
    facet_wrap(
      as.formula(paste("~", contrast)),
      scales = scales
    ) +
    labs(y = colname)
  if (!is.null(xintercept)) {
    p <- p +
      geom_vline(
        xintercept = xintercept,
        linetype = "dashed",
        colour = "red"
      )
  }
  if (!is.null(yintercept)) {
    p <- p +
      geom_hline(
        yintercept = -log10(yintercept),
        linetype = "dashed",
        color = "blue",
        show.legend = FALSE
      )
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
#' @examples
#' x <- matrix(rnorm(20), ncol=4)
#' rownames(x) <- LETTERS[seq_len(nrow(x))]
#' matrix_to_tibble(x)
#' !(is.matrix(x) || is.data.frame(x))
#'
matrix_to_tibble <- function(x, preserve_row_names = "row.names", ...) {
  if (!(is.matrix(x) || is.data.frame(x))) {
    stop("Error: `x` is not a matrix or data.frame object.")
  }
  if (!is.null(preserve_row_names)) {
    row.names <- rownames(x)
    if (!is.null(row.names)) {
      dplyr::bind_cols(
        tibble::tibble(!!preserve_row_names := row.names),
        tibble::as_tibble(x, ...)
      )
    } else {
      warning("Warning: No row names to preserve. ", "Object otherwise converted to tibble successfully.")
      tibble::as_tibble(x, ...)
    }
  } else {
    tibble::as_tibble(x, ...)
  }
}

#' correlation panel for pairs plot function (used as default in pairs_smooth)
#' @export
#' @family utilities
#' @keywords internal
#' @param x numeric data
#' @param y numeric data
#' @param ... not used
#' @param digits number of digits to display
panel_cor <- function(x, y, digits = 2, ...) {
  usr <- par("usr")
  on.exit(par(usr = usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y, use = "pairwise.complete.obs")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.7, txt)

  txt <- format(c(r^2, 0.123456789), digits = digits)[1]
  txt <- bquote(
    R^{
      2
    } ~ "=" ~ .(txt)
  )
  text(0.5, 0.5, txt)

  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if (p < 0.01) {
    txt2 <- paste("p= ", "<0.01", sep = "")
  }
  text(0.5, 0.3, txt2)
}

#' smoothScatter pairs
#' @param data data matrix or data.frame as normally passed to pairs
#' @param legend  add legend to plots
#' @param ... params usually passed to pairs
#' @export
#' @family utilities
#' @keywords internal
#' @examples
#' tmp = matrix(rep((1:100),times = 4) + rnorm(100*4,0,3),ncol=4)
#' pairs_smooth(tmp,main="small data", legend=TRUE)
#' pairs_smooth(tmp,main="small data")
#' pairs_smooth(tmp,log="xy",main="small data", legend=TRUE)
#' @seealso also \code{\link{pairs}}
pairs_smooth <- function(data, legend = FALSE, ...) {
  non_na_counts <- apply(data, 2, function(x) sum(!is.na(x)))
  min_col <- which(non_na_counts < 2)
  if (length(min_col) >= 1) {
    data <- data[, -min_col, drop = FALSE]
  }
  pairs(
    data,
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
        cR2 <- stats::cor(x, y, use = "pairwise.complete.obs")^2
        graphics::legend("topleft", legend = paste("R^2=", round(cR2, digits = 2), sep = ""), text.col = 3)
      }
    },
    lower.panel = panel_cor,
    ...
  )
}


#' table facade to easily switch implementations
#' @export
#' @family utilities
#' @param df data.frame to display
#' @param caption table caption
#' @param digits number of digits (default from options)
#' @param kable if TRUE use knitr::kable
table_facade <- function(df, caption, digits = getOption("digits"), kable = TRUE) {
  if (kable) {
    knitr::kable(df, digits = digits, caption = caption)
  }
}


.volcano <- function(
  data,
  contrast = NULL,
  effect = "log2_EFCs",
  significance = "BFDR",
  x_annot = min(data[[effect]], na.rm = TRUE),
  y_annot = min(data[[significance]], na.rm = TRUE),
  proteinID = "Prey",
  color = "modelName",
  palette = NULL,
  xintercept = c(-1, 1),
  yintercept = 0.1,
  title_size = 25
) {
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

  p <- plotly::plot_ly(
    data,
    x = as.formula(paste0("~", effect)),
    y = as.formula(paste0("~ I(-log10(", significance, "))")),
    type = "scatter",
    mode = "markers",
    color = as.formula(paste0("~", color)),
    colors = palette,
    text = as.formula(paste0("~", proteinID)),
    showlegend = FALSE
  )

  sh2 <- c(lapply(xintercept, vline), list(hline(-log10(yintercept))))
  p2 <- p |>
    plotly::layout(shapes = sh2) |>
    plotly::add_annotations(
      contrast,
      x = x_annot,
      y = -log10(y_annot),
      showarrow = FALSE,
      xanchor = "left",
      font = list(size = title_size)
    )
  return(p2)
}

#' volcano plotly
#' @param .data data frame
#' @param effect effect size x-axis
#' @param significance significance
#' @param contrast column with contrast labels
#' @param proteinID column with protein ids
#' @param color column used for colouring points
#' @param palette named colour vector for the colour aesthetic
#' @param xintercept vertical abline at x
#' @param yintercept horizontal abline at y
#' @param minsignificance minimum significance value (floor for -log10 axis)
#' @param title_size font size of the subplot title annotation
#' @param group crosstalk group name for linked brushing
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
volcano_plotly <- function(
  .data,
  effect = "fc",
  significance = "BFDR",
  contrast = "condition",
  proteinID = "Prey",
  color = "modelName",
  palette = NULL,
  xintercept = c(-2, 2),
  yintercept = 0.1,
  minsignificance = 1e-4,
  title_size = 25,
  group = "BB"
) {
  .data[[significance]] <- pmax(.data[[significance]], minsignificance)

  xx <- .data |>
    dplyr::group_by(!!dplyr::sym(contrast)) |>
    tidyr::nest()

  makeshared <- function(x, proteinID = "Prey") {
    crosstalk::SharedData$new(x, as.formula(paste0("~", proteinID)), group = group)
  }
  xx <- dplyr::mutate(xx, shared_data = purrr::map(data, makeshared, proteinID = proteinID))

  x_annot <- min(.data[[effect]], na.rm = TRUE)
  y_annot <- max(minsignificance, min(.data[[significance]], na.rm = TRUE))

  xd <- purrr::map2(
    xx$shared_data,
    xx[[contrast]],
    .volcano,
    effect = effect,
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


.scatter <- function(
  data,
  contrast = NULL,
  dx = "diff.protein",
  dy = "diff.site",
  x_annot = min(data[[dx]], na.rm = TRUE),
  y_annot = min(data[[dy]], na.rm = TRUE),
  proteinID = "Prey",
  color = "modelName.site",
  palette = NULL,
  title_size = 25
) {
  p <- plotly::plot_ly(
    data,
    x = as.formula(paste0("~", dx)),
    y = as.formula(paste0("~", dy)),
    type = "scatter",
    mode = "markers",
    color = as.formula(paste0("~", color)),
    colors = palette,
    text = as.formula(paste0("~", proteinID)),
    showlegend = FALSE
  )

  p <- p |>
    plotly::add_annotations(
      contrast,
      x = x_annot,
      y = y_annot,
      showarrow = FALSE,
      xanchor = "left",
      font = list(size = title_size)
    )
  return(p)
}


#' scatter plotly
#' @param .data data frame
#' @param dx column name for the x-axis difference
#' @param dy column name for the y-axis difference
#' @param contrast column with contrast labels
#' @param proteinID column with protein ids
#' @param color column used for colouring points
#' @param palette named colour vector for the colour aesthetic
#' @param title_size font size of the subplot title annotation
#' @param group crosstalk group name for linked brushing
#' @export
#' @examples
#'
#' data <- data.frame(diff.protein = c(-1,0,1,2,8), diff.site = c(0.01,1, 0.01, 0.005,0),
#' condition = rep("A",5), protein_Id = LETTERS[1:5], modelName = c("A","A","B","A","A"))
#'
#' dataB <- data.frame(diff.protein = c(-1,0,1,2,8), diff.site = c(0.01,1, 0.01, 0.005,0),
#' condition = rep("B",5), protein_Id = LETTERS[1:5],modelName = c("A","A","B","B","B"))
#' data <- dplyr::bind_rows(data, dataB)
#' bc <- scatter_plotly(data, palette = c(A = "black" , B = "red"))
#' bc[[1]]
#' bc[[2]]
#' bc |> plotly::subplot()
#'
scatter_plotly <- function(
  .data,
  dx = "diff.protein",
  dy = "diff.site",
  contrast = "condition",
  proteinID = "protein_Id",
  color = "modelName",
  palette = NULL,
  title_size = 25,
  group = "BB"
) {
  xx <- .data |>
    dplyr::group_by(!!dplyr::sym(contrast)) |>
    tidyr::nest()

  makeshared <- function(x, proteinID = "Prey") {
    crosstalk::SharedData$new(x, as.formula(paste0("~", proteinID)), group = group)
  }
  xx <- dplyr::mutate(
    xx,
    shared_data = purrr::map(
      data,
      makeshared,
      proteinID = proteinID
    )
  )

  x_annot <- min(.data[[dx]], na.rm = TRUE)
  y_annot <- min(.data[[dy]], na.rm = TRUE)

  xd <- purrr::map2(
    xx$shared_data,
    xx[[contrast]],
    .scatter,
    dx = dx,
    dy = dy,
    x_annot = x_annot,
    y_annot = y_annot,
    proteinID = proteinID,
    color = color,
    palette = palette,
    title_size = title_size
  )
  return(xd)
}


#' load data from prolfqua
#' @param datastr name of dataset
#' @param package default prolfqua
#' @export
#' @family data
prolfqua_data <- function(datastr, package = "prolfqua") {
  e <- new.env(parent = emptyenv())
  data(list = datastr, package = package, envir = e)
  return(e[[datastr]])
}


#' create interaction column from factors
#' @export
#' @keywords internal
#' @family configuration
#' @examples
#' xx <- data.frame(A = c("a","a","a"), B = c("d","d","e"))
#' x <- make_interaction_column(xx, c("B","A"))
#' x <- make_interaction_column(xx, c("A"))
#' bb <- prolfqua::sim_lfq_data_protein_config()
#' config <- bb$config
#' analysis <- bb$data
#'
#' config$factorDepth <- 1
#' make_interaction_column(analysis,
#'    config$factor_keys_depth())
#'
make_interaction_column <- function(data, columns, sep = ".") {
  intr <- dplyr::select(data, dplyr::all_of(columns))
  intr <- purrr::map_dfc(intr, factor)
  names(columns) <- columns
  newlev <- purrr::map2(columns, intr, function(x, y) {
    paste0(x, levels(y))
  })
  intr <- purrr::map2_dfc(columns, intr, paste0)
  intr <- purrr::map2_dfc(intr, newlev, forcats::fct_relevel)

  colnames(intr) <- paste0("interaction_", columns)
  colname <- "interaction"
  data <- data |> dplyr::mutate(!!colname := interaction(intr, sep = sep))
  return(data)
}
