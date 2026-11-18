#' Internal error handler for tryCatch in model fitting and plotting: the error
#'
#' becomes a warning of class prolfqua_fit_failed, and the failed fit its message.
#' @param e condition raised by a failed fit
#' @noRd
.error_handler <- function(e) {
  rlang::warn(conditionMessage(e), class = "prolfqua_fit_failed")
  as.character(e)
}

#' Whether an optional (Suggests) package is installed; a function so tests can mock it.
#' @param pkg optional package name
#' @noRd
.has_package <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

#' Close graphics device `which` and make the previously active device current
#'
#' again, as withr::local_pdf() does when its frame exits.
#' @param which graphics device number
#' @noRd
.dev_close <- function(which) {
  prev <- grDevices::dev.prev(which)
  grDevices::dev.off(which)
  if (prev != which) {
    grDevices::dev.set(prev)
  }
  invisible(prev)
}

#' Internal helper for rendering ggplot/ComplexHeatmap objects to an open graphics device.
#' @param x ggplot or ComplexHeatmap object to draw
#' @noRd
.render_plot_to_device <- function(x) {
  if (methods::is(x, "Heatmap") || methods::is(x, "HeatmapList")) {
    ComplexHeatmap::draw(x)
  } else {
    plot(x)
  }
  invisible(x)
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
#' df <- data.frame(
#'   protein_Id = c("sp|P12345|PROT1_HUMAN", "sp|Q67890|PROT2_HUMAN",
#'                  "tr|A0A0A0|PROT3_HUMAN", "NOFASTA_ID"),
#'   value = 1:4
#' )
#' result <- get_uniprot_id_from_fasta_header(df, idcolumn = "protein_Id")
#' stopifnot("UniprotID" %in% colnames(result))
#' stopifnot(sum(!is.na(result$UniprotID)) == 3)
#'
get_uniprot_id_from_fasta_header <- function(df, idcolumn = "protein_Id") {
  map <- df |>
    dplyr::select(dplyr::all_of(idcolumn)) |>
    distinct() |>
    dplyr::filter(grepl(pattern = "^[a-z_]{2,7}\\|", !!sym(idcolumn))) |>
    tidyr::separate_wider_delim(
      dplyr::all_of(idcolumn),
      delim = "_",
      names = c("begin", "end"),
      too_few = "align_start",
      too_many = "drop",
      cols_remove = FALSE
    ) |>
    tidyr::separate_wider_delim(
      "begin",
      delim = "|",
      names = c("prefix", "UniprotID", "Symbol"),
      too_few = "align_start",
      too_many = "drop"
    ) |>
    dplyr::select(dplyr::all_of(c(idcolumn, "UniprotID")))
  if (nrow(map) > 0) {
    res <- dplyr::right_join(map, df, by = idcolumn)
  } else {
    warning("No SwissProt (sp) or Trembl (tr) id's were extracted.")
    res <- df
  }
  return(res)
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
#' @param max_nr_significant_text maximum number of significant labels to display
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
#' max_nr_significant_text = 300)
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
  max_nr_significant_text = 20
) {
  labeled_data <- tidyr::unite(.data, "label", dplyr::all_of(label))

  p <- .multigroup_volcano(
    labeled_data,
    effect = effect,
    significance = significance,
    contrast = contrast,
    colour = colour,
    xintercept = xintercept,
    yintercept = yintercept,
    text = "label",
    scales = scales
  )
  colname <- paste0("-log10(", significance, ")")

  if (!is.null(label)) {
    effect_values <- labeled_data[, effect]
    significance_values <- labeled_data[, significance]
    significant_subset <- subset(
      labeled_data,
      (effect_values < xintercept[1] | xintercept[2] < effect_values) & significance_values < yintercept
    ) |>
      head(n = max_nr_significant_text)
    if (nrow(significant_subset) > 0) {
      p <- p +
        ggrepel::geom_text_repel(
          data = significant_subset,
          aes(x = .data[[effect]], y = !!rlang::parse_expr(colname), label = .data$label),
          size = size,
          segment.size = segment.size,
          segment.alpha = segment.alpha
        )
    }
  }
  return(p)
}

#' Bound significance scores away from zero before log transformation
#' @param values numeric significance scores to bound away from zero
#' @param min_score smallest displayed significance score
#' @noRd
.floor_significance_values <- function(values, min_score = NULL) {
  if (!is.null(min_score)) {
    return(pmax(values, min_score))
  }

  positive_values <- values[is.finite(values) & values > 0]
  zero_floor <- if (length(positive_values) > 0) {
    min(positive_values, na.rm = TRUE)
  } else {
    .Machine$double.xmin
  }
  ifelse(values <= 0, zero_floor, values)
}

#' Choose colours for the contrast estimate types
#' @param colour estimate-type colour specification
#' @noRd
.contrast_colour_values <- function(colour) {
  if (identical(colour, "estimate_type")) {
    return(c(observed = "black", lod_imputed = "green", missing_fallback = "blue"))
  }
  c("black", "green", "blue", "red")
}

#' Draw one volcano panel per contrast with shared score definitions
#' @param data long-format measurement data frame
#' @param effect effect-size column name
#' @param significance significance-score column name
#' @param contrast name of the column identifying contrast panels
#' @param colour estimate-type colour specification
#' @param xintercept effect-size reference lines
#' @param yintercept significance reference lines
#' @param text name of the optional tooltip column
#' @param scales facet-axis scaling passed to ggplot2::facet_wrap
#' @noRd
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
  colname <- paste0("-log10(", significance, ")")
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
      values = .contrast_colour_values(colour)
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
#' @return The computed result.
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
    abort_bad_argument("x", "be a matrix or data.frame")
  }
  if (!is.null(preserve_row_names)) {
    row.names <- rownames(x)
    if (!is.null(row.names)) {
      dplyr::bind_cols(
        tibble::tibble(!!preserve_row_names := row.names),
        tibble::as_tibble(x, ...)
      )
    } else {
      warning("No row names to preserve. Object otherwise converted to tibble successfully.", call. = FALSE)
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
#' @examples
#' graphics::plot(1:3, 1:3, type = "n")
#' panel_cor(1:3, 1:3)
panel_cor <- function(x, y, digits = 2, ...) {
  usr <- par("usr")
  on.exit(par(usr = usr), add = TRUE)
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y, use = "pairwise.complete.obs")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0("r= ", txt)
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
  txt2 <- paste0("p= ", txt2)
  if (p < 0.01) {
    txt2 <- "p= <0.01"
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
        r_squared <- stats::cor(x, y, use = "pairwise.complete.obs")^2
        graphics::legend("topleft", legend = paste0("R^2=", round(r_squared, digits = 2)), text.col = 3)
      }
    },
    lower.panel = panel_cor,
    ...
  )
}


#' table facade to easily switch implementations
#' @return The requested plot, table, or transformed object.
#' @export
#' @family utilities
#' @param df data.frame to display
#' @param caption table caption
#' @param digits number of digits (default from options)
#' @param kable if TRUE use knitr::kable
#' @examples
#' table_facade(head(iris, 3), caption = "Iris preview")
table_facade <- function(df, caption, digits = getOption("digits"), kable = TRUE) {
  if (kable) {
    knitr::kable(df, digits = digits, caption = caption)
  }
}


#' Draw a volcano panel with effect and significance thresholds
#' @param data long-format measurement data frame
#' @param contrast label displayed as a plot annotation
#' @param effect effect-size column name
#' @param significance significance-score column name
#' @param x_annot horizontal position of the contrast annotation
#' @param y_annot significance value locating the annotation before the -log10 transform
#' @param proteinID feature identifier column
#' @param color colour column
#' @param palette named mapping from groups to colours
#' @param xintercept effect-size reference lines
#' @param yintercept significance reference lines
#' @param title_size plot title text size
#' @noRd
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
#' @param minsignificance optional minimum significance value (floor for
#'   -log10 axis). If NULL, only exact zero and negative values are replaced
#'   with the smallest positive observed value in the same score column.
#' @param title_size font size of the subplot title annotation
#' @param group crosstalk group name for linked brushing
#' @return The requested plot, table, or transformed object.
#' @export
#' @examplesIf requireNamespace("plotly", quietly = TRUE) && requireNamespace("crosstalk", quietly = TRUE)
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
  minsignificance = NULL,
  title_size = 25,
  group = "BB"
) {
  rlang::check_installed(c("plotly", "crosstalk"), reason = "for volcano_plotly().")
  .data[[significance]] <- .floor_significance_values(.data[[significance]], minsignificance)

  xx <- .data |>
    dplyr::group_by(!!dplyr::sym(contrast)) |>
    tidyr::nest()

  makeshared <- function(x, proteinID = "Prey") {
    crosstalk::SharedData$new(x, as.formula(paste0("~", proteinID)), group = group)
  }
  xx <- dplyr::mutate(xx, shared_data = purrr::map(data, makeshared, proteinID = proteinID))

  x_annot <- min(.data[[effect]], na.rm = TRUE)
  y_annot <- min(.data[[significance]], na.rm = TRUE)

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


#' load data from prolfqua
#' @param datastr name of dataset
#' @param package default prolfqua
#' @return The computed result.
#' @export
#' @family data
#' @examples
#' ionstar <- prolfqua_data("data_ionstar")
#' names(ionstar)
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
#' config$factor_depth <- 1
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
  intr <- purrr::map2_dfc(intr, newlev, function(x, lev) factor(x, levels = lev))

  colnames(intr) <- paste0("interaction_", columns)
  colname <- "interaction"
  data <- data |> dplyr::mutate(!!colname := interaction(intr, sep = sep))
  return(data)
}


# Decoy / contaminant identifier detection --------------------------------
# Shared, database-agnostic prefix detection used by both the protein
# annotation (deduplication) and the quant layer (flagging), so both sides
# agree on what a decoy / contaminant looks like.

# Anchored default prefixes. Kept conservative so a normal id is not flagged.
.default_decoy_prefixes <- "^REV_|^rev_|^DECOY|^decoy_|^XXX_|^reverse_|^##"
.default_contaminant_prefixes <- "^zz|^CON__|^CON_|^Cont_|^contam_"

#' defaults, unioned with an optional configured pattern. An empty / NULL /
#'
#' no-op ("a^") pattern yields the defaults only -- never grepl("", x), which
#' would match every id (the exact failure mode behind WU347806).
#' @param pattern optional user-specified identifier pattern
#' @param defaults built-in identifier patterns
#' @noRd
.effective_prefix_pattern <- function(pattern, defaults) {
  has_pat <- !is.null(pattern) &&
    length(pattern) == 1 &&
    !is.na(pattern) &&
    nzchar(pattern) &&
    !identical(pattern, "a^")
  if (has_pat) paste(pattern, defaults, sep = "|") else defaults
}

#' Detect decoy / reverse-database identifiers
#'
#' Flags identifiers that look like decoy (reversed) database entries. Built-in
#' anchored default prefixes are always considered, unioned with an optional
#' configured \code{pattern}. An empty string or \code{NULL} pattern uses the
#' defaults only.
#' @export
#' @family utilities
#' @param ids character vector of (prefixed) identifiers
#' @param pattern optional additional decoy regex, unioned with the defaults
#' @return logical vector, \code{TRUE} where an id looks like a decoy
#' @examples
#' is_decoy(c("REV_sp|P1|X", "sp|P2|X", "decoy_3", "normalProtein"))
#' is_decoy(c("shuffled_1", "P2"), pattern = "^shuffled_")
is_decoy <- function(ids, pattern = NULL) {
  grepl(
    .effective_prefix_pattern(pattern, .default_decoy_prefixes),
    as.character(ids)
  )
}

#' Detect contaminant identifiers
#'
#' As \code{\link{is_decoy}}, but for contaminant entries (keratin, trypsin,
#' BSA, ...). Built-in anchored default prefixes unioned with an optional
#' configured \code{pattern}; empty / \code{NULL} uses the defaults only.
#' @export
#' @family utilities
#' @param ids character vector of (prefixed) identifiers
#' @param pattern optional additional contaminant regex, unioned with defaults
#' @return logical vector, \code{TRUE} where an id looks like a contaminant
#' @examples
#' is_contaminant(c("zz|Cont00001|X", "sp|P2|X", "CON__ALBU"))
is_contaminant <- function(ids, pattern = NULL) {
  grepl(
    .effective_prefix_pattern(pattern, .default_contaminant_prefixes),
    as.character(ids)
  )
}
