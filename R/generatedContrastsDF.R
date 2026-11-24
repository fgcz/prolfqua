#' Group labels of a two-factor design
#'
#' Names the cell of a two-factor design as `G_<primary>_<secondary>`, the group names that the contrast generators
#' ([generate_contrasts()] and its parts) use in their contrast strings.
#'
#' @param primary character, primary factor level
#' @param secondary character, secondary factor level
#' @return A character vector with combined group labels.
#' @export
#' @family modelling
#' @examples
#' group_label("a", "b")
group_label <- function(primary, secondary) paste0("G_", primary, "_", secondary)


#' Main-effect contrasts of a two-factor design
#'
#' For every pair of primary levels, the difference of their means averaged over all secondary levels.
#'
#' @param primary_levels character vector of primary factor levels
#' @param secondary_levels character vector of secondary factor levels
#' @return A named list of contrast strings, one per pair of primary levels, named `<high>_vs_<low>`; groups are
#'   named by [group_label()].
#' @export
#' @family contrasts
#' @examples
#' primary_levels <- c("KO", "WT")
#' secondary_levels <- c("T150", "T0", "T300")
#' main_effect_contrasts(primary_levels, secondary_levels)
#' main_effect_contrasts(secondary_levels, primary_levels)
main_effect_contrasts <- function(primary_levels, secondary_levels) {
  contrasts <- list()

  for (pair in combn(primary_levels, 2, simplify = FALSE)) {
    low <- pair[1]
    high <- pair[2]
    avg_high <- paste0(group_label(high, secondary_levels), collapse = " + ")
    avg_low <- paste0(group_label(low, secondary_levels), collapse = " + ")
    key <- paste0(high, "_vs_", low)
    contrasts[[key]] <- sprintf(
      "(%s)/%d - (%s)/%d",
      avg_high,
      length(secondary_levels),
      avg_low,
      length(secondary_levels)
    )
  }
  contrasts
}

#' Level-specific contrasts of a two-factor design
#'
#' For every pair of primary levels, their difference within each secondary level.
#'
#' @param primary_levels character vector of primary factor levels
#' @param secondary_levels character vector of secondary factor levels
#' @return A named list of contrast strings, one per pair of primary levels and secondary level, named
#'   `<high>_vs_<low>_at_<secondary>`.
#' @export
#' @family contrasts
#' @examples
#' primary_levels <- c("KO", "WT")
#' secondary_levels <- c("T0", "T150", "T300")
#' x2 <- level_specific_contrasts(primary_levels, secondary_levels)
#' x3 <- level_specific_contrasts(secondary_levels, primary_levels)
level_specific_contrasts <- function(primary_levels, secondary_levels) {
  contrasts <- list()
  for (pair in combn(primary_levels, 2, simplify = FALSE)) {
    low <- pair[1]
    high <- pair[2]
    for (sec in secondary_levels) {
      key <- paste0(high, "_vs_", low, "_at_", sec)
      contrasts[[key]] <- sprintf(
        "%s - %s",
        group_label(high, sec),
        group_label(low, sec)
      )
    }
  }
  contrasts
}


#' Interaction contrasts of a two-factor design
#'
#' For every pair of primary levels and every pair of secondary levels, the difference of differences: does the
#' primary difference change between the two secondary levels?
#'
#' @param primary_levels character vector of primary factor levels
#' @param secondary_levels character vector of secondary factor levels
#' @return A named list of contrast strings, named `interaction_<high>_vs_<low>_at_<s2>_vs_<s1>`.
#' @export
#' @family contrasts
#' @examples
#' primary_levels <- c("KO", "WT")
#' secondary_levels <- c("T0", "T150", "T300")
#' interaction_contrasts(primary_levels, secondary_levels)
#' interaction_contrasts(secondary_levels, primary_levels)
interaction_contrasts <- function(primary_levels, secondary_levels) {
  contrasts <- list()
  for (pair in combn(primary_levels, 2, simplify = FALSE)) {
    low <- pair[1]
    high <- pair[2]
    for (sp in combn(secondary_levels, 2, simplify = FALSE)) {
      s1 <- sp[1]
      s2 <- sp[2]
      key <- paste0("interaction_", high, "_vs_", low, "_at_", s2, "_vs_", s1)
      contrasts[[key]] <- sprintf(
        "(%s - %s) - (%s - %s)",
        group_label(high, s2),
        group_label(low, s2),
        group_label(high, s1),
        group_label(low, s1)
      )
    }
  }
  contrasts
}

#' Pairwise contrasts of one factor
#'
#' One contrast for every pair of levels: the later level minus the earlier one.
#'
#' @param levels character vector of factor levels
#' @return A named list of contrast strings `"<high> - <low>"`, named `<high>_vs_<low>`.
#' @export
#' @family contrasts
#' @examples
#' primary_levels <- c("KO", "WT")
#' secondary_levels <- c("T0", "T150", "T300")
#' x2 <- level_specific_contrasts(primary_levels, secondary_levels)
#' x3 <- level_specific_contrasts(secondary_levels, primary_levels)
#' generate_contrasts_for_factor(names(x2))
generate_contrasts_for_factor <- function(levels) {
  contrasts <- list()
  for (pair in combn(levels, 2, simplify = FALSE)) {
    low <- pair[1]
    high <- pair[2]
    key <- paste0(high, "_vs_", low)
    contrasts[[key]] <- paste0(high, " - ", low)
  }
  contrasts
}

#' All contrasts of a two-factor design
#'
#' Combines [main_effect_contrasts()], [level_specific_contrasts()] and, optionally, [interaction_contrasts()].
#'
#' @param primary_levels character vector of primary factor levels
#' @param secondary_levels character vector of secondary factor levels
#' @param interactions logical, if TRUE include interaction contrasts
#' @return A data.frame with the columns `ContrastName` and `Contrast` (the contrast string), one row per contrast.
#' @export
#' @family contrasts
#' @examples
#' primary_levels <- c("KO", "WT")
#' secondary_levels <- c("T0", "T150", "T300")
#' generate_contrasts(primary_levels, secondary_levels)
#' generate_contrasts(secondary_levels, primary_levels)
generate_contrasts <- function(primary_levels, secondary_levels, interactions = TRUE) {
  res <- c(
    unlist(main_effect_contrasts(primary_levels, secondary_levels)),
    unlist(level_specific_contrasts(primary_levels, secondary_levels)),
    if (interactions) {
      unlist(interaction_contrasts(primary_levels, secondary_levels))
    } else {
      NULL
    }
  )
  res <- data.frame(ContrastName = names(res), Contrast = res)
  res
}

#' Add the contrasts of a two-factor design to a sample annotation
#'
#' Unites the two factor columns into a `Group` column and appends the contrasts of [generate_contrasts()] as the
#' columns `ContrastName` and `Contrast`, in the layout of a prolfquapp dataset annotation. The shorter side is padded
#' with `NA` rows: with fewer contrasts than samples the contrast columns, with more contrasts than samples the sample
#' columns, so the extra rows carry only a contrast. Nothing is written to disk.
#'
#' @param df data.frame with annotation data
#' @param primary_col character, name of the primary factor column
#' @param secondary_col character, name of the secondary factor column
#' @param prefix character, prefix for output naming (default "G_")
#' @param dataset_id character, dataset identifier for output naming
#' @param decreasing logical, if TRUE sort factor levels in decreasing order
#' @param interactions logical, if TRUE include interaction contrasts
#' @return A list with `annot`, the annotation with the `Group`, `ContrastName` and `Contrast` columns (the two factor
#'   columns are dropped), and `name`, a suggested file name `DEA_<prefix>_<dataset_id>.csv`.
#' @export
#' @family contrasts
#' @examples
#'
#' # 2 genotypes x 3 time points x 6 replicates
#' design <- expand.grid(
#'   rep = 1:6, factor_B = c("T0", "T150", "T300"), factor_A = c("WT", "KO"),
#'   stringsAsFactors = FALSE
#' )
#' annotation <- data.frame(
#'   Name = paste(design$factor_A, design$factor_B, design$rep, sep = "_"),
#'   factor_A = design$factor_A,
#'   factor_B = design$factor_B
#' )
#' annotation_add_contrasts(annotation, "factor_A", "factor_B", "primary")
#' annotation_add_contrasts(annotation, "factor_B", "factor_A", "secondary")
annotation_add_contrasts <- function(
  df,
  primary_col,
  secondary_col,
  prefix = "G_",
  dataset_id = "dataset",
  decreasing = FALSE,
  interactions = TRUE
) {
  levels_p <- sort(unique(df[[primary_col]]), decreasing = decreasing)
  levels_s <- sort(unique(df[[secondary_col]]), decreasing = decreasing)

  df <- df |> tidyr::unite("Group", dplyr::all_of(c(primary_col, secondary_col)), sep = "_", remove = FALSE)
  cells <- unique(df[c("Group", primary_col, secondary_col)])
  shared <- cells[cells$Group %in% cells$Group[duplicated(cells$Group)], ]
  if (nrow(shared) > 0) {
    abort_bad_argument(
      "df",
      paste0("give each combination of `", primary_col, "` and `", secondary_col, "` its own group name"),
      not = toString(paste0(shared[[primary_col]], " + ", shared[[secondary_col]], " as ", shared$Group))
    )
  }
  df_contrast <- generate_contrasts(levels_p, levels_s, interactions = interactions)

  n_rows <- max(nrow(df), nrow(df_contrast))
  annot <- dplyr::bind_cols(
    pad_rows(df |> dplyr::select(-dplyr::all_of(c(primary_col, secondary_col))), n_rows),
    pad_rows(df_contrast, n_rows)
  )
  name <- paste0("DEA_", prefix, "_", dataset_id, ".csv")
  list(annot = annot, name = name)
}

pad_rows <- function(x, n) {
  x <- x[c(seq_len(nrow(x)), rep(NA_integer_, n - nrow(x))), , drop = FALSE]
  rownames(x) <- NULL
  x
}
