#' group label function
#' @param primary character, primary factor level
#' @param secondary character, secondary factor level
#' @return A character vector with combined group labels.
#' @export
#' @family modelling
#' @examples
#' group_label("a", "b")
group_label <- function(primary, secondary) paste0("G_", primary, "_", secondary)


# Generic function to build pairwise contrasts for any two factors
# Main effect contrasts (averaged across all secondary levels)
#' main effects contrasts
#' @param primary_levels character vector of primary factor levels
#' @param secondary_levels character vector of secondary factor levels
#' @return Contrast definitions or contrast results.
#' @export
#' @family contrasts
#' @examples
#' primary_levels <- c("MI", "MINOCAM")
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

#' Level-specific contrasts (per secondary level)
#' @param primary_levels character vector of primary factor levels
#' @param secondary_levels character vector of secondary factor levels
#' @return Contrast definitions or contrast results.
#' @export
#' @family contrasts
#' @examples
#' # example code
#'
#' primary_levels <- c("MI", "MINOCAM")
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


#' Interaction contrasts (difference of differences)
#' @param primary_levels character vector of primary factor levels
#' @param secondary_levels character vector of secondary factor levels
#' @return Contrast definitions or contrast results.
#' @export
#' @family contrasts
#' @examples
#' primary_levels <- c("MI", "MINOCAM")
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

#' Single-factor contrasts (pairwise comparisons)
#' @param levels character vector of factor levels
#' @return Contrast definitions or contrast results.
#' @export
#' @family contrasts
#' @examples
#' # example code
#' primary_levels <- c("MI", "MINOCAM")
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

# Combined generate_contrasts

#' Combined generate_contrasts
#' @param primary_levels character vector of primary factor levels
#' @param secondary_levels character vector of secondary factor levels
#' @param interactions logical, if TRUE include interaction contrasts
#' @return Contrast definitions or contrast results.
#' @export
#' @family contrasts
#' @examples
#' primary_levels <- c("MI", "MINOCAM")
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
  return(res)
}

#' DRY function: process and export annotated contrasts
#' @param df data.frame with annotation data
#' @param primary_col character, name of the primary factor column
#' @param secondary_col character, name of the secondary factor column
#' @param prefix character, prefix for output naming (default "G_")
#' @param dataset_id character, dataset identifier for output naming
#' @param decreasing logical, if TRUE sort factor levels in decreasing order
#' @param interactions logical, if TRUE include interaction contrasts
#' @return Contrast definitions or contrast results.
#' @export
#' @family contrasts
#' @examples
#'
#' data("x5463yzwer453bbb", package = "prolfqua")
#' annotation_add_contrasts(x5463yzwer453bbb, "factor_A", "factor_B", "primary")
#' annotation_add_contrasts(x5463yzwer453bbb, "factor_B", "factor_A", "secondary")
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

  df <- df |> tidyr::unite("Group", primary_col, secondary_col, sep = "_", remove = FALSE)
  df_contrast <- generate_contrasts(levels_p, levels_s, interactions = interactions)

  if (nrow(df_contrast) < nrow(df)) {
    pad <- tibble::tibble(
      ContrastName = rep(NA_character_, nrow(df) - nrow(df_contrast)),
      Contrast = rep(NA_character_, nrow(df) - nrow(df_contrast))
    )
    df_contrast <- bind_rows(df_contrast, pad)
  }
  # Annotate and export
  annot <- dplyr::bind_cols(
    df %>% dplyr::select(-dplyr::all_of(c(primary_col, secondary_col))),
    df_contrast
  )
  name <- paste0("DEA_", prefix, "_", dataset_id, ".csv")
  return(list(annot = annot, name = name))
}
