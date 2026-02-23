#' group label function
#' @param primary character, primary factor level
#' @param secondary character, secondary factor level
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
      "( (%s)/%d - (%s)/%d )",
      avg_high, length(secondary_levels),
      avg_low,  length(secondary_levels)
    )
  }
  contrasts
}

#' Level-specific contrasts (per secondary level)
#' @param primary_levels character vector of primary factor levels
#' @param secondary_levels character vector of secondary factor levels
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
        group_label(high, s2), group_label(low, s2),
        group_label(high, s1), group_label(low, s1)
      )
    }
  }
  contrasts
}

#' Single-factor contrasts (pairwise comparisons)
#' @param levels character vector of factor levels
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
    if(interactions){unlist(interaction_contrasts(primary_levels, secondary_levels))} else {NULL}
  )
  res <- data.frame(ContrastName = names(res), Contrast = res)
  return(res)
}

x5463yzwer453bbb <- structure(
  list(
    `Relative Path` = c(
      "20250429_019_C38312-12_S935004_MI_150_6_S1-D2.d.zip",
      "20250429_003_C38312-21_S935013_MINOCA_0_3_S1-E3.d.zip",
      "20250429_015_C38312-35_S935027_MINOCA_300_5_S1-C5.d.zip",
      "20250429_016_C38312-11_S935003_MI_150_5_S1-C2.d.zip",
      "20250429_014_C38312-7_S934999_MI_150_1_S1-G1.d.zip",
      "20250429_035_C38312-3_S934995_MI_0_3_S1-C1.d.zip", "20250429_022_C38312-25_S935017_MINOCA_150_1_S1-A4.d.zip",
      "20250429_043_C38312-14_S935006_MI_300_2_S1-F2.d.zip", "20250429_004_C38312-29_S935021_MINOCA_150_5_S1-E4.d.zip",
      "20250429_034_C38312-15_S935007_MI_300_3_S1-G2.d.zip", "20250429_017_C38312-8_S935000_MI_150_2_S1-H1.d.zip",
      "20250429_032_C38312-13_S935005_MI_300_1_S1-E2.d.zip", "20250429_028_C38312-23_S935015_MINOCA_0_5_S1-G3.d.zip",
      "20250429_027_C38312-34_S935026_MINOCA_300_4_S1-B5.d.zip", "20250429_009_C38312-4_S934996_MI_0_4_S1-D1.d.zip",
      "20250429_044_C38312-32_S935024_MINOCA_300_2_S1-H4.d.zip",
      "20250429_036_C38312-24_S935016_MINOCA_0_6_S1-H3.d.zip",
      "20250429_012_C38312-1_S934993_MI_0_1_S1-A1.d.zip", "20250429_006_C38312-9_S935001_MI_150_3_S1-A2.d.zip",
      "20250429_005_C38312-18_S935010_MI_300_6_S1-B3.d.zip", "20250429_031_C38312-16_S935008_MI_300_4_S1-H2.d.zip",
      "20250429_026_C38312-19_S935011_MINOCA_0_1_S1-C3.d.zip",
      "20250429_010_C38312-28_S935020_MINOCA_150_4_S1-D4.d.zip",
      "20250429_038_C38312-20_S935012_MINOCA_0_2_S1-D3.d.zip", "20250429_033_C38312-17_S935009_MI_300_5_S1-A3.d.zip",
      "20250429_037_C38312-26_S935018_MINOCA_150_2_S1-B4.d.zip", "20250429_023_C38312-5_S934997_MI_0_5_S1-E1.d.zip",
      "20250429_018_C38312-33_S935025_MINOCA_300_3_S1-A5.d.zip",
      "20250429_041_C38312-36_S935028_MINOCA_300_6_S1-D5.d.zip",
      "20250429_029_C38312-30_S935022_MINOCA_150_6_S1-F4.d.zip",
      "20250429_024_C38312-31_S935023_MINOCA_300_1_S1-G4.d.zip",
      "20250429_042_C38312-2_S934994_MI_0_2_S1-B1.d.zip", "20250429_013_C38312-22_S935014_MINOCA_0_4_S1-F3.d.zip",
      "20250429_007_C38312-6_S934998_MI_0_6_S1-F1.d.zip", "20250429_025_C38312-10_S935002_MI_150_4_S1-B2.d.zip",
      "20250429_008_C38312-27_S935019_MINOCA_150_3_S1-C4.d.zip"
    ), Name = c(
      "MI_150_6",
      "MINOCA_0_3", "MINOCA_300_5", "MI_150_5", "MI_150_1", "MI_0_3",
      "MINOCA_150_1", "MI_300_2", "MINOCA_150_5", "MI_300_3", "MI_150_2",
      "MI_300_1", "MINOCA_0_5", "MINOCA_300_4", "MI_0_4", "MINOCA_300_2",
      "MINOCA_0_6", "MI_0_1", "MI_150_3", "MI_300_6", "MI_300_4", "MINOCA_0_1",
      "MINOCA_150_4", "MINOCA_0_2", "MI_300_5", "MINOCA_150_2", "MI_0_5",
      "MINOCA_300_3", "MINOCA_300_6", "MINOCA_150_6", "MINOCA_300_1",
      "MI_0_2", "MINOCA_0_4", "MI_0_6", "MI_150_4", "MINOCA_150_3"
    ),
    Group = c(
      "MI_T150", "MINOCA_T0", "MINOCA_T300", "MI_T150",
      "MI_T150", "MI_T0", "MINOCA_T150", "MI_T300", "MINOCA_T150",
      "MI_T300", "MI_T150", "MI_T300", "MINOCA_T0", "MINOCA_T300",
      "MI_T0", "MINOCA_T300", "MINOCA_T0", "MI_T0", "MI_T150",
      "MI_T300", "MI_T300", "MINOCA_T0", "MINOCA_T150", "MINOCA_T0",
      "MI_T300", "MINOCA_T150", "MI_T0", "MINOCA_T300", "MINOCA_T300",
      "MINOCA_T150", "MINOCA_T300", "MI_T0", "MINOCA_T0", "MI_T0",
      "MI_T150", "MINOCA_T150"
    ), factor_A = c(
      "MI", "MINOCA", "MINOCA",
      "MI", "MI", "MI", "MINOCA", "MI", "MINOCA", "MI", "MI", "MI",
      "MINOCA", "MINOCA", "MI", "MINOCA", "MINOCA", "MI", "MI",
      "MI", "MI", "MINOCA", "MINOCA", "MINOCA", "MI", "MINOCA",
      "MI", "MINOCA", "MINOCA", "MINOCA", "MINOCA", "MI", "MINOCA",
      "MI", "MI", "MINOCA"
    ), factor_B = c(
      "T150", "T0", "T300",
      "T150", "T150", "T0", "T150", "T300", "T150", "T300", "T150",
      "T300", "T0", "T300", "T0", "T300", "T0", "T0", "T150", "T300",
      "T300", "T0", "T150", "T0", "T300", "T150", "T0", "T300",
      "T300", "T150", "T300", "T0", "T0", "T0", "T150", "T150"
    )
  ),
  row.names = c(
    NA,
    -36L
  ), class = c("tbl_df", "tbl", "data.frame")
)

#' DRY function: process and export annotated contrasts
#' @param df data.frame with annotation data
#' @param primary_col character, name of the primary factor column
#' @param secondary_col character, name of the secondary factor column
#' @param prefix character, prefix for output naming (default "G_")
#' @param dataset_id character, dataset identifier for output naming
#' @param decreasing logical, if TRUE sort factor levels in decreasing order
#' @param interactions logical, if TRUE include interaction contrasts
#' @export
#' @family contrasts
#' @examples
#'
#' annotation_add_contrasts(prolfqua::x5463yzwer453bbb, "factor_A", "factor_B", "primary")
#' annotation_add_contrasts(prolfqua::x5463yzwer453bbb, "factor_B", "factor_A", "secondary")
annotation_add_contrasts <- function(
    df,
    primary_col,
    secondary_col,
    prefix = "G_",
    dataset_id = "dataset",
    decreasing = FALSE,
    interactions = TRUE) {
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
