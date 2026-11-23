#' Facets "factor1 + factor2" used by the missingness plots.
#' @param factors named list of explanatory-variable columns
#' @noRd
.factors_facet <- function(factors) {
  ggplot2::facet_grid(cols = ggplot2::vars(!!!rlang::syms(factors)))
}

#' Encode intensities as a binary detection response
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param name column name for binary response
#' @return a tibble of the long-format data with an added integer column `name`: 1 if the response is observed,
#'   0 if it is `NA`
#' @noRd
.encode_bin_resp <- function(lfqdata, name = "bin_resp") {
  pdata <- lfqdata$data_long()
  pdata[[name]] <- as.integer(!is.na(pdata[[lfqdata$response()]]))
  pdata
}

# Functions - Missingness ----

#' Split a contrast "LHS - RHS" into its left/right side expressions, used by
#'
#' get_contrast() to derive group_1 / group_2. Returns NULL for anything that is
#' not a top-level binary minus (get_contrast() rejects those).
#' @param contrast one contrast expression string
#' @noRd
.contrast_sides_expr <- function(contrast) {
  expr <- rlang::parse_expr(contrast)
  if (rlang::is_call(expr, "-") && length(expr) == 3L) {
    list(lhs = expr[[2]], rhs = expr[[3]], full = expr)
  } else {
    NULL
  }
}


#' Histogram summarizing missingness
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param showempty show empty values
#' @param factors factor columns to use
#' @param alpha transparency
#' @return a ggplot object with density plots of the mean intensity per number of missing values, faceted by
#'   factor levels
#' @noRd
.missingness_histogram <- function(lfqdata, showempty = FALSE, factors = lfqdata$relevant_factor_keys(), alpha = 0.1) {
  missing_percent <- .summarize_stats(lfqdata, factor_key = factors)
  missing_percent <- missing_percent |>
    dplyr::ungroup() |>
    dplyr::mutate(nrNAs = as.factor(.data$nrNAs))

  if (showempty) {
    offset <- if (lfqdata$is_transformed()) 1 else 20
    missing_percent <- missing_percent |>
      dplyr::mutate(
        meanAbundance = ifelse(
          is.na(.data$meanAbundance),
          min(.data$meanAbundance, na.rm = TRUE) - offset,
          .data$meanAbundance
        )
      )
  }

  mean_abundance <- paste0("mean_", lfqdata$response())
  missing_percent <- dplyr::rename(missing_percent, !!mean_abundance := "meanAbundance")

  p <- ggplot2::ggplot(
    missing_percent,
    ggplot2::aes(x = .data[[mean_abundance]], fill = .data$nrNAs, colour = .data$nrNAs)
  ) +
    theme_prolfqua() +
    # proteins without any measured value have no mean abundance
    ggplot2::geom_density(alpha = alpha, position = "identity", na.rm = TRUE) +
    .factors_facet(factors) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1))

  if (!lfqdata$is_transformed()) {
    p <- p + ggplot2::scale_x_log10()
  }
  p
}

#' Cumulative counts of missing values per condition
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param factors factor columns to use
#' @return a list with `data`, a wide tibble with one row per condition and one column per number of missing
#'   values holding the cumulative feature count, and `figure`, a ggplot bar chart of these counts
#' @noRd
.missingness_per_condition_cumsum <- function(lfqdata, factors = lfqdata$relevant_factor_keys()) {
  missing_percent <- .summarize_stats(lfqdata, factor_key = factors)

  xx <- missing_percent |>
    group_by(across(all_of(c(factors, "nrNAs", "nrReplicates")))) |>
    dplyr::summarize(nrTransitions = n(), .groups = "drop")

  xxcs <- xx |>
    arrange(.data$nrNAs) |>
    dplyr::mutate(cumulative_sum = cumsum(.data$nrTransitions), .by = all_of(factors))
  res <- xxcs |> dplyr::select(-dplyr::all_of("nrTransitions"))

  nudgeval <- mean(res$cumulative_sum) * 0.05
  p <- ggplot(res, aes(x = .data$nrNAs, y = .data$cumulative_sum)) +
    theme_prolfqua() +
    geom_col(color = "black", fill = "white") +
    geom_text(aes(label = .data$cumulative_sum), nudge_y = nudgeval, angle = -45) +
    .factors_facet(factors)

  res <- res |> tidyr::pivot_wider(names_from = "nrNAs", values_from = "cumulative_sum")
  list(data = res, figure = p)
}

#' Summarize missing in condition as barplot
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param factors factor columns to use
#' @return a list with `data`, a wide tibble with one row per condition and one column per number of missing
#'   values holding the feature count, and `figure`, a ggplot bar chart of these counts
#' @noRd
.missingness_per_condition <- function(lfqdata, factors = lfqdata$relevant_factor_keys()) {
  missing_percent <- .summarize_stats(lfqdata, factor_key = factors)
  hierarchy_key <- tail(lfqdata$hierarchy_keys(), 1)
  hierarchy_key <- paste0("nr_", hierarchy_key)
  xx <- missing_percent |>
    group_by(across(all_of(c(factors, "nrNAs", "nrReplicates")))) |>
    dplyr::summarize(!!hierarchy_key := n(), .groups = "drop")

  nudgeval <- max(xx[[hierarchy_key]]) * 0.05

  p <- ggplot(xx, aes(x = .data$nrNAs, y = .data[[hierarchy_key]])) +
    theme_prolfqua() +
    geom_col(color = "black", fill = "white") +
    geom_text(aes(label = .data[[hierarchy_key]]), nudge_y = nudgeval, angle = 45) +
    .factors_facet(factors)
  xx <- tidyr::pivot_wider(xx, names_from = "nrNAs", values_from = dplyr::all_of(hierarchy_key))

  list(data = xx, figure = p)
}


#' Prepare data for an UpSet plot of observations per group
#'
#' Marks, for each protein and factor level, whether it has at least `tr` observations, in the input format of
#' `UpSetR::upset()`, to show which proteins are observed in which groups.
#'
#' @return a list with `data`, a data.frame with one row per hierarchy element and one 0/1 column per factor
#'   level (1 when at least `tr` observations), and `nsets`, the number of factor-level columns.
#' @export
#' @family plotting
#' @family imputation
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param tr if less than tr observations in condition then missing
#' @examplesIf requireNamespace("UpSetR", quietly = TRUE)
#' istar <- sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' pups <- upset_interaction_missing_stats(lfq)
#' stopifnot(ncol(pups$data) == 5)
#' UpSetR::upset(pups$data, order.by = "freq", nsets = pups$nsets)
upset_interaction_missing_stats <- function(lfqdata, tr = 2) {
  tmp <- .summarize_stats(lfqdata)
  nr_missing <- tmp |>
    tidyr::pivot_wider(
      id_cols = dplyr::all_of(lfqdata$hierarchy_keys()),
      names_from = dplyr::all_of(lfqdata$relevant_factor_keys()),
      values_from = "nrMeasured"
    )

  hl <- length(lfqdata$hierarchy_keys())
  nr_missing[, -seq_len(hl)][nr_missing[, -seq_len(hl)] < tr] <- 0
  nr_missing[, -seq_len(hl)][nr_missing[, -seq_len(hl)] >= tr] <- 1
  list(data = as.data.frame(nr_missing), nsets = ncol(nr_missing) - hl)
}

#' Prepare data for an UpSetR plot of all samples
#'
#' @param lfqdata [LFQData] object with the data and its configuration
#' @return a list with `data`, a data.frame with one row per hierarchy element and one 0/1 presence column per
#'   sample, and `nsets`, the number of sample columns
#' @noRd
.upset_missing_stats <- function(lfqdata) {
  data <- lfqdata$data_long()
  data <- data |>
    dplyr::mutate(
      isThere = dplyr::case_when(
        !is.na(.data[[lfqdata$response()]]) ~ 1,
        TRUE ~ 0
      )
    )
  pups2 <- data |>
    tidyr::pivot_wider(
      id_cols = dplyr::all_of(lfqdata$hierarchy_keys()),
      names_from = dplyr::all_of(lfqdata$sample_name()),
      values_from = "isThere"
    )
  res <- list(data = as.data.frame(pups2), nsets = ncol(pups2) - length(lfqdata$hierarchy_keys()))
  res
}
