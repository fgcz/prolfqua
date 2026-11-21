#' Facet formula "isotope label ~ factor1 + factor2" used by the missingness plots.
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param factors named list of explanatory-variable columns
#' @noRd
.isotope_by_factors_formula <- function(lfqdata, factors) {
  paste(lfqdata$isotope_label(), "~", paste(factors, collapse = "+"))
}

#' Encode intensities as a binary detection response
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param name column name for binary response
#' @return a tibble of the long-format data with an added integer column `name`: 1 if the response is observed,
#'   0 if it is `NA`
#' @export
#' @keywords internal
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' lfq$set_data(encode_bin_resp(lfq))
#' lfq$set_config_value("bin_resp", "bin_resp")
#' lfq$data_long()[["bin_resp"]]
encode_bin_resp <- function(lfqdata, name = "bin_resp") {
  pdata <- lfqdata$data_long()
  pdata[[name]] <- as.integer(!is.na(pdata[[lfqdata$response()]]))
  return(pdata)
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

#' Compute fold changes given Contrasts
#' @keywords internal
#' @family imputation
#' @param data data.frame
#' @param hierarchy_keys hierarchy key columns of the AnalysisConfiguration
#' @param contrasts list of contrasts to compute
#' @return a tibble with one row per hierarchy element and contrast and columns of the `hierarchy_keys`,
#'   `group_1`, `group_2`, `estimate`, `group_1_name`, `group_2_name` and `contrast`
#' @export
#' @examples
#'
#'
#' istar <- sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' lfq$complete_cases()
#'
#' Contrasts <- c("dilution.b-a" = "group_A - group_B", "dilution.c-e" = "group_A - group_Ctrl")
#'
#' var <- summarize_stats(lfq)
#' var <- prolfqua::make_interaction_column(var, columns = lfq$relevant_factor_keys())
#'
#' imp <- var |> tidyr::pivot_wider(id_cols = lfq$hierarchy_keys(),
#'                         names_from = interaction,
#'                         values_from = !!rlang::sym("meanAbundance"))
#'
#' imputed <- get_contrast(imp, lfq$hierarchy_keys(), Contrasts)
#'
#'
get_contrast <- function(data, hierarchy_keys, contrasts) {
  res <- vector(mode = "list", length(contrasts))
  names(res) <- names(contrasts)
  for (i in seq_along(contrasts)) {
    cname <- names(contrasts)[i]

    # A contrast must be a difference "LHS - RHS" (e.g. "group_A - group_B" or
    # "(group_A + group_B)/2 - group_Ctrl"). group_1 / group_2 are the evaluated
    # left / right sides; estimate is their difference.
    sx <- .contrast_sides_expr(contrasts[i])
    if (is.null(sx)) {
      abort_prolfqua(
        "get_contrast: contrast '",
        cname,
        "' (",
        contrasts[i],
        ") is not of the required form 'LHS - RHS'. A contrast must be a ",
        "difference of two group expressions, e.g. 'group_A - group_B' or ",
        "'(group_A + group_B)/2 - group_Ctrl'."
      )
    }

    # A contrast referencing a level absent from `data` errors here at mutate
    # ("object '<level>' not found") -- once per contrast, not per row.
    dd <- dplyr::mutate(data, group_1 = !!sx$lhs, group_2 = !!sx$rhs, estimate = !!sx$full)
    # Materialize this contrast under its own name so a *later* nested contrast
    # can reference it as a group, e.g. "T_C_gv_WT - T_C_gv_KO" after
    # T_C_gv_WT / T_C_gv_KO have been defined. Sides may therefore be a level
    # present in `data` or the name of an earlier contrast.
    data[[cname]] <- dd[["estimate"]]
    df <- dplyr::select(dd, dplyr::all_of(hierarchy_keys), "group_1", "group_2", "estimate")
    df$group_1_name <- rlang::as_label(sx$lhs)
    df$group_2_name <- rlang::as_label(sx$rhs)
    df$contrast <- cname
    res[[cname]] <- df
  }
  dplyr::ungroup(dplyr::bind_rows(res))
}


#' Histogram summarizing missingness
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param showempty show empty values
#' @param factors factor columns to use
#' @param alpha transparency
#' @return a ggplot object with density plots of the mean intensity per number of missing values, faceted by
#'   factor levels
#' @export
#' @keywords internal
#' @family plotting
#' @family imputation
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' pl <- missigness_histogram(lfq, showempty = FALSE)
#' stopifnot("ggplot" %in% class(pl))
#' pl <- missigness_histogram(lfq, showempty = TRUE)
#' stopifnot("ggplot" %in% class(pl))
#'
missigness_histogram <- function(lfqdata, showempty = FALSE, factors = lfqdata$relevant_factor_keys(), alpha = 0.1) {
  missing_percent <- summarize_stats(lfqdata, factor_key = factors)
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

  factors <- lfqdata$relevant_factor_keys()
  formula <- .isotope_by_factors_formula(lfqdata, factors)
  mean_abundance <- paste0("mean_", lfqdata$response())
  missing_percent <- dplyr::rename(missing_percent, !!sym(mean_abundance) := "meanAbundance")

  p <- ggplot2::ggplot(
    missing_percent,
    ggplot2::aes(x = !!sym(mean_abundance), fill = .data$nrNAs, colour = .data$nrNAs)
  ) +
    theme_prolfqua() +
    ggplot2::geom_density(alpha = alpha, position = "identity") +
    ggplot2::facet_grid(as.formula(formula)) +
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
#' @export
#' @keywords internal
#' @family plotting
#' @family imputation
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' res <- missingness_per_condition_cumsum(lfq)
#' stopifnot("ggplot" %in% class(res$figure))
#' stopifnot(ncol(res$data) >= 6)
#'
missingness_per_condition_cumsum <- function(lfqdata, factors = lfqdata$relevant_factor_keys()) {
  missing_percent <- summarize_stats(lfqdata, factor_key = factors)

  xx <- missing_percent |>
    group_by(across(all_of(c(lfqdata$isotope_label(), factors, "nrNAs", "nrReplicates")))) |>
    dplyr::summarize(nrTransitions = n(), .groups = "drop")

  xxcs <- xx |>
    group_by(across(all_of(c(lfqdata$isotope_label(), factors)))) |>
    arrange(.data$nrNAs) |>
    dplyr::mutate(cumulative_sum = cumsum(.data$nrTransitions)) |>
    dplyr::ungroup()
  res <- xxcs |> dplyr::select(-dplyr::all_of("nrTransitions"))

  formula <- .isotope_by_factors_formula(lfqdata, factors)

  nudgeval <- mean(res$cumulative_sum) * 0.05
  p <- ggplot(res, aes(x = .data$nrNAs, y = .data$cumulative_sum)) +
    theme_prolfqua() +
    geom_col(color = "black", fill = "white") +
    geom_text(aes(label = .data$cumulative_sum), nudge_y = nudgeval, angle = -45) +
    facet_grid(as.formula(formula))

  res <- res |> tidyr::pivot_wider(names_from = "nrNAs", values_from = "cumulative_sum")
  return(list(data = res, figure = p))
}

#' Summarize missing in condition as barplot
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param factors factor columns to use
#' @return a list with `data`, a wide tibble with one row per condition and one column per number of missing
#'   values holding the feature count, and `figure`, a ggplot bar chart of these counts
#' @export
#' @keywords internal
#' @family plotting
#' @family imputation
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' res <- missingness_per_condition(lfq)
#' stopifnot("ggplot" %in% class(res$figure))
#' stopifnot(ncol(res$data) >= 5)
#'
missingness_per_condition <- function(lfqdata, factors = lfqdata$relevant_factor_keys()) {
  missing_percent <- summarize_stats(lfqdata, factor_key = factors)
  hierarchy_key <- tail(lfqdata$hierarchy_keys(), 1)
  hierarchy_key <- paste0("nr_", hierarchy_key)
  xx <- missing_percent |>
    group_by(across(all_of(c(lfqdata$isotope_label(), factors, "nrNAs", "nrReplicates")))) |>
    dplyr::summarize(!!sym(hierarchy_key) := n(), .groups = "drop")

  formula <- .isotope_by_factors_formula(lfqdata, factors)

  nudgeval <- max(xx[[hierarchy_key]]) * 0.05

  p <- ggplot(xx, aes(x = .data$nrNAs, y = .data[[hierarchy_key]])) +
    theme_prolfqua() +
    geom_col(color = "black", fill = "white") +
    geom_text(aes(label = !!sym(hierarchy_key)), nudge_y = nudgeval, angle = 45) +
    facet_grid(as.formula(formula))
  xx <- tidyr::pivot_wider(xx, names_from = "nrNAs", values_from = dplyr::all_of(hierarchy_key))

  return(list(data = xx, figure = p))
}


#' UpSetR plot from interaction_missing_stats
#'
#' @return The computed result.
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
  tmp <- prolfqua::summarize_stats(lfqdata)
  nr_missing <- tmp |>
    tidyr::pivot_wider(
      id_cols = lfqdata$hierarchy_keys(),
      names_from = lfqdata$relevant_factor_keys(),
      values_from = !!rlang::sym("nrMeasured")
    )

  hl <- length(lfqdata$hierarchy_keys())
  nr_missing[, -seq_len(hl)][nr_missing[, -seq_len(hl)] < tr] <- 0
  nr_missing[, -seq_len(hl)][nr_missing[, -seq_len(hl)] >= tr] <- 1
  return(list(data = as.data.frame(nr_missing), nsets = ncol(nr_missing) - hl))
}

#' Prepare data for an UpSetR plot of all samples
#'
#' @param lfqdata [LFQData] object with the data and its configuration
#' @return a list with `data`, a data.frame with one row per hierarchy element and one 0/1 presence column per
#'   sample, and `nsets`, the number of sample columns
#' @export
#' @keywords internal
#' @family plotting
#' @family imputation
#' @examplesIf requireNamespace("UpSetR", quietly = TRUE)
#' istar <- sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' pups <- upset_missing_stats(lfq)
#' UpSetR::upset(pups$data, order.by = "freq", nsets = pups$nsets)
upset_missing_stats <- function(lfqdata) {
  data <- lfqdata$data_long()
  data <- data |>
    dplyr::mutate(
      isThere = dplyr::case_when(
        !is.na(!!rlang::sym(lfqdata$response())) ~ 1,
        TRUE ~ 0
      )
    )
  pups2 <- data |>
    tidyr::pivot_wider(
      id_cols = lfqdata$hierarchy_keys(),
      names_from = lfqdata$sample_name(),
      values_from = !!rlang::sym("isThere")
    )
  res <- list(data = as.data.frame(pups2), nsets = ncol(pups2) - length(lfqdata$hierarchy_keys()))
  return(res)
}
