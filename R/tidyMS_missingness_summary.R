#' convert to binary response
#' @param lfqdata LFQData object
#' @param name column name for binary response
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

# Functions - Missigness ----

.get_sides <- function(contrast) {
  get_ast <- function(ee) purrr::map_if(as.list(ee), is.call, get_ast)

  ast_list <- get_ast(rlang::parse_expr(contrast))
  ast_array <- array(as.character(unlist(ast_list)))
  ast_array <- gsub("`", "", ast_array)
  return(ast_array)
}

#' Compute fold changes given Contrasts
#' @keywords internal
#' @family imputation
#' @param data data.frame
#' @param data hierarchy_keys of Analysis Configuration
#' @param contrasts list of contrasts to compute
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
  for (i in seq_along(contrasts)) {
    message(names(contrasts)[i], "=", contrasts[i], "\n")
    data <- dplyr::mutate(data, !!names(contrasts)[i] := !!rlang::parse_expr(contrasts[i]))
  }
  res <- vector(mode = "list", length(contrasts))
  names(res) <- names(contrasts)
  for (i in seq_along(contrasts)) {
    sides <- .get_sides(contrasts[i])
    sides <- intersect(sides, colnames(data))

    df <- dplyr::select(
      data,
      dplyr::all_of(hierarchy_keys),
      group_1 = dplyr::all_of(sides[1]),
      group_2 = dplyr::all_of(sides[2]),
      estimate = dplyr::all_of(names(contrasts)[i])
    )

    df$group_1_name <- sides[1]
    df$group_2_name <- sides[2]
    df$contrast <- names(contrasts)[i]

    res[[names(contrasts)[i]]] <- df
  }
  res <- dplyr::bind_rows(res)
  return(dplyr::ungroup(res))
}


#' Histogram summarizing missigness
#' @param lfqdata LFQData object
#' @param showempty show empty values
#' @param factors factor columns to use
#' @param alpha transparency
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
    if (lfqdata$is_transformed()) {
      missing_percent <- missing_percent |>
        dplyr::mutate(
          meanAbundance = ifelse(
            is.na(.data$meanAbundance),
            min(.data$meanAbundance, na.rm = TRUE) - 1,
            .data$meanAbundance
          )
        )
    } else {
      missing_percent <- missing_percent |>
        dplyr::mutate(
          meanAbundance = ifelse(
            is.na(.data$meanAbundance),
            min(.data$meanAbundance, na.rm = TRUE) - 20,
            .data$meanAbundance
          )
        )
    }
  }

  factors <- lfqdata$relevant_factor_keys()
  formula <- paste(lfqdata$isotope_label(), "~", paste(factors, collapse = "+"))
  message(formula)
  mean_abundance <- paste0("mean_", lfqdata$response())
  missing_percent <- dplyr::rename(missing_percent, !!sym(mean_abundance) := .data$meanAbundance)

  p <- ggplot2::ggplot(
    missing_percent,
    ggplot2::aes(x = !!sym(mean_abundance), fill = .data$nrNAs, colour = .data$nrNAs)
  ) +
    ggplot2::geom_density(alpha = alpha, position = "identity") +
    ggplot2::facet_grid(as.formula(formula)) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1))

  if (!lfqdata$is_transformed()) {
    p <- p + ggplot2::scale_x_log10()
  }
  p
}

#' cumulative sums of missing
#' @param lfqdata LFQData object
#' @param factors factor columns to use
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
    dplyr::summarize(nrTransitions = n())

  xxcs <- xx |>
    group_by(across(all_of(c(lfqdata$isotope_label(), factors)))) |>
    arrange(.data$nrNAs) |>
    dplyr::mutate(cumulative_sum = cumsum(.data$nrTransitions))
  res <- xxcs |> dplyr::select(-dplyr::all_of("nrTransitions"))

  formula <- paste(lfqdata$isotope_label(), "~", paste(factors, collapse = "+"))
  message(formula)

  nudgeval <- mean(res$cumulative_sum) * 0.05
  p <- ggplot(res, aes(x = .data$nrNAs, y = .data$cumulative_sum)) +
    geom_bar(stat = "identity", color = "black", fill = "white") +
    geom_text(aes(label = .data$cumulative_sum), nudge_y = nudgeval, angle = -45) +
    facet_grid(as.formula(formula))

  res <- res |> tidyr::pivot_wider(names_from = "nrNAs", values_from = "cumulative_sum")
  return(list(data = res, figure = p))
}

#' Summarize missing in condition as barplot
#' @param lfqdata LFQData object
#' @param factors factor columns to use
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
    dplyr::summarize(!!sym(hierarchy_key) := n())

  formula <- paste(lfqdata$isotope_label(), "~", paste(factors, collapse = "+"))

  nudgeval <- max(xx[[hierarchy_key]]) * 0.05

  p <- ggplot(xx, aes(x = .data$nrNAs, y = .data[[hierarchy_key]])) +
    geom_bar(stat = "identity", color = "black", fill = "white") +
    geom_text(aes(label = !!sym(hierarchy_key)), nudge_y = nudgeval, angle = 45) +
    facet_grid(as.formula(formula))
  xx <- tidyr::pivot_wider(xx, names_from = "nrNAs", values_from = hierarchy_key)

  return(list(data = xx, figure = p))
}


#' UpSetR plot from interaction_missing_stats
#'
#' @return The computed result.
#' @export
#' @family plotting
#' @family imputation
#' @param lfqdata LFQData object
#' @param tr if less than tr observations in condition then missing
#' @examples
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

#' prepare dataframe for UpSetR plot for all samples
#'
#' @param lfqdata LFQData object
#' @export
#' @keywords internal
#' @family plotting
#' @family imputation
#' @examples
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
