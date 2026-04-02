#' convert to binary response
#' @export
#' @keywords internal
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' istar$data <- encode_bin_resp(istar$data, istar$config)
#' istar$config$bin_resp == "bin_resp"
#' istar$data[["bin_resp"]]
encode_bin_resp <- function(pdata, config, name = "bin_resp") {
  config$bin_resp <- "bin_resp"
  pdata <- complete_cases(pdata, config)
  pdata[[config$bin_resp]] <- as.integer(!is.na(pdata[[config$get_response()]]))
  return(pdata)
}

# Functions - Missigness ----

.get_sides <- function(contrast) {
  getAST <- function(ee) purrr::map_if(as.list(ee), is.call, getAST)

  ast_list <- getAST(rlang::parse_expr(contrast))
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
#' config <- istar$config
#' analysis <- istar$data
#' data <- complete_cases(analysis, config)
#'
#' Contrasts <- c("dilution.b-a" = "group_A - group_B", "dilution.c-e" = "group_A - group_Ctrl")
#'
#' var = summarize_stats(data, config)
#' var <- prolfqua::make_interaction_column(var, columns = config$factor_keys_depth())
#'
#' imp <- var |> tidyr::pivot_wider(id_cols = config$hierarchy_keys(),
#'                         names_from = interaction,
#'                         values_from = !!rlang::sym("meanAbundance"))
#'
#' imputed <- get_contrast(imp, config$hierarchy_keys(), Contrasts)
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
#' @export
#' @keywords internal
#' @family plotting
#' @family imputation
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' config <- istar$config
#' analysis <- istar$data
#' xx <- complete_cases(analysis, config)
#' pl <- missigness_histogram(xx, config)
#'
#' pl <- missigness_histogram(analysis, config, showempty=FALSE)
#' stopifnot("ggplot" %in% class(pl))
#' pl <- missigness_histogram(analysis, config, showempty=TRUE)
#' stopifnot("ggplot" %in% class(pl))
#'
missigness_histogram <- function(x, config, showempty = FALSE, factors = config$factor_keys_depth(), alpha = 0.1) {
  missing_percent <- summarize_stats(x, config, factor_key = factors)
  missing_percent <- missing_percent |>
    dplyr::ungroup() |>
    dplyr::mutate(nrNAs = as.factor(.data$nrNAs))

  if (showempty) {
    if (config$is_response_transformed) {
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

  factors <- config$factor_keys_depth()
  formula <- paste(config$isotope_label, "~", paste(factors, collapse = "+"))
  message(formula)
  mean_abundance <- paste0("mean_", config$get_response())
  missing_percent <- dplyr::rename(missing_percent, !!sym(mean_abundance) := .data$meanAbundance)

  p <- ggplot2::ggplot(
    missing_percent,
    ggplot2::aes(x = !!sym(mean_abundance), fill = .data$nrNAs, colour = .data$nrNAs)
  ) +
    ggplot2::geom_density(alpha = alpha, position = "identity") +
    ggplot2::facet_grid(as.formula(formula)) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1))

  if (!config$is_response_transformed) {
    p <- p + ggplot2::scale_x_log10()
  }
  p
}

#' cumulative sums of missing
#' @export
#' @keywords internal
#' @family plotting
#' @family imputation
#' @examples
#'
#'
#' istar <- sim_lfq_data_peptide_config()
#' config <- istar$config
#' analysis <- istar$data
#'
#' res <- missingness_per_condition_cumsum(analysis,config)
#' stopifnot("ggplot" %in% class(res$figure))
#' stopifnot(ncol(res$data) >= 6)
#'
missingness_per_condition_cumsum <- function(x, config, factors = config$factor_keys_depth()) {
  missing_percent <- summarize_stats(x, config, factor_key = factors)

  xx <- missing_percent |>
    group_by(across(all_of(c(config$isotope_label, factors, "nrNAs", "nrReplicates")))) |>
    dplyr::summarize(nrTransitions = n())

  xxcs <- xx |>
    group_by(across(all_of(c(config$isotope_label, factors)))) |>
    arrange(.data$nrNAs) |>
    dplyr::mutate(cumulative_sum = cumsum(.data$nrTransitions))
  res <- xxcs |> dplyr::select(-dplyr::all_of("nrTransitions"))

  formula <- paste(config$isotope_label, "~", paste(factors, collapse = "+"))
  message(formula)

  nudgeval <- mean(res$cumulative_sum) * 0.05
  p <- ggplot(res, aes(x = .data$nrNAs, y = .data$cumulative_sum)) +
    geom_bar(stat = "identity", color = "black", fill = "white") +
    geom_text(aes(label = .data$cumulative_sum), nudge_y = nudgeval, angle = -45) +
    facet_grid(as.formula(formula))

  res <- res |> tidyr::pivot_wider(names_from = "nrNAs", values_from = "cumulative_sum")
  return(list(data = res, figure = p))
}

#' Summarize missing in condtion as barplot
#' @export
#' @keywords internal
#' @family plotting
#' @family imputation
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' config <- istar$config
#' analysis <- istar$data
#'
#' res <- missingness_per_condition(analysis, config)
#' stopifnot("ggplot" %in% class(res$figure))
#'
#' stopifnot(ncol(res$data) >= 5)
#'
missingness_per_condition <- function(x, config, factors = config$factor_keys_depth()) {
  missing_percent <- summarize_stats(x, config, factor_key = factors)
  hierarchy_key <- tail(config$hierarchy_keys(), 1)
  hierarchy_key <- paste0("nr_", hierarchy_key)
  xx <- missing_percent |>
    group_by(across(all_of(c(config$isotope_label, factors, "nrNAs", "nrReplicates")))) |>
    dplyr::summarize(!!sym(hierarchy_key) := n())

  formula <- paste(config$isotope_label, "~", paste(factors, collapse = "+"))

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
#' @export
#' @family plotting
#' @family imputation
#' @param data data.frame
#' @param cf AnalysisConfiguration
#' @param tr if less than tr observations in condition then missing
#' @examples
#' istar <- sim_lfq_data_peptide_config()
#' config <- istar$config
#' analysis <- istar$data
#'
#' pups <- UpSet_interaction_missing_stats(analysis, config)
#' stopifnot(ncol(pups$data) == 5)
#' UpSetR::upset(pups$data, order.by = "freq", nsets = pups$nsets)
UpSet_interaction_missing_stats <- function(data, cf, tr = 2) {
  tmp <- prolfqua::summarize_stats(data, cf)
  nr_missing <- tmp |>
    tidyr::pivot_wider(
      id_cols = cf$hierarchy_keys(),
      names_from = cf$factor_keys_depth(),
      values_from = !!rlang::sym("nrMeasured")
    )

  hl <- length(cf$hierarchy_keys())
  nr_missing[, -(1:hl)][nr_missing[, -(1:hl)] < tr] <- 0
  nr_missing[, -(1:hl)][nr_missing[, -(1:hl)] >= tr] <- 1
  return(list(data = as.data.frame(nr_missing), nsets = ncol(nr_missing) - length(cf$hierarchy_keys())))
}

#' prepare dataframe for UpSetR plot for all samples
#'
#' @export
#' @keywords internal
#' @family plotting
#' @family imputation
#' @examples
#' istar <- sim_lfq_data_peptide_config()
#' config <- istar$config
#' analysis <- istar$data
#' pups <- UpSet_missing_stats(analysis, config)
#' UpSetR::upset(pups$data , order.by = "freq", nsets = pups$nsets)
UpSet_missing_stats <- function(data, config) {
  data <- prolfqua::complete_cases(data, config)
  data <- data |>
    dplyr::mutate(
      isThere = dplyr::case_when(
        !is.na(!!rlang::sym(config$get_response())) ~ 1,
        TRUE ~ 0
      )
    )
  pups2 <- data |>
    tidyr::pivot_wider(
      id_cols = config$hierarchy_keys(),
      names_from = config$sample_name,
      values_from = !!rlang::sym("isThere")
    )
  res <- list(data = as.data.frame(pups2), nsets = ncol(pups2) - length(config$hierarchy_keys()))
  return(res)
}
