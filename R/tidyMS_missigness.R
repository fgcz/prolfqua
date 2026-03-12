#' convert to binary response
#' @export
#' @keywords internal
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' istar$data <- encode_bin_resp(istar$data, istar$config)
#' istar$config$bin_resp == "bin_resp"
#' istar$data[["bin_resp"]]
encode_bin_resp <- function(pdata, config, name = "bin_resp"){
  config$bin_resp = "bin_resp"
  pdata <- complete_cases(pdata,config)
  pdata[[config$bin_resp]] <- as.integer(!is.na(pdata[[config$get_response()]]))
  return(pdata)
}

# Functions - Missigness ----

#' compute missingness statistics per hierarchy and factor level
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param factors factor to include (default up to factor depth)
#' @param hierarchy hierarchy to include (default up to hierarchy depth)
#' @param workIntensity work intensity column
#' @export
#' @keywords internal
#' @examples
#'
#'
#' istar <- sim_lfq_data_peptide_config()
#' config <- istar$config
#' analysis <- istar$data
#'
#' xx <- complete_cases(analysis, config)
#' x <- interaction_missing_stats(xx, config)$data |> dplyr::arrange(desc(nrNAs))
#' nrow(x)
#' tmp <- interaction_missing_stats(xx, config,
#'  factors= character(),
#'   hierarchy = config$hierarchy_keys()[1])$data
#' stopifnot(nrow(tmp) == 10)
#' tmp <- interaction_missing_stats(xx, config,
#'   hierarchy = config$hierarchy_keys()[1])$data
#' stopifnot(nrow(tmp) == length(unique(xx$protein_Id))* length(unique(xx$group_)))
#' stopifnot(sum(is.na(tmp$nrMeasured))==0)
#'
#' tmp <- interaction_missing_stats(xx, config, factors = NULL)
interaction_missing_stats <- function(pdata,
                                      config,
                                      factors = config$factor_keys_depth(),
                                      hierarchy = config$hierarchy_keys(),
                                      workIntensity = config$get_response())
{
  warning(">>>> deprecated! <<<< \n
          use summarize_stats_factors instead.")
  pdata <- complete_cases(pdata, config)
  missingPrec <- pdata |> group_by(across(all_of(c(factors,
                                        hierarchy,
                                        config$isotopeLabel
  ))))
  missingPrec <- missingPrec |>
    dplyr::summarize(nrReplicates = n(),
                     nrNAs = sum(is.na(!!sym(workIntensity))),
                     meanAbundance = mean(!!sym(workIntensity), na.rm = TRUE),
                     medianAbundance = median(!!sym(workIntensity), na.rm = TRUE)) |>
    dplyr::mutate(nrMeasured = .data$nrReplicates - .data$nrNAs) |> dplyr::ungroup()
  return(list(data = missingPrec,
              summaries = c("nrReplicates","nrNAs","nrMeasured","meanAbundance", "medianAbundance")))
}







.get_sides <- function(contrast) {
  getAST <- function(ee) purrr::map_if(as.list(ee), is.call, getAST)

  ast_list <- getAST(rlang::parse_expr(contrast))
  ast_array <- array(as.character(unlist(ast_list)))
  ast_array <- gsub("`","",ast_array)
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
get_contrast <- function(data,
                         hierarchy_keys,
                         contrasts)
{
  for (i in seq_along(contrasts)) {
    message(names(contrasts)[i], "=", contrasts[i],"\n")
    data <- dplyr::mutate(data, !!names(contrasts)[i] := !!rlang::parse_expr(contrasts[i]))
  }
  res <- vector(mode = "list", length(contrasts))
  names(res) <- names(contrasts)
  for (i in seq_along(contrasts)) {
    sides <- .get_sides(contrasts[i] )
    sides <- intersect(sides,colnames(data))

    df  <- dplyr::select(
      data,
      dplyr::all_of(hierarchy_keys),
      group_1 = dplyr::all_of(sides[1]),
      group_2 = dplyr::all_of(sides[2]),
      estimate = dplyr::all_of(names(contrasts)[i])
    )

    df$group_1_name <- sides[1]
    df$group_2_name <- sides[2]
    df$contrast <-  names(contrasts)[i]

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
missigness_histogram <- function(x,
                                 config,
                                 showempty = FALSE,
                                 factors = config$factor_keys_depth(),
                                 alpha = 0.1){
  missingPrec <- interaction_missing_stats(x, config , factors)$data
  missingPrec <- missingPrec |>
    dplyr::ungroup() |> dplyr::mutate(nrNAs = as.factor(.data$nrNAs))

  if (showempty) {
    if (config$is_response_transformed) {
      missingPrec <- missingPrec |>
        dplyr::mutate(meanAbundance = ifelse(is.na(.data$meanAbundance), min(.data$meanAbundance, na.rm = TRUE) - 1,
                                             .data$meanAbundance))
    }else{
      missingPrec <- missingPrec |>
        dplyr::mutate(meanAbundance = ifelse(is.na(.data$meanAbundance),min(.data$meanAbundance, na.rm = TRUE) - 20,.data$meanAbundance))
    }

  }

  factors <- config$factor_keys_depth()
  formula <- paste(config$isotopeLabel, "~", paste(factors, collapse = "+"))
  message(formula)
  meanAbundance <- paste0("mean_", config$get_response())
  missingPrec <- dplyr::rename(missingPrec, !!sym(meanAbundance) := .data$meanAbundance )

  p <- ggplot2::ggplot(missingPrec, ggplot2::aes(x = !!sym(meanAbundance), fill = .data$nrNAs, colour = .data$nrNAs)) +
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
missingness_per_condition_cumsum <- function(x,
                                             config,
                                             factors = config$factor_keys_depth()){
  missingPrec <- interaction_missing_stats(x, config,factors)$data

  xx <- missingPrec |> group_by(across(all_of(c(config$isotopeLabel, factors,"nrNAs","nrReplicates")))) |>
    dplyr::summarize(nrTransitions = n())

  xxcs <- xx |> group_by(across(all_of(c(config$isotopeLabel,factors)))) |> arrange(.data$nrNAs) |>
    dplyr::mutate(cumulative_sum = cumsum(.data$nrTransitions))
  res <- xxcs  |> dplyr::select(-nrTransitions)

  formula <- paste(config$isotopeLabel, "~", paste(factors, collapse = "+"))
  message(formula)

  nudgeval = mean(res$cumulative_sum) * 0.05
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
missingness_per_condition <- function(x, config, factors = config$factor_keys_depth()){
  missingPrec <- interaction_missing_stats(x, config, factors)$data
  hierarchyKey <- tail(config$hierarchy_keys(),1)
  hierarchyKey <- paste0("nr_",hierarchyKey)
  xx <- missingPrec |> group_by(across(all_of(c(config$isotopeLabel,
                                     factors,"nrNAs","nrReplicates")))) |>
    dplyr::summarize( !!sym(hierarchyKey) := n())

  formula <- paste(config$isotopeLabel, "~", paste(factors, collapse = "+"))
  #message(formula)

  nudgeval = max(xx[[hierarchyKey]]) * 0.05

  p <- ggplot(xx, aes(x = .data$nrNAs, y = .data[[hierarchyKey]])) +
    geom_bar(stat = "identity", color = "black", fill = "white") +
    geom_text(aes(label = !!sym(hierarchyKey)), nudge_y = nudgeval, angle = 45) +
    facet_grid(as.formula(formula))
  xx <- tidyr::pivot_wider(xx, names_from = "nrNAs", values_from = hierarchyKey)

  return(list(data = xx ,figure = p))
}


#' UpSetR plot from interaction_missing_stats
#'
#' @export
#' @keywords internal
#' @family plotting
#' @family imputation
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
  tmp <- prolfqua::interaction_missing_stats(data, cf)
  nrMiss <- tmp$data |> tidyr::pivot_wider(id_cols = cf$hierarchy_keys(),
                                           names_from = cf$factor_keys_depth(),
                                           values_from = !!rlang::sym("nrMeasured"))

  hl <- length(cf$hierarchy_keys())
  nrMiss[,-(1:hl)][nrMiss[,-(1:hl)] < tr] <- 0
  nrMiss[,-(1:hl)][nrMiss[,-(1:hl)] >= tr] <- 1
  return(list(data = as.data.frame(nrMiss), nsets = ncol(nrMiss) - length(cf$hierarchy_keys())))
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
UpSet_missing_stats <- function(data, config){
  data <- prolfqua::complete_cases(data, config)
  responseName <- config$get_response()
  data <- data |> dplyr::mutate(isThere =
                                  dplyr::case_when(
                                    !is.na(!!rlang::sym(responseName)) ~ 1,
                                    TRUE ~ 0
                                  )
  )
  pups2 <- data |> tidyr::pivot_wider(id_cols = config$hierarchy_keys(),
                                      names_from = config$sampleName,
                                      values_from = !!rlang::sym("isThere"))
  #colnames(pups2) <- make.names(colnames(pups2))
  res <- list(data = as.data.frame(pups2), nsets = ncol(pups2) - length(config$hierarchy_keys()) )
  return(res)
}

