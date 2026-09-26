# Functions - Plotting ----
# Plot peptide and fragments
plot_hierarchies_line_default <- function(
  data,
  protein_name,
  sample,
  intensity,
  peptide,
  fragment,
  factor,
  isotope_label,
  log_y = FALSE,
  show.legend = FALSE
) {
  if (length(isotope_label)) {
    data <- tidyr::unite(data, "fragment_label", dplyr::all_of(c(fragment, isotope_label)), remove = FALSE)
    fragment <- "fragment_label"
    points <- geom_point(aes(shape = .data[[isotope_label]]), show.legend = show.legend)
    lines <- geom_line(aes(linetype = .data[[isotope_label]]), show.legend = show.legend)
  } else {
    points <- geom_point(show.legend = show.legend)
    lines <- geom_line(show.legend = show.legend)
  }
  p <- ggplot(
    data,
    aes(x = .data[[sample]], y = .data[[intensity]], group = .data[[fragment]], color = .data[[peptide]])
  )
  p <- p + points + lines
  p <- p + facet_grid(as.formula(sprintf("~%s", paste(factor, collapse = " + "))), scales = "free_x")
  p <- p + ggtitle(protein_name)
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top")
  if (log_y) {
    p <- p + scale_y_continuous(trans = "log10")
  }
  return(p)
}

#' Plot peptide intensities of protein as a function of the sample and factor
#'
#' @export
#' @param res data.frame
#' @param protein_name title of plot
#' @param lfqdata LFQData object
#' @param show.legend logical, show legend in plot
#' @family aggregation
#' @family plotting
#'
#' @keywords internal
#' @examples
#'
#'
#' istar <- sim_lfq_data_peptide_config()
#' config <- istar$config
#' analysis <- istar$data
#'
#' xnested <- analysis |>
#'   dplyr::group_by(across(all_of(config$hierarchy_keys_depth()))) |>
#'   tidyr::nest()
#'
#' lfq <- LFQData$new(analysis, config)
#' prolfqua::plot_hierarchies_line(xnested$data[[1]], xnested$protein_Id[[1]], lfq)
#'
plot_hierarchies_line <- function(res, protein_name, lfqdata, show.legend = FALSE) {
  rev_hnames <- lfqdata$hierarchy_keys()
  rev_hnames <- rev(rev_hnames)
  fragment <- rev_hnames[1]
  peptide <- rev_hnames[1]

  if (length(rev_hnames) > 2) {
    peptide <- rev_hnames[2]
  }
  res <- plot_hierarchies_line_default(
    res,
    protein_name = protein_name,
    sample = lfqdata$sample_name(),
    intensity = lfqdata$response(),
    peptide = peptide,
    fragment = fragment,
    factor = lfqdata$relevant_factor_keys(),
    isotope_label = lfqdata$isotope_label(),
    log_y = !lfqdata$is_transformed(),
    show.legend = show.legend
  )
  return(res)
}


.reestablish_condition <- function(data, medpolishRes, sample_name, factor_keys, file_name, isotope_label) {
  xx <- data |>
    dplyr::select(c(sample_name, factor_keys, file_name, isotope_label)) |>
    dplyr::distinct()
  res <- dplyr::inner_join(xx, medpolishRes, by = sample_name)
  res
}


# Tukey's median polish estimate of e.g. a protein from its feature (peptide) intensities, one row per sample.
.medpolish_estimate <- function(pdata, response, feature, sample_name) {
  wide <- pdata |>
    dplyr::select(all_of(c(sample_name, feature, response))) |>
    tidyr::pivot_wider(names_from = all_of(sample_name), values_from = all_of(response))
  mat <- as.matrix(wide[, vapply(wide, is.numeric, logical(1))])
  X <- medpolish(mat, na.rm = TRUE, trace.iter = FALSE, maxiter = 10)
  tibble(!!sample_name := names(X$col), medpolish = X$col + X$overall)
}

.rlm_estimate <- function(pdata, response, feature, samples = "samples", maxIt = 20) {
  data <- pdata |>
    select(all_of(c(samples, feature, response))) |>
    na.omit()
  expname <- paste0("mean.", response)

  ## All branches hand off here so the output schema is identical: a single, fixed column order
  ## (samples, mean.<response>, lmrob, weights) with one row per sample present in the input.
  ## The final left_join re-establishes samples dropped by na.omit, so downstream
  ## unnest()/.reestablish_condition() always bind the same columns and never lose a sample.
  all_samples <- pdata |>
    dplyr::select(all_of(samples)) |>
    distinct()
  finalize <- function(per_sample) {
    per_sample <- per_sample |>
      dplyr::select(all_of(c(samples, expname, "lmrob", "weights")))
    dplyr::left_join(all_samples, per_sample, by = samples)
  }

  ## If there is only 1 peptide for all samples use the response of that peptide directly
  if (length(unique(data[[feature]])) == 1L) {
    per_sample <- data |>
      dplyr::mutate(lmrob = !!sym(response), weights = 1) |>
      rename(!!expname := !!sym(response))
    return(finalize(per_sample))
  }

  ## model-matrix breaks on factors with 1 level so make vector of ones (will be intercept)
  if (length(unique(data[[samples]])) == 1L) {
    per_sample <- data |>
      group_by(across(all_of(samples))) |>
      summarize(
        lmrob = mean(!!sym(response)),
        !!expname := mean(!!sym(response)),
        weights = 1,
        .groups = "drop"
      )
    return(finalize(per_sample))
  }

  ## sum contrast on peptide level so sample effect will be mean over all peptides instead of reference level
  formula <- as.formula(paste0("~ -1 + ", samples, " + ", feature))
  contr.arg <- list("contr.sum")
  names(contr.arg) <- feature
  X <- model.matrix(formula, data = data, contrasts.arg = contr.arg)
  ## MASS::rlm breaks on singular values.
  ## check with base lm if singular values are present.
  ## if so, these coefficients will be zero, remove this column from model matrix
  ## rinse and repeat on reduced model-matrix untill no singular values are present
  y <- data[[response]]
  repeat {
    fit <- .lm.fit(X, y)
    id <- fit$coefficients != 0
    X <- X[, id, drop = FALSE]
    if (!any(!id)) break
  }
  ## Last step is always rlm if X > has some columns left
  if (ncol(X) > 0) {
    fit <- MASS::rlm(X, y, maxit = maxIt)
    data$residuals <- fit$residuals
    usamples <- unique(data[[samples]])
    coef_names <- paste0(samples, usamples)
    lmrob <- tibble(!!samples := usamples, lmrob = fit$coefficients[coef_names])

    sumdata <- data |>
      select(-!!sym(feature)) |>
      group_by(across(all_of(samples))) |>
      dplyr::summarize(!!expname := mean(!!sym(response)), weights = 1 / mean(residuals^2), .groups = "drop")
    if (any(is.infinite(sumdata$weights) | is.na(sumdata$weights) | sumdata$weights > 10e6)) {
      sumdata$weights <- 1
    }

    res <- inner_join(sumdata, lmrob, by = samples)
    if (sum(is.na(res[[expname]])) < sum(is.na(res$lmrob))) {
      res$lmrob <- res[[expname]]
    }
  } else {
    res <- data |>
      select(-!!sym(feature)) |>
      group_by(across(all_of(samples))) |>
      dplyr::summarize(!!expname := mean(!!sym(response)), lmrob = mean(!!sym(response)), .groups = "drop")
    res$weights <- 1
  }
  return(finalize(res))
}

.add_nr_children <- function(
  data,
  aggregated_data,
  hierarchy_keys_depth,
  file_name,
  response,
  nr_children_col,
  newconfig
) {
  new_child <- paste(c("nr_children", hierarchy_keys_depth), collapse = "_")
  res_nr_children <- nr_obs_sample(
    data,
    response,
    hierarchy_keys_depth,
    file_name,
    nr_children_col,
    new_child = new_child
  )
  result <- inner_join(aggregated_data, res_nr_children, by = c(hierarchy_keys_depth, file_name))
  newconfig$nr_children <- new_child
  return(list(data = result, config = newconfig))
}

#' Aggregates e.g. protein abundances from peptide abundances
#'
#' @param lfqdata LFQData object
#' @param method "medpolish" (Tukey's median polish, column `medpolish`) or "rlm" (robust regression with
#'   `MASS::rlm`, column `lmrob`)
#' @return returns list with data (data.frame) and config (AnalysisConfiguration)
#' @family aggregation
#' @keywords internal
#' @importFrom MASS rlm
#' @export
#' @examples
#'
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(dd$data, dd$config)
#' lfq <- lfq$get_Transformer()$log2()$lfq
#' bbMed <- estimate_intensity(lfq, method = "medpolish")
#' bbRob <- estimate_intensity(lfq, method = "rlm")
#' nrow(bbMed$data)
#' nrow(bbRob$data)
#' xt <- dplyr::inner_join(bbMed$data, bbRob$data)
#' plot(xt$medpolish, xt$lmrob, log = "xy", pch = "*")
#' abline(0, 1, col = 2)
#'
estimate_intensity <- function(lfqdata, method = c("medpolish", "rlm")) {
  method <- match.arg(method)
  make_name <- if (method == "medpolish") "medpolish" else "lmrob"
  config <- lfqdata$get_config()$clone(deep = TRUE)
  data <- lfqdata$data_long()

  # Extract column names once for passing to sub-functions
  response <- lfqdata$response()
  hkeys <- lfqdata$hierarchy_keys()
  hkeysd <- lfqdata$relevant_hierarchy_keys()
  sname <- lfqdata$sample_name()
  fname <- lfqdata$file_name()
  fkeys <- lfqdata$factor_keys()
  iso <- lfqdata$isotope_label()
  nrc <- lfqdata$nr_children_col()
  feature <- base::setdiff(hkeys, hkeysd)
  estimate <- if (method == "medpolish") {
    function(d) .medpolish_estimate(d, response, feature, sname)
  } else {
    function(d) .rlm_estimate(unite(d, "feature", all_of(feature)), response, "feature", sname)
  }

  xnested <- data |>
    group_by(across(all_of(hkeysd))) |>
    nest()

  pb <- progress::progress_bar$new(total = nrow(xnested))
  message("starting aggregation")
  res <- purrr::map(xnested$data, function(d) {
    pb$tick()
    .reestablish_condition(d, estimate(d), sname, fkeys, fname, iso)
  })

  xnested[[make_name]] <- res
  newconfig <- make_reduced_hierarchy_config(
    config,
    work_intensity = make_name,
    hierarchy = config$hierarchy_keys_depth(names = FALSE)
  )

  unnested <- xnested |>
    dplyr::select(all_of(c(hkeysd, make_name))) |>
    tidyr::unnest(cols = dplyr::all_of(make_name)) |>
    dplyr::ungroup()

  return(.add_nr_children(data, unnested, hkeysd, fname, response, nrc, newconfig))
}


#' Plot feature data and result of aggregation
#'
#' @param lfqdata LFQData object (original data)
#' @param lfqdata_agg LFQData object (aggregated data)
#' @param show.legend logical, show legend in plot
#' @family plotting
#' @family aggregation
#' @keywords internal
#' @export
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' lfq <- lfq$get_Transformer()$log2()$lfq
#' bbMed <- estimate_intensity(lfq, method = "medpolish")
#' lfq_med <- LFQData$new(bbMed$data, bbMed$config)
#' tmpMed <- plot_estimate(lfq, lfq_med)
#' stopifnot("ggplot" %in% class(tmpMed$plots[[1]]))
#'
plot_estimate <- function(lfqdata, lfqdata_agg, show.legend = FALSE) {
  hierarchy_id <- "hierarchy_id"
  hkeysd <- lfqdata$relevant_hierarchy_keys()

  xnested <- lfqdata$data_long() |>
    group_by(!!!syms(hkeysd)) |>
    nest()
  xnested <- xnested |> tidyr::unite(hierarchy_id, !!!syms(hkeysd))
  xnested_aggr <- lfqdata_agg$data_long() |>
    group_by(!!!syms(lfqdata_agg$relevant_hierarchy_keys())) |>
    nest_by(.key = "other")
  xnested_aggr <- xnested_aggr |> tidyr::unite(hierarchy_id, !!!syms(hkeysd))
  xnested_all <- inner_join(xnested, xnested_aggr, by = hierarchy_id)

  plots <- vector(mode = "list", length = nrow(xnested_all))
  # protein estimate drawn on top of the peptide intensities
  sample_name <- lfqdata$sample_name()
  aes_y <- lfqdata_agg$response()
  quant_aes <- aes(x = .data[[sample_name]], y = .data[[aes_y]], group = 1)

  pb <- progress::progress_bar$new(total = nrow(xnested_all))
  for (i in seq_len(nrow(xnested_all))) {
    p1 <- plot_hierarchies_line(
      xnested_all$data[[i]],
      xnested_all[[hierarchy_id]][i],
      lfqdata = lfqdata,
      show.legend = show.legend
    )
    plots[[i]] <- p1 +
      geom_line(data = xnested_all$other[[i]], quant_aes, linewidth = 1.3, color = "black", linetype = "solid") +
      geom_point(data = xnested_all$other[[i]], quant_aes, color = "black", shape = 10)
    pb$tick()
  }
  xnested_all$plots <- plots

  return(xnested_all)
}


#' Aggregates top N intensities
#'
#' run \link{rank_peptide_by_intensity} first
#' @param ranked_data data.frame with ranked peptides
#' @param lfqdata LFQData object
#' @param .func function to use for aggregation
#' @param N default 3 top intensities.
#' @return list with data and new reduced configuration (config)
#' @family aggregation
#' @export
#' @keywords internal
#' @examples
#'
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(dd$data, dd$config)
#' ranked <- rank_peptide_by_intensity(lfq$data_long(), lfq$response(), lfq$hierarchy_keys())
#'
#' mean_f <- function(x, name = FALSE) {
#'   if (name) return("mean")
#'   mean(x, na.rm = TRUE)
#' }
#'
#' resTOPN <- aggregate_intensity_top_n(ranked, lfq, .func = mean_f, N = 3)
#' stopifnot(names(resTOPN) %in% c("data", "config"))
#' lfq_agg <- LFQData$new(resTOPN$data, resTOPN$config)
#' tmpRob <- plot_estimate(lfq, lfq_agg, show.legend = TRUE)
#' stopifnot("ggplot" %in% class(tmpRob$plots[[4]]))
#'
aggregate_intensity_top_n <- function(ranked_data, lfqdata, .func, N = 3) {
  newcol <- make.names(paste0("srm_", .func(name = TRUE), "_", N))
  config <- lfqdata$get_config()$clone(deep = TRUE)

  top_intensities <-
    ranked_data |> dplyr::filter(!!sym("srm_meanIntRank") <= N)

  top_intensities <- top_intensities |>
    dplyr::group_by(across(all_of(c(
      lfqdata$relevant_hierarchy_keys(),
      lfqdata$sample_name(),
      lfqdata$file_name(),
      lfqdata$isotope_label(),
      lfqdata$factor_keys()
    ))))
  sum_top_intensities <- top_intensities |>
    dplyr::summarize(
      !!newcol := .func(!!sym(lfqdata$response())),
      ident_qValue = min(!!sym(config$ident_q_value)),
      .groups = "drop"
    )

  newconfig <- make_reduced_hierarchy_config(
    config,
    work_intensity = newcol,
    hierarchy = config$hierarchy[seq_len(config$hierarchy_depth)]
  )

  return(.add_nr_children(
    ranked_data,
    sum_top_intensities,
    lfqdata$relevant_hierarchy_keys(),
    lfqdata$file_name(),
    lfqdata$response(),
    lfqdata$nr_children_col(),
    newconfig
  ))
}


#' Aggregates e.g. protein abundances from peptide abundances
#'
#' @param data data.frame
#' @param response character — intensity column name
#' @param hierarchy_keys_depth character vector — hierarchy columns at current depth
#' @param file_name character — file name column
#' @param nr_children_col character — nr_children column name
#' @param new_child character — output column name
#' @return The computed result.
#' @export
#' @examples
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' dd$data <- na.omit(dd$data)
#'
#' xd <- nr_obs_sample(na.omit(dd$data), dd$config$get_response(),
#'   dd$config$hierarchy_keys_depth(), dd$config$file_name, dd$config$nr_children)
#' xd$nr_children |> table()
#' # xd |> pivot_wider(id_cols = protein_Id, names_from = sample, values_from = nr_children)
#'
#' dp <- prolfqua::sim_lfq_data_protein_config()
#' xp <- nr_obs_sample(dp$data, dp$config$get_response(),
#'   dp$config$hierarchy_keys_depth(), dp$config$file_name, dp$config$nr_children)
#' # xp
#' # xp |> pivot_wider(id_cols = protein_Id, names_from = sample, values_from = nrPeptides)
#' xp$nrPeptides |> table()
#'
nr_obs_sample <- function(
  data,
  response,
  hierarchy_keys_depth,
  file_name,
  nr_children_col,
  new_child = nr_children_col
) {
  data <- data[!is.na(data[[response]]), ]
  nr_children <- data |>
    group_by(!!!rlang::syms(c(hierarchy_keys_depth, file_name))) |>
    summarize(!!new_child := sum(!!sym(nr_children_col), na.rm = TRUE), .groups = "drop")
  return(nr_children)
}


#' Max nr_children across samples per hierarchy unit
#'
#' Aggregates nr_children per sample (via nr_obs_sample), then takes the max per hierarchy unit.
#'
#' @param data data.frame
#' @param response character — intensity column name
#' @param hierarchy_keys_depth character vector — hierarchy columns at current depth
#' @param file_name character — file name column
#' @param nr_children_col character — nr_children column name
#' @param name_nr_child character — output column name
#' @return The computed result.
#' @export
#' @examples
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(dd$data, dd$config)
#' xd <- nr_children_experiment(lfq$data_long(), lfq$response(),
#'   lfq$relevant_hierarchy_keys(), lfq$file_name(), lfq$nr_children_col())
#' stopifnot(min(xd$nr_child_exp) == 1)
#'
nr_children_experiment <- function(
  data,
  response,
  hierarchy_keys_depth,
  file_name,
  nr_children_col,
  name_nr_child = "nr_child_exp"
) {
  per_sample_col <- "nr_children_per_sample"
  # children (e.g. peptides) observed per hierarchy unit in EACH sample
  per_sample <- nr_obs_sample(
    data,
    response,
    hierarchy_keys_depth,
    file_name,
    nr_children_col = nr_children_col,
    new_child = per_sample_col
  )
  # experiment-wide count = the largest per-sample count across all samples
  per_sample |>
    dplyr::group_by(!!!syms(hierarchy_keys_depth)) |>
    dplyr::summarize(!!name_nr_child := max(!!sym(per_sample_col)), .groups = "drop")
}

#' Count distinct child features per hierarchy unit
#'
#' Counts the number of distinct child-level entries (e.g. peptides per protein).
#'
#' @param data data.frame
#' @param hierarchy_keys character vector — all hierarchy column names
#' @param hierarchy_keys_depth character vector — hierarchy columns at current depth
#' @param name_nr_child character — output column name
#' @return The computed result.
#' @export
#' @examples
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(dd$data, dd$config)
#' xd <- nr_features_experiment(lfq$data_long(), lfq$hierarchy_keys(),
#'   lfq$relevant_hierarchy_keys())
#' stopifnot(min(xd$nr_child_exp) == 1)
#'
nr_features_experiment <- function(data, hierarchy_keys, hierarchy_keys_depth, name_nr_child = "nr_child_exp") {
  data |>
    dplyr::select(hierarchy_keys) |>
    distinct() |>
    dplyr::group_by(!!!syms(hierarchy_keys_depth)) |>
    dplyr::summarize(!!name_nr_child := dplyr::n(), .groups = "drop")
}

# Summarize Intensities by Intensity or NAs ----
#' ranks precursor - peptide by intensity.
#' @param pdata data.frame
#' @param response character — intensity column name
#' @param hierarchy_keys character vector — all hierarchy column names
#' @return data.frame
#' @export
#' @keywords internal
#' @examples
#'
#' bb <- prolfqua::sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(bb$data, bb$config)
#' res <- rank_peptide_by_intensity(lfq$data_long(), lfq$response(), lfq$hierarchy_keys())
#' X <- res |> dplyr::select(c(lfq$hierarchy_keys(),
#'  srm_meanInt, srm_meanIntRank)) |> dplyr::distinct()
#' X |> dplyr::arrange(!!!rlang::syms(c(lfq$hierarchy_keys()[1], "srm_meanIntRank")))
rank_peptide_by_intensity <- function(pdata, response, hierarchy_keys) {
  summary_column <- "srm_meanInt"
  rank_column <- "srm_meanIntRank"
  # mean intensity of each precursor, ranked within its protein (hierarchy_keys[1])
  ranked <- pdata |>
    dplyr::group_by(!!!syms(hierarchy_keys)) |>
    dplyr::summarize(!!summary_column := mean(!!sym(response), na.rm = TRUE)) |>
    dplyr::arrange(!!sym(hierarchy_keys[1])) |>
    dplyr::group_by(!!sym(hierarchy_keys[1])) |>
    dplyr::mutate(!!rank_column := min_rank(desc(!!sym(summary_column))))
  pdata <- dplyr::inner_join(pdata, ranked)

  message("Columns added : ", summary_column, " ", rank_column)
  return(pdata)
}
