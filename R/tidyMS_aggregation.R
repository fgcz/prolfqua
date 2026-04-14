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
  separate = FALSE,
  log_y = FALSE,
  show.legend = FALSE
) {
  if (length(isotope_label)) {
    if (separate) {
      formula <- paste(paste(isotope_label, collapse = "+"), "~", paste(factor, collapse = "+"))
      p <- ggplot(
        data,
        aes(
          x = .data[[sample]],
          y = .data[[intensity]],
          group = .data[[fragment]],
          color = .data[[peptide]]
        )
      )
    } else {
      formula <- sprintf("~%s", paste(factor, collapse = " + "))
      data <- tidyr::unite(data, "fragment_label", dplyr::all_of(c(fragment, isotope_label)), remove = FALSE)
      p <- ggplot(
        data,
        aes(
          x = .data[[sample]],
          y = .data[[intensity]],
          group = .data[["fragment_label"]],
          color = .data[[peptide]]
        )
      )
    }
    p <- p +
      geom_point(aes(shape = .data[[isotope_label]]), show.legend = show.legend) +
      geom_line(aes(linetype = .data[[isotope_label]]), show.legend = show.legend)
  } else {
    formula <- sprintf("~%s", paste(factor, collapse = " + "))
    p <- ggplot(
      data,
      aes(x = .data[[sample]], y = .data[[intensity]], group = .data[[fragment]], color = .data[[peptide]])
    )
    p <- p + geom_point(show.legend = show.legend) + geom_line(show.legend = show.legend)
  }

  # p <- ggplot(
  #   data,
  #   aes(x = .data[[sample]], y = .data[[intensity]], group = .data[[fragment]],
  #       color = .data[[peptide]], linetype = .data[[isotope_label]])
  # )
  p <- p + facet_grid(as.formula(formula), scales = "free_x")
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
#' @param separate if heavy and light show in one plot or with separate y axis?
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
plot_hierarchies_line <- function(res, protein_name, lfqdata, separate = FALSE, show.legend = FALSE) {
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
    separate = separate,
    log_y = !lfqdata$is_transformed(),
    show.legend = show.legend
  )
  return(res)
}


#' Generates peptide level plots for all proteins
#'
#' @seealso \code{\link{plot_hierarchies_line}}
#' @export
#' @param pdata data.frame
#' @param lfqdata LFQData object
#' @family aggregation
#' @family plotting
#' @keywords internal
#' @examples
#'
#'
#' istar <- sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' res <- plot_hierarchies_line_df(lfq$data_long(), lfq)
#' res[[1]]
#'
plot_hierarchies_line_df <- function(pdata, lfqdata, show.legend = FALSE) {
  hierarchy_id <- "hierarchy_id"
  pdata <- pdata |> tidyr::unite(hierarchy_id, !!!syms(lfqdata$relevant_hierarchy_keys()), remove = FALSE)

  xnested <- pdata |>
    dplyr::group_by(across(all_of(hierarchy_id))) |>
    tidyr::nest()

  figs <- xnested |>
    dplyr::mutate(
      plot = map2(data, !!sym(hierarchy_id), plot_hierarchies_line, lfqdata = lfqdata, show.legend = show.legend)
    )
  return(figs$plot)
}


#' Add protein estimate to plot of peptide intensities
#' @export
#' @family aggregation
#' @family plotting
#' @keywords internal
#' @examples
#' # todo
plot_hierarchies_add_quantline <- function(p, data, aes_y, sample_name) {
  p +
    geom_line(
      data = data,
      aes(x = .data[[sample_name]], y = .data[[aes_y]], group = 1),
      linewidth = 1.3,
      color = "black",
      linetype = "solid"
    ) +
    geom_point(
      data = data,
      aes(x = .data[[sample_name]], y = .data[[aes_y]], group = 1),
      color = "black",
      shape = 10
    )
}


.reestablish_condition <- function(data, medpolishRes, sample_name, factor_keys, file_name, isotope_label) {
  xx <- data |>
    dplyr::select(c(sample_name, factor_keys, file_name, isotope_label)) |>
    dplyr::distinct()
  res <- dplyr::inner_join(xx, medpolishRes, by = sample_name)
  res
}


#' Median polish estimate of e.g. protein from peptide intensities
#'
#' Compute Tukey's median polish estimate of a protein from peptide or precursor intensities
#'
#' @family aggregation
#' @param x a matrix
#' @param name if TRUE returns the name of the summary column
#' @return data.frame with number of rows equal to number of columns of input matrix.
#' @export
#' @keywords internal
#'
#' @examples
#'
#' medpolish_estimate(name = TRUE)
#' gg <- matrix(runif(20), 4, 5)
#' rownames(gg) <- paste0("A", 1:4)
#' colnames(gg) <- make.names(1:5)
#' gg
#' mx <- medpolish_estimate(gg)
#'
medpolish_estimate <- function(x, name = FALSE, sample_name = "sample_name") {
  if (name) {
    return("medpolish")
  }
  X <- medpolish(x, na.rm = TRUE, trace.iter = FALSE, maxiter = 10)
  res <- tibble(!!sample_name := names(X$col), medpolish = X$col + X$overall)
  res
}

.extractInt <- function(pdata, response, feature, sample_name) {
  pdata <- pdata |>
    dplyr::select(all_of(c(
      sample_name,
      feature,
      response
    ))) |>
    tidyr::pivot_wider(names_from = all_of(sample_name), values_from = all_of(response)) |>
    .ExtractMatrix()
  return(pdata)
}

#' Extract response column of a protein into matrix
#' Median polish estimates of e.g. protein abundances for entire data.frame
#'
#' @seealso \code{\link{medpolish_estimate}}
#' @param pdata data.frame
#' @param response column name with intensities
#' @param feature column name e.g. peptide ids
#' @param sample_name column name e.g. sample_name
#' @return data.frame
#' @export
#' @keywords internal
#' @family aggregation
#' @family plotting
#' @examples
#'
#' bb <- sim_lfq_data_peptide_config(Nprot = 20)
#' conf <- bb$config
#' data <- bb$data
#'
#' conf$hierarchy_depth <- 1
#' xnested <- data |>
#'   dplyr::group_by(across(all_of(conf$hierarchy_keys_depth()))) |>
#'   tidyr::nest()
#'
#' feature <- base::setdiff(
#'   conf$hierarchy_keys(),
#'   conf$hierarchy_keys_depth()
#' )
#' x <- xnested$data[[1]]
#' bb <- medpolish_estimate_df(x,
#'   response = conf$get_response(),
#'   feature = feature,
#'   sample_name = conf$sample_name
#' )
#' prolfqua:::.reestablish_condition(x, bb, conf$sample_name,
#'   conf$factor_keys(), conf$file_name, conf$isotope_label)
#'
medpolish_estimate_df <- function(pdata, response, feature, sample_name) {
  bb <- .extractInt(pdata, response = response, feature = feature, sample_name = sample_name)
  medpolish_estimate(bb, sample_name = sample_name)
}


#' Median polish estimates of e.g. protein abundances for entire data.frame
#'
#' @seealso \code{\link{medpolish_estimate_df}}
#' @param pdata data.frame
#' @param response character — intensity column name
#' @param hierarchy_keys character vector — all hierarchy column names
#' @param hierarchy_keys_depth character vector — hierarchy columns at current depth
#' @param sample_name character — sample name column
#' @param name if TRUE return function name only
#' @family aggregation
#' @keywords internal
#' @export
#' @examples
#'
#' bb <- sim_lfq_data_peptide_config(Nprot = 20)
#' conf <- bb$config
#' data <- bb$data
#' conf$hierarchy_depth <- 1
#' xnested <- data |>
#'   dplyr::group_by(across(all_of(conf$hierarchy_keys_depth()))) |>
#'   tidyr::nest()
#'
#' feature <- base::setdiff(conf$hierarchy_keys(), conf$hierarchy_keys_depth())
#' x <- xnested$data[[1]]
#' bb <- medpolish_estimate_dfconfig(x, conf$get_response(),
#'   conf$hierarchy_keys(), conf$hierarchy_keys_depth(), conf$sample_name)
#' prolfqua:::.reestablish_condition(x, bb, conf$sample_name,
#'   conf$factor_keys(), conf$file_name, conf$isotope_label)
#'
medpolish_estimate_dfconfig <- function(
  pdata,
  response,
  hierarchy_keys,
  hierarchy_keys_depth,
  sample_name,
  name = FALSE
) {
  if (name) {
    return("medpolish")
  }

  feature <- base::setdiff(hierarchy_keys, hierarchy_keys_depth)
  res <- medpolish_estimate_df(
    pdata,
    response = response,
    feature = feature,
    sample_name = sample_name
  )
  return(res)
}


.rlm_estimate <- function(pdata, response, feature, samples = "samples", maxIt = 20) {
  data <- pdata |>
    select(all_of(c(samples, feature, response))) |>
    na.omit()
  ## If there is only one 1 peptide for all samples return response of that peptide
  expname <- paste0("mean.", response)

  if (length(unique(data[[feature]])) == 1L) {
    data$lmrob <- data[[response]]
    data$weights <- 1
    data <- rename(data, !!expname := !!sym(response))
    data <- data |> select(-!!sym(feature))
    return(data)
  }

  ## model-matrix breaks on factors with 1 level so make vector of ones (will be intercept)
  if (length(unique(data[[samples]])) == 1L) {
    data <- data |>
      group_by(across(all_of(samples))) |>
      summarize(
        lmrob = mean(!!sym(response)),
        !!expname := mean(!!sym(response)),
        .groups = "drop"
      )
    data$weights <- 1
    return(data)
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
    sumdata <- data |>
      select(-!!sym(feature)) |>
      group_by(across(all_of(samples))) |>
      dplyr::summarize(!!expname := mean(!!sym(response)), lmrob = mean(!!sym(response)), .groups = "drop")
    sumdata$weights <- 1
    res <- sumdata
  }
  pdata <- pdata |>
    dplyr::select(all_of(samples)) |>
    distinct()
  res <- left_join(pdata, res, by = samples)
  return(res)
}


#' Estimate e.g. protein abundance from peptides using MASS:rlm
#'
#' @keywords internal
#' @family aggregation
#' @param pdata data
#' @param response intensities
#' @param feature e.g. peptideIDs.
#' @param samples e.g. sample_name
#' @importFrom MASS rlm
#' @export
#' @examples
#'
#' xx <- data.frame(response = rnorm(20, 0, 10), feature = rep(LETTERS[1:5], 4),
#'   samples = rep(letters[1:4], 5))
#'
#' bb <- rlm_estimate(xx, "response", "feature", "samples", maxIt = 20)
#'
#' xx2 <- data.frame(log2Area = rnorm(20, 0, 10), peptide_Id = rep(LETTERS[1:5], 4),
#'   sample_name = rep(letters[1:4], 5))
#' rlm_estimate(xx2, "log2Area", "peptide_Id", "sample_name")
#' rlm_estimate(prolfqua_data("data_checksummarizationrobust87"),
#'   "log2Area", "peptide_Id", "sampleName")
#' rlm_estimate(prolfqua_data("data_checksummarizerobust69"),
#'   "log2Area", "peptide_Id", "sampleName")
#' res <- vector(100, mode = "list")
#' for (i in seq_len(100)) {
#'   xx3 <- xx2
#'   xx3$log2Area[sample(1:20, sample(1:15, 1))] <- NA
#'   res[[i]] <- list(data = xx3, summary = rlm_estimate(xx3, "log2Area", "peptide_Id", "sample_name"))
#' }
#' rlm_estimate(xx2[xx2$peptide_Id == "A", ], "log2Area", "peptide_Id", "sample_name")
#' rlm_estimate(xx2[xx2$sample_name == "a", ], "log2Area", "peptide_Id", "sample_name")
#'
#'
#' bb <- sim_lfq_data_peptide_config(Nprot = 20)
#' conf <- bb$config
#' data <- bb$data
#' conf$hierarchy_depth <- 1
#' xnested <- data |>
#'   dplyr::group_by(across(all_of(conf$hierarchy_keys_depth()))) |>
#'   tidyr::nest()
#'
#' feature <- base::setdiff(conf$hierarchy_keys(), conf$hierarchy_keys_depth())
#' x <- xnested$data[[1]]
#' bb <- rlm_estimate(x,
#'   response = conf$get_response(),
#'   feature = feature,
#'   samples = conf$sample_name
#' )
#'
rlm_estimate <- function(pdata, response, feature, samples, maxIt = 20) {
  pdata <- unite(pdata, "feature", all_of(feature))

  res <- .rlm_estimate(pdata, response, feature = "feature", samples, maxIt = maxIt)
  return(res)
}

#' Estimate protein abundance from peptide abundances using MASS::rlm
#'
#'
#'
#' @seealso \code{\link{rlm_estimate}}
#' @export
#' @param pdata data.frame
#' @param response character — intensity column name
#' @param hierarchy_keys character vector — all hierarchy column names
#' @param hierarchy_keys_depth character vector — hierarchy columns at current depth
#' @param sample_name character — sample name column
#' @param name if TRUE return function name only
#' @family aggregation
#' @keywords internal
#' @examples
#'
#' bb <- sim_lfq_data_peptide_config(Nprot = 20)
#' conf <- bb$config
#' data <- bb$data
#' conf$hierarchy_depth <- 1
#' xnested <- data |>
#'   dplyr::group_by(across(all_of(conf$hierarchy_keys_depth()))) |>
#'   tidyr::nest()
#'
#' feature <- base::setdiff(conf$hierarchy_keys(), conf$hierarchy_keys_depth())
#' x <- xnested$data[[1]]
#' bb <- rlm_estimate_dfconfig(x, conf$get_response(),
#'   conf$hierarchy_keys(), conf$hierarchy_keys_depth(), conf$sample_name)
#' prolfqua:::.reestablish_condition(x, bb, conf$sample_name,
#'   conf$factor_keys(), conf$file_name, conf$isotope_label)
#'
rlm_estimate_dfconfig <- function(pdata, response, hierarchy_keys, hierarchy_keys_depth, sample_name, name = FALSE) {
  if (name) {
    return("lmrob")
  }

  feature <- base::setdiff(hierarchy_keys, hierarchy_keys_depth)
  rlm_estimate(pdata, response = response, feature = feature, samples = sample_name, maxIt = 20)
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
#' @seealso \code{\link{medpolish_estimate_dfconfig}} \code{\link{rlm_estimate_dfconfig}}
#' @param func - a function working on a matrix of intensities for each protein.
#' @return returns list with data (data.frame) and config (AnalysisConfiguration)
#' @family aggregation
#' @keywords internal
#' @export
#' @examples
#'
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(dd$data, dd$config)
#' lfq <- lfq$get_Transformer()$log2()$lfq
#' bbMed <- estimate_intensity(lfq, .func = medpolish_estimate_dfconfig)
#' bbRob <- estimate_intensity(lfq, .func = rlm_estimate_dfconfig)
#' nrow(bbMed$data)
#' nrow(bbRob$data)
#' xt <- dplyr::inner_join(bbMed$data, bbRob$data)
#' plot(xt$medpolish, xt$lmrob, log = "xy", pch = "*")
#' abline(0, 1, col = 2)
#'
estimate_intensity <- function(lfqdata, .func) {
  make_name <- .func(name = TRUE)
  config <- lfqdata$config$clone(deep = TRUE)
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

  xnested <- data |>
    group_by(across(all_of(hkeysd))) |>
    nest()

  pb <- progress::progress_bar$new(total = nrow(xnested))
  message("starting aggregation")
  res <- purrr::map(xnested$data, function(d) {
    pb$tick()
    aggr <- .func(d, response, hkeys, hkeysd, sname)
    .reestablish_condition(d, aggr, sname, fkeys, fname, iso)
  })

  xnested[[make_name]] <- res
  newconfig <- make_reduced_hierarchy_config(
    config,
    work_intensity = .func(name = TRUE),
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
#' bbMed <- estimate_intensity(lfq, .func = medpolish_estimate_dfconfig)
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

  pb <- progress::progress_bar$new(total = nrow(xnested_all))
  for (i in seq_len(nrow(xnested_all))) {
    p1 <- plot_hierarchies_line(
      xnested_all$data[[i]],
      xnested_all[[hierarchy_id]][i],
      lfqdata = lfqdata,
      show.legend = show.legend
    )
    p2 <- plot_hierarchies_add_quantline(
      p1,
      xnested_all$other[[i]],
      lfqdata_agg$response(),
      lfqdata$sample_name()
    )
    plots[[i]] <- p2
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
  config <- lfqdata$config$clone(deep = TRUE)

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


.ExtractMatrix <- function(x, sep = "~lfq~") {
  idx <- sapply(x, is.numeric)
  xmat <- as.matrix(x[, idx])
  idcols <- x |> dplyr::select(which(!idx == TRUE))
  if (ncol(idcols) > 0) {
    rownames(xmat) <- x |> dplyr::select(which(!idx == TRUE)) |> tidyr::unite(x, sep = sep) |> dplyr::pull(x)
  }
  xmat
}


#' Aggregates e.g. protein abundances from peptide abundances
#'
#' @param data data.frame
#' @param response character — intensity column name
#' @param hierarchy_keys_depth character vector — hierarchy columns at current depth
#' @param file_name character — file name column
#' @param nr_children_col character — nr_children column name
#' @param new_child character — output column name
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
#' # xp |> pivot_wider(id_cols = protein_Id, names_from = sample, values_from = nr_peptides)
#' xp$nr_peptides |> table()
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


#' Count observations per experiment (max nr_children across samples)
#'
#' @export
#' @param data data.frame
#' @param hierarchy_keys character vector — all hierarchy column names
#' @param hierarchy_keys_depth character vector — hierarchy columns at current depth
#' @param nr_children_col character — nr_children column name
#' @param response character — intensity column name (needed when from_children = TRUE)
#' @param file_name character — file name column (needed when from_children = TRUE)
#' @param from_children TRUE compute from nr_children column, FALSE count distinct hierarchy entries
#' @param name_nr_child how to name output column
#' @examples
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(dd$data, dd$config)
#' xd <- nr_obs_experiment(lfq$data_long(), lfq$hierarchy_keys(),
#'   lfq$relevant_hierarchy_keys(), lfq$nr_children_col(),
#'   response = lfq$response(), file_name = lfq$file_name())
#' stopifnot(min(xd$nr_child_exp) == 1)
#' xd2 <- nr_obs_experiment(lfq$data_long(), lfq$hierarchy_keys(),
#'   lfq$relevant_hierarchy_keys(), lfq$nr_children_col(), from_children = FALSE)
#'
nr_obs_experiment <- function(
  data,
  hierarchy_keys,
  hierarchy_keys_depth,
  nr_children_col,
  response = NULL,
  file_name = NULL,
  from_children = TRUE,
  name_nr_child = "nr_child_exp"
) {
  if (from_children) {
    stopifnot(!is.null(response), !is.null(file_name))
    xz <- nr_obs_sample(data, response, hierarchy_keys_depth, file_name, nr_children_col)
    xz <- xz |>
      group_by(!!!syms(hierarchy_keys_depth)) |>
      dplyr::summarize(!!name_nr_child := max(!!sym(nr_children_col)), .groups = "drop")
    return(xz)
  } else {
    xq <- data |>
      dplyr::select(hierarchy_keys) |>
      distinct() |>
      dplyr::group_by(!!!syms(hierarchy_keys_depth)) |>
      dplyr::summarize(!!name_nr_child := dplyr::n(), .groups = "drop")
    return(xq)
  }
}


# Summarize Intensities by Intensity or NAs ----
.rankProteinPrecursors <- function(
  data,
  hierarchy_keys,
  column,
  fun = function(x) {
    sum(x, na.rm = TRUE)
  },
  summary_column = "srm_meanInt",
  rank_column = "srm_meanIntRank",
  rank_function = function(x) {
    min_rank(desc(x))
  }
) {
  summary_per_precursor <- data |>
    dplyr::group_by(!!!syms(hierarchy_keys)) |>
    dplyr::summarize(!!summary_column := fun(!!sym(column)))

  grouped_by_protein <- summary_per_precursor |>
    dplyr::arrange(!!sym(hierarchy_keys[1])) |>
    dplyr::group_by(!!sym(hierarchy_keys[1]))
  ranked_by_summary <- grouped_by_protein |>
    dplyr::mutate(!!rank_column := rank_function(!!sym(summary_column)))

  data <- dplyr::inner_join(data, ranked_by_summary)
  return(data)
}

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
  pdata <- .rankProteinPrecursors(
    pdata,
    hierarchy_keys,
    column = response,
    fun = function(x) {
      mean(x, na.rm = TRUE)
    },
    summary_column = summary_column,
    rank_column = rank_column,
    rank_function = function(x) {
      min_rank(desc(x))
    }
  )

  message("Columns added : ", summary_column, " ", rank_column)
  return(pdata)
}
