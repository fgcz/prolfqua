# Functions - Plotting ----
# Plot peptide and fragments
plot_hierarchies_line_default <- function(
  data,
  proteinName,
  sample,
  intensity,
  peptide,
  fragment,
  factor,
  isotopeLabel,
  separate = FALSE,
  log_y = FALSE,
  show.legend = FALSE
) {
  if (length(isotopeLabel)) {
    if (separate) {
      formula <- paste(paste(isotopeLabel, collapse = "+"), "~", paste(factor, collapse = "+"))
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
      data <- tidyr::unite(data, "fragment_label", dplyr::all_of(c(fragment, isotopeLabel)), remove = FALSE)
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
      geom_point(aes(shape = .data[[isotopeLabel]]), show.legend = show.legend) +
      geom_line(aes(linetype = .data[[isotopeLabel]]), show.legend = show.legend)
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
  #       color = .data[[peptide]], linetype = .data[[isotopeLabel]])
  # )
  p <- p + facet_grid(as.formula(formula), scales = "free_x")
  p <- p + ggtitle(proteinName)
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
#' @param proteinName title of plot
#' @param config AnalysisConfiguration
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
#' prolfqua::plot_hierarchies_line(xnested$data[[1]], xnested$protein_Id[[1]], config)
#'
#' bb <- prolfqua_data("data_skylineSRM_HL_A")
#' conf <- bb$config_f()
#' analysis <- bb$analysis(bb$data, conf)
#'
#' nest <- analysis |>
#'   dplyr::group_by(conf$hierarchy_keys_depth()) |>
#'   tidyr::nest()
#' prolfqua::plot_hierarchies_line(nest$data[[1]],
#'   "DUM",
#'   conf,
#'   separate = TRUE
#' )
#' prolfqua::plot_hierarchies_line(nest$data[[1]],
#'   "DUM",
#'   conf,
#'   separate = TRUE,
#'   show.legend = TRUE
#' )
#'
plot_hierarchies_line <- function(res, proteinName, config, separate = FALSE, show.legend = FALSE) {
  rev_hnames <- config$hierarchy_keys(TRUE)
  fragment <- rev_hnames[1]
  peptide <- rev_hnames[1]

  if (length(rev_hnames) > 2) {
    peptide <- rev_hnames[2]
  }
  res <- plot_hierarchies_line_default(
    res,
    proteinName = proteinName,
    sample = config$sampleName,
    intensity = config$get_response(),
    peptide = peptide,
    fragment = fragment,
    factor = config$factor_keys_depth(),
    isotopeLabel = config$isotopeLabel,
    separate = separate,
    log_y = !config$is_response_transformed,
    show.legend = show.legend
  )
  return(res)
}


#' Generates peptide level plots for all proteins
#'
#' @seealso \code{\link{plot_hierarchies_line}}
#' @export
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @family aggregation
#' @family plotting
#' @keywords internal
#' @examples
#'
#'
#' istar <- sim_lfq_data_peptide_config()
#'
#' istar$config$is_response_transformed <- FALSE
#' res <- plot_hierarchies_line_df(istar$data, istar$config)
#' res[[1]]
#' istar$config$is_response_transformed <- TRUE
#' res <- plot_hierarchies_line_df(istar$data, istar$config)
#' res[[2]]
#'
#' # TODO make it work for other hiearachy levels.
#'
plot_hierarchies_line_df <- function(pdata, config, show.legend = FALSE) {
  hierarchy_ID <- "hierarchy_ID"
  pdata <- pdata |> tidyr::unite(hierarchy_ID, !!!syms(config$hierarchy_keys_depth()), remove = FALSE)

  xnested <- pdata |>
    dplyr::group_by(across(all_of(hierarchy_ID))) |>
    tidyr::nest()

  figs <- xnested |>
    dplyr::mutate(
      plot = map2(data, !!sym(hierarchy_ID), plot_hierarchies_line, config = config, show.legend = show.legend)
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
plot_hierarchies_add_quantline <- function(p, data, aes_y, configuration) {
  p +
    geom_line(
      data = data,
      aes(x = .data[[configuration$sampleName]], y = .data[[aes_y]], group = 1),
      linewidth = 1.3,
      color = "black",
      linetype = "solid"
    ) +
    geom_point(
      data = data,
      aes(x = .data[[configuration$sampleName]], y = .data[[aes_y]], group = 1),
      color = "black",
      shape = 10
    )
}


.reestablish_condition <- function(data, medpolishRes, config) {
  xx <- data |>
    dplyr::select(c(
      config$sampleName,
      config$factor_keys(),
      config$fileName,
      config$isotopeLabel
    )) |>
    dplyr::distinct()
  res <- dplyr::inner_join(xx, medpolishRes, by = config$sampleName)
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
medpolish_estimate <- function(x, name = FALSE, sampleName = "sampleName") {
  if (name) {
    return("medpolish")
  }
  X <- medpolish(x, na.rm = TRUE, trace.iter = FALSE, maxiter = 10)
  res <- tibble(!!sampleName := names(X$col), medpolish = X$col + X$overall)
  res
}

.extractInt <- function(pdata, response, feature, sampleName) {
  pdata <- pdata |>
    dplyr::select(all_of(c(
      sampleName,
      feature,
      response
    ))) |>
    tidyr::pivot_wider(names_from = all_of(sampleName), values_from = all_of(response)) |>
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
#' @param sampleName column name e.g. sampleName
#' @return data.frame
#' @export
#' @keywords internal
#' @family aggregation
#' @family plotting
#' @examples
#'
#' bb <- prolfqua_data("data_ionstar")$filtered()
#' bb$config <- old2new(bb$config)
#' stopifnot(nrow(bb$data) == 25780)
#' conf <- bb$config
#' data <- bb$data
#'
#' conf$hierarchyDepth <- 1
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
#'   sampleName = conf$sampleName
#' )
#' prolfqua:::.reestablish_condition(x, bb, conf)
#'
medpolish_estimate_df <- function(pdata, response, feature, sampleName) {
  bb <- .extractInt(pdata, response = response, feature = feature, sampleName = sampleName)
  medpolish_estimate(bb, sampleName = sampleName)
}


#' Median polish estimates of e.g. protein abundances for entire data.frame
#'
#' @seealso \code{\link{medpolish_estimate_df}}
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @family aggregation
#' @keywords internal
#' @export
#' @examples
#'
#' bb <- prolfqua_data("data_ionstar")$filtered()
#' bb$config <- old2new(bb$config)
#' stopifnot(nrow(bb$data) == 25780)
#' conf <- bb$config
#' data <- bb$data
#' conf$hierarchyDepth <- 1
#' xnested <- data |>
#'   dplyr::group_by(across(all_of(conf$hierarchy_keys_depth()))) |>
#'   tidyr::nest()
#'
#' feature <- base::setdiff(conf$hierarchy_keys(), conf$hierarchy_keys_depth())
#' x <- xnested$data[[1]]
#' bb <- medpolish_estimate_dfconfig(x, conf)
#' prolfqua:::.reestablish_condition(x, bb, conf)
#'
medpolish_estimate_dfconfig <- function(pdata, config, name = FALSE) {
  if (name) {
    return("medpolish")
  }

  feature <- base::setdiff(config$hierarchy_keys(), config$hierarchy_keys_depth())
  res <- medpolish_estimate_df(
    pdata,
    response = config$get_response(),
    feature = feature,
    sampleName = config$sampleName
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
    coefNames <- paste0(samples, usamples)
    lmrob <- tibble(!!samples := usamples, lmrob = fit$coefficients[coefNames])

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
#' @param samples e.g. sampleName
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
#'   sampleName = rep(letters[1:4], 5))
#' rlm_estimate(xx2, "log2Area", "peptide_Id", "sampleName")
#' rlm_estimate(prolfqua_data("data_checksummarizationrobust87"),
#'   "log2Area", "peptide_Id", "sampleName")
#' rlm_estimate(prolfqua_data("data_checksummarizerobust69"),
#'   "log2Area", "peptide_Id", "sampleName")
#' res <- vector(100, mode = "list")
#' for (i in seq_len(100)) {
#'   xx3 <- xx2
#'   xx3$log2Area[sample(1:20, sample(1:15, 1))] <- NA
#'   res[[i]] <- list(data = xx3, summary = rlm_estimate(xx3, "log2Area", "peptide_Id", "sampleName"))
#' }
#' rlm_estimate(xx2[xx2$peptide_Id == "A", ], "log2Area", "peptide_Id", "sampleName")
#' rlm_estimate(xx2[xx2$sampleName == "a", ], "log2Area", "peptide_Id", "sampleName")
#'
#'
#' bb <- prolfqua_data("data_ionstar")$filtered()
#' bb$config <- old2new(bb$config)
#' stopifnot(nrow(bb$data) == 25780)
#' conf <- bb$config
#' data <- bb$data
#' conf$hierarchyDepth <- 1
#' xnested <- data |>
#'   dplyr::group_by(across(all_of(conf$hierarchy_keys_depth()))) |>
#'   tidyr::nest()
#'
#' feature <- base::setdiff(conf$hierarchy_keys(), conf$hierarchy_keys_depth())
#' x <- xnested$data[[1]]
#' bb <- rlm_estimate(x,
#'   response = conf$get_response(),
#'   feature = feature,
#'   samples = conf$sampleName
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
#' @param config \code{\link{AnalysisConfiguration}}
#' @family aggregation
#' @keywords internal
#' @examples
#'
#' bb <- prolfqua_data("data_ionstar")$filtered()
#' conf <- old2new(bb$config)
#' data <- bb$data
#' conf$hierarchyDepth <- 1
#' xnested <- data |>
#'   dplyr::group_by(across(all_of(conf$hierarchy_keys_depth()))) |>
#'   tidyr::nest()
#'
#' feature <- base::setdiff(conf$hierarchy_keys(), conf$hierarchy_keys_depth())
#' x <- xnested$data[[1]]
#' bb <- rlm_estimate_dfconfig(x, conf)
#'
#' prolfqua:::.reestablish_condition(x, bb, conf)
#'
rlm_estimate_dfconfig <- function(pdata, config, name = FALSE) {
  if (name) {
    return("lmrob")
  }

  feature <- base::setdiff(config$hierarchy_keys(), config$hierarchy_keys_depth())
  rlm_estimate(pdata, response = config$get_response(), feature = feature, samples = config$sampleName, maxIt = 20)
}

#' Convert old proflqua configurations (prolfqua 0.4) to new Analysis configurations
#' prolfqua 0.5
#' @param config an old (as stored in the example data) AnalysisConfiguration
#' @return  a new AnalysisConfiguration
#' @export
#' @examples
#' dd <- prolfqua_data("data_ionstar")$filtered()
#' dd$config <- old2new(dd$config)
#'
old2new <- function(config) {
  newconfig <- AnalysisConfiguration$new()
  newconfig$factors <- config$factors
  newconfig$factorDepth <- config$factorDepth
  newconfig$workIntensity <- config$workIntensity
  newconfig$hierarchyDepth <- config$hierarchyDepth
  newconfig$hierarchy <- config$hierarchy
  newconfig$isotopeLabel <- config$isotopeLabel
  newconfig$is_response_transformed <- if (!is.null(config$is_response_transformed)) {
    config$is_response_transformed
  } else if (!is.null(config$is_intensity_transformed)) {
    config$is_intensity_transformed
  } else {
    FALSE
  }
  newconfig$ident_qValue <- config$ident_qValue
  newconfig$sampleName <- config$sampleName
  newconfig$isotopeLabel <- config$isotopeLabel
  newconfig$fileName <- config$fileName
  return(newconfig)
}

.add_nr_children <- function(data, aggregated_data, config, newconfig) {
  new_child <- paste(c("nr_children", config$hierarchy_keys_depth()), collapse = "_")
  res_nr_children <- nr_obs_sample(data, config, new_child = new_child)
  result <- inner_join(aggregated_data, res_nr_children, by = c(config$hierarchy_keys_depth(), config$fileName))
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
#' config <- dd$config
#' data <- dd$data
#' data <- prolfqua::transform_work_intensity(data, config, log2)
#' bbMed <- estimate_intensity(data, config, .func = medpolish_estimate_dfconfig)
#' bbRob <- estimate_intensity(data, config, .func = rlm_estimate_dfconfig)
#' nrow(bbMed$data)
#' nrow(bbRob$data)
#' length(bbMed$data$medpolish)
#' length(bbRob$data$lmrob)
#' xt <- dplyr::inner_join(bbMed$data, bbRob$data)
#' plot(xt$medpolish, xt$lmrob, log = "xy", pch = "*")
#' abline(0, 1, col = 2)
#'
estimate_intensity <- function(data, config, .func) {
  makeName <- .func(name = TRUE)
  config <- config$clone(deep = TRUE)

  xnested <- data |>
    group_by(across(all_of(config$hierarchy_keys_depth()))) |>
    nest()

  loopOverNested <- function(xnested, .func, config) {
    pb <- progress::progress_bar$new(total = nrow(xnested))
    message("starting aggregation")
    purrr::map(xnested$data, function(d) {
      pb$tick()
      aggr <- .func(d, config)
      .reestablish_condition(d, aggr, config)
    })
  }

  res <- loopOverNested(xnested, .func = .func, config = config)

  xnested[[makeName]] <- res
  newconfig <- make_reduced_hierarchy_config(
    config,
    workIntensity = .func(name = TRUE),
    hierarchy = config$hierarchy_keys_depth(names = FALSE)
  )

  unnested <- xnested |>
    dplyr::select(all_of(c(config$hierarchy_keys_depth(), makeName))) |>
    tidyr::unnest(cols = dplyr::all_of(makeName)) |>
    dplyr::ungroup()

  return(.add_nr_children(data, unnested, config, newconfig))
}


#' Plot feature data and result of aggregation
#'
#' @param data data.frame before aggregation
#' @param config AnalysisConfiguration
#' @param data_aggr data.frame after aggregation
#' @param config_reduced AnalysisConfiguration of aggregated data
#' @param show.legend logical, show legend in plot
#' @family plotting
#' @family aggregation
#' @keywords internal
#' @export
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' config <- istar$config
#' analysis <- istar$data
#'
#' analysis <- prolfqua::transform_work_intensity(analysis, config, log2)
#' bbMed <- estimate_intensity(analysis, config, .func = medpolish_estimate_dfconfig)
#' tmpMed <- plot_estimate(analysis, config, bbMed$data, bbMed$config)
#' stopifnot("ggplot" %in% class(tmpMed$plots[[1]]))
#' stopifnot("ggplot" %in% class(tmpMed$plots[[2]]))
#'
#' bbRob <- estimate_intensity(analysis, config, .func = rlm_estimate_dfconfig)
#' tmpRob <- plot_estimate(analysis, config, bbRob$data, bbRob$config)
#' stopifnot("ggplot" %in% class(tmpRob$plots[[1]]))
#' stopifnot("ggplot" %in% class(tmpRob$plots[[2]]))
#'
plot_estimate <- function(data, config, data_aggr, config_reduced, show.legend = FALSE) {
  hierarchy_ID <- "hierarchy_ID"
  xnested <- data |>
    group_by(!!!syms(config$hierarchy_keys_depth())) |>
    nest()
  xnested <- xnested |> tidyr::unite(hierarchy_ID, !!!syms(config$hierarchy_keys_depth()))
  xnested_aggr <- data_aggr |>
    group_by(!!!syms(config_reduced$hierarchy_keys_depth())) |>
    nest_by(.key = "other")
  xnested_aggr <- xnested_aggr |> tidyr::unite(hierarchy_ID, !!!syms(config$hierarchy_keys_depth()))
  xnested_all <- inner_join(xnested, xnested_aggr, by = hierarchy_ID)

  plots <- vector(mode = "list", length = nrow(xnested_all))

  pb <- progress::progress_bar$new(total = nrow(xnested_all))
  for (i in seq_len(nrow(xnested_all))) {
    p1 <- plot_hierarchies_line(
      xnested_all$data[[i]],
      xnested_all[[hierarchy_ID]][i],
      config = config,
      show.legend = show.legend
    )
    p2 <- plot_hierarchies_add_quantline(
      p1,
      xnested_all$other[[i]],
      config_reduced$get_response(),
      config
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
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param func function to use for aggregation
#' @param N default 3 top intensities.
#' @return list with data and new reduced configuration (config)
#' @family aggregation
#' @export
#' @keywords internal
#' @examples
#'
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' config <- dd$config
#' res <- dd$data
#' ranked <- rank_peptide_by_intensity(res, config)
#'
#' mean_f <- function(x, name = FALSE) {
#'   if (name) {
#'     return("mean")
#'   }
#'   mean(x, na.rm = TRUE)
#' }
#' sum_f <- function(x, name = FALSE) {
#'   if (name) {
#'     return("sum")
#'   }
#'   sum(x, na.rm = TRUE)
#' }
#'
#' resTOPN <- aggregate_intensity_topN(
#'   ranked,
#'   config,
#'   .func = mean_f,
#'   N = 3
#' )
#'
#' print(dim(resTOPN$data))
#' # stopifnot(dim(resTOPN$data) == c(3260, 8))
#' stopifnot(names(resTOPN) %in% c("data", "config"))
#' config$get_response()
#' tmpRob <- plot_estimate(ranked,
#'   config,
#'   resTOPN$data,
#'   resTOPN$config,
#'   show.legend = TRUE
#' )
#' stopifnot("ggplot" %in% class(tmpRob$plots[[4]]))
#'
aggregate_intensity_topN <- function(pdata, config, .func, N = 3) {
  xcall <- as.list(match.call())
  newcol <- make.names(paste0("srm_", .func(name = TRUE), "_", xcall$N))

  topInt <-
    pdata |> dplyr::filter(!!sym("srm_meanIntRank") <= N)

  topInt <- topInt |>
    dplyr::group_by(across(all_of(c(
      config$hierarchy_keys_depth(),
      config$sampleName,
      config$fileName,
      config$isotopeLabel,
      config$factor_keys()
    ))))
  sumTopInt <- topInt |>
    dplyr::summarize(
      !!newcol := .func(!!sym(config$get_response())),
      ident_qValue = min(!!sym(config$ident_qValue)),
      .groups = "drop"
    )

  newconfig <- make_reduced_hierarchy_config(
    config,
    workIntensity = newcol,
    hierarchy = config$hierarchy[seq_len(config$hierarchyDepth)]
  )

  return(.add_nr_children(pdata, sumTopInt, config, newconfig))
}
