#' Functions - Plotting ----
#'
#' Plot peptide and fragments
#' @param data long-format measurement data frame
#' @param protein_name protein identifier displayed as the plot title
#' @param sample sample identifier column
#' @param intensity intensity column
#' @param peptide peptide identifier column
#' @param fragment fragment identifier column
#' @param factor explanatory-variable column
#' @param isotope_label column identifying isotope labels
#' @param log_y use a logarithmic intensity axis
#' @param show.legend show the plot legend
#' @noRd
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
  p <- p + .theme_prolfqua_shared(grid = TRUE) + points + lines
  p <- p + facet_grid(cols = vars(!!!syms(factor)), scales = "free_x")
  p <- p + ggtitle(protein_name)
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top")
  if (log_y) {
    p <- p + scale_y_log10()
  }
  p
}

#' Plot peptide intensities of protein as a function of the sample and factor
#'
#' @export
#' @param res data.frame
#' @param protein_name title of plot
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param show.legend logical, show legend in plot
#' @return a ggplot object with one line per feature across samples, coloured by peptide, faceted by factor
#'   levels and titled `protein_name`
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
  res
}


#' Attach sample annotations to the per-sample estimates of all subjects
#'
#' Keeps the column order subject keys, sample annotation, estimates.
#' @param data long-format measurement data frame
#' @param estimates per-subject, per-sample estimates with the subject keys first
#' @param subject_keys columns identifying the subjects
#' @param sample_name column identifying samples
#' @param factor_keys configured explanatory-variable columns
#' @param file_name column identifying the input measurement files
#' @param isotope_label column identifying isotope labels
#' @noRd
.reestablish_condition <- function(data, estimates, subject_keys, sample_name, factor_keys, file_name, isotope_label) {
  annotation <- data |>
    dplyr::select(dplyr::all_of(c(sample_name, factor_keys, file_name, isotope_label))) |>
    dplyr::distinct()
  # several isotope labels share a sample, so one estimate can match several annotation rows
  res <- dplyr::inner_join(estimates, annotation, by = sample_name, relationship = "many-to-many")
  dplyr::select(res, dplyr::all_of(c(subject_keys, colnames(annotation))), dplyr::everything())
}


#' Tukey's median polish estimate of e.g. a protein from its feature (peptide) intensities, one row per sample.
#' @param pdata long-format measurement table
#' @param response name of the numeric response column
#' @param feature child-feature identifier column
#' @param sample_name column identifying samples
#' @param control `maxiter` and `eps` of \code{stats::medpolish}
#' @noRd
.medpolish_estimate <- function(pdata, response, feature, sample_name, control = .aggregation_control("medpolish")) {
  # feature x sample matrix filled by index; the polish does not depend on row or column order
  feature_id <- if (length(feature) == 1L) {
    pdata[[feature]]
  } else {
    do.call(paste, c(unname(as.list(pdata[feature])), sep = "\r"))
  }
  features <- unique(feature_id)
  samples <- unique(pdata[[sample_name]])
  mat <- matrix(NA_real_, length(features), length(samples), dimnames = list(NULL, as.character(samples)))
  mat[cbind(match(feature_id, features), match(pdata[[sample_name]], samples))] <- pdata[[response]]
  X <- medpolish(mat, na.rm = TRUE, trace.iter = FALSE, maxiter = control$maxiter, eps = control$eps)
  tibble(!!sample_name := names(X$col), medpolish = unname(X$col + X$overall))
}

#' Estimate protein intensities by robust peptide regression
#' @param pdata long-format measurement table
#' @param response name of the numeric response column
#' @param feature child-feature identifier column
#' @param samples name of the sample identifier column
#' @param control `maxit` and `acc` of \code{MASS::rlm}
#' @noRd
.rlm_estimate <- function(pdata, response, feature, samples = "samples", control = .aggregation_control("rlm")) {
  data <- pdata |>
    select(all_of(c(samples, feature, response))) |>
    na.omit() |>
    droplevels()
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
    dplyr::left_join(all_samples, per_sample, by = samples, relationship = "one-to-many")
  }

  ## If there is only 1 peptide for all samples use the response of that peptide directly
  if (length(unique(data[[feature]])) == 1L) {
    per_sample <- data |>
      dplyr::mutate(lmrob = .data[[response]], weights = 1) |>
      rename(!!expname := all_of(response))
    return(finalize(per_sample))
  }

  ## model-matrix breaks on factors with 1 level so make vector of ones (will be intercept)
  if (length(unique(data[[samples]])) == 1L) {
    per_sample <- data |>
      summarize(
        lmrob = mean(.data[[response]]),
        !!expname := mean(.data[[response]]),
        weights = 1,
        .by = all_of(samples)
      )
    return(finalize(per_sample))
  }

  ## sum contrast on peptide level so sample effect will be mean over all peptides instead of reference level
  formula <- as.formula(paste0("~ -1 + ", samples, " + ", feature))
  contr.arg <- list("contr.sum")
  names(contr.arg) <- feature
  X <- model.matrix(formula, data = data, contrasts.arg = contr.arg)
  ## MASS::rlm breaks on aliased columns: keep the columns of full rank found by the pivoted QR of the
  ## least-squares fit. .lm.fit() reports coefficients in pivoted order, so select by pivot and rank.
  y <- data[[response]]
  fit <- .lm.fit(X, y)
  X <- X[, sort(fit$pivot[seq_len(fit$rank)]), drop = FALSE]
  ## Last step is always rlm if X > has some columns left
  if (ncol(X) > 0) {
    fit <- MASS::rlm(X, y, maxit = control$maxit, acc = control$acc)
    data$residuals <- fit$residuals
    usamples <- unique(data[[samples]])
    coef_names <- paste0(samples, usamples)
    lmrob <- tibble(!!samples := usamples, lmrob = unname(fit$coefficients[coef_names]))

    sumdata <- data |>
      select(-all_of(feature)) |>
      dplyr::summarize(!!expname := mean(.data[[response]]), weights = 1 / mean(residuals^2), .by = all_of(samples))
    if (any(is.infinite(sumdata$weights) | is.na(sumdata$weights) | sumdata$weights > 10e6)) {
      sumdata$weights <- 1
    }

    res <- inner_join(sumdata, lmrob, by = samples, relationship = "one-to-one")
    if (sum(is.na(res[[expname]])) < sum(is.na(res$lmrob))) {
      res$lmrob <- res[[expname]]
    }
  } else {
    res <- data |>
      select(-all_of(feature)) |>
      dplyr::summarize(!!expname := mean(.data[[response]]), lmrob = mean(.data[[response]]), .by = all_of(samples))
    res$weights <- 1
  }
  finalize(res)
}

#' Count measured children per parent and sample after aggregation
#' @param data long-format measurement data frame
#' @param aggregated_data protein-level estimates returned by aggregation
#' @param hierarchy_keys_depth identifier columns at the output depth
#' @param file_name column identifying the input measurement files
#' @param response name of the numeric response column
#' @param nr_children_col column tracking measured child-feature counts
#' @param newconfig configuration of the aggregated output
#' @noRd
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
  result <- inner_join(
    aggregated_data,
    res_nr_children,
    by = c(hierarchy_keys_depth, file_name),
    relationship = "many-to-one"
  )
  newconfig$nr_children <- new_child
  list(data = result, config = newconfig)
}

#' Aggregates e.g. protein abundances from peptide abundances
#'
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param method "medpolish" (Tukey's median polish, column `medpolish`) or "rlm" (robust regression with
#'   `MASS::rlm`, column `lmrob`)
#' @param control convergence controls merged over the defaults. medpolish: `maxiter` (default 10) and `eps`
#'   (default 0.01, the change of the sum of absolute residuals that stops the polish, on the log2 scale of the
#'   input); rlm: `maxit` (default 20) and `acc` (default 1e-4, relative change of the IRLS residuals).
#' @return returns list with data (data.frame), config (AnalysisConfiguration) and failures (tibble with the
#'   subject `id` and `warning` text of every warning raised while estimating one subject; the warnings are not
#'   re-raised)
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
estimate_intensity <- function(lfqdata, method = c("medpolish", "rlm"), control = list()) {
  method <- match.arg(method)
  control <- .aggregation_control(method, control)
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
    function(d) .medpolish_estimate(d, response, feature, sname, control)
  } else {
    function(d) .rlm_estimate(unite(d, "feature", all_of(feature)), response, "feature", sname, control)
  }

  xnested <- data |>
    group_by(across(all_of(hkeysd))) |>
    nest()

  pb <- .make_progress(nrow(xnested))
  # per-subject warnings (e.g. rlm or medpolish not converging) are collected with the subject id and returned
  ids <- .key_ids(xnested, hkeysd)
  failures <- list()
  res <- purrr::map(seq_len(nrow(xnested)), function(i) {
    pb$tick()
    withCallingHandlers(
      estimate(xnested$data[[i]]),
      warning = function(w) {
        failures[[length(failures) + 1L]] <<- tibble::tibble(id = ids[i], warning = conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
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
  unnested <- .reestablish_condition(data, unnested, hkeysd, sname, fkeys, fname, iso)

  result <- .add_nr_children(data, unnested, hkeysd, fname, response, nrc, newconfig)
  result$failures <- dplyr::bind_rows(tibble::tibble(id = character(0), warning = character(0)), failures)
  result
}


#' Plot feature data and result of aggregation
#'
#' @param lfqdata LFQData object (original data)
#' @param lfqdata_agg LFQData object (aggregated data)
#' @param show.legend logical, show legend in plot
#' @return a tibble with one row per aggregated element and the columns `hierarchy_id`, `data` (the feature
#'   intensities), `other` (the aggregated intensities) and `plots` (ggplot objects of the feature intensities with
#'   the aggregate drawn in black)
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
  plot_data <- .estimate_plot_data(lfqdata, lfqdata_agg)
  estimate_aes <- .estimate_aes(lfqdata$sample_name(), lfqdata_agg$response())
  plots <- vector(mode = "list", length = nrow(plot_data))
  pb <- .make_progress(nrow(plot_data))
  for (i in seq_len(nrow(plot_data))) {
    plots[[i]] <- .estimate_plot(plot_data, i, lfqdata, estimate_aes, show.legend)
    pb$tick()
  }
  plot_data$plots <- plots
  plot_data
}

#' Feature data of each aggregated element next to its aggregate, one row per element
#' @param lfqdata LFQData with the feature data
#' @param lfqdata_agg LFQData with the aggregated data
#' @return tibble with `hierarchy_id`, `data` (feature intensities) and `other` (aggregated intensities)
#' @noRd
.estimate_plot_data <- function(lfqdata, lfqdata_agg) {
  hkeysd <- lfqdata$relevant_hierarchy_keys()
  xnested <- lfqdata$data_long() |>
    group_by(across(all_of(hkeysd))) |>
    nest()
  xnested <- xnested |> tidyr::unite("hierarchy_id", all_of(hkeysd))
  xnested_aggr <- lfqdata_agg$data_long() |>
    group_by(across(all_of(lfqdata_agg$relevant_hierarchy_keys()))) |>
    nest_by(.key = "other")
  xnested_aggr <- xnested_aggr |> tidyr::unite("hierarchy_id", all_of(hkeysd))
  inner_join(xnested, xnested_aggr, by = "hierarchy_id", relationship = "one-to-one")
}

#' Mapping of the aggregate line; built here so its environment holds the two column names and nothing else
#' @param sample_name sample column
#' @param response aggregated response column
#' @noRd
.estimate_aes <- function(sample_name, response) {
  aes(x = .data[[sample_name]], y = .data[[response]], group = 1)
}

#' Plot of element `i`: its feature intensities with the aggregate drawn in black
#' @param plot_data result of .estimate_plot_data()
#' @param i row of `plot_data`
#' @param lfqdata LFQData with the feature data
#' @param estimate_aes result of .estimate_aes()
#' @param show.legend show the legend
#' @noRd
.estimate_plot <- function(plot_data, i, lfqdata, estimate_aes, show.legend) {
  p <- plot_hierarchies_line(
    plot_data$data[[i]],
    plot_data$hierarchy_id[i],
    lfqdata = lfqdata,
    show.legend = show.legend
  )
  p +
    geom_line(data = plot_data$other[[i]], estimate_aes, linewidth = 1.3, color = "black", linetype = "solid") +
    geom_point(data = plot_data$other[[i]], estimate_aes, color = "black", shape = 10)
}


#' Aggregates top N intensities
#'
#' run \link{rank_peptide_by_intensity} first
#' @param ranked_data data.frame with ranked peptides
#' @param lfqdata [LFQData] object with the data and its configuration
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
    ranked_data |> dplyr::filter(.data$srm_meanIntRank <= .env$N)

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
      !!newcol := .func(.data[[lfqdata$response()]]),
      # an unknown q-value of a top feature leaves the aggregate's q-value unknown
      ident_qValue = min(.data[[config$ident_q_value]], na.rm = FALSE),
      .groups = "drop"
    )

  newconfig <- make_reduced_hierarchy_config(
    config,
    work_intensity = newcol,
    hierarchy = config$hierarchy[seq_len(config$hierarchy_depth)]
  )

  .add_nr_children(
    ranked_data,
    sum_top_intensities,
    lfqdata$relevant_hierarchy_keys(),
    lfqdata$file_name(),
    lfqdata$response(),
    lfqdata$nr_children_col(),
    newconfig
  )
}


#' Sum of observed children per hierarchy unit and sample
#'
#' For each hierarchy unit and sample, sums the `nr_children_col` column over
#' the rows with a non-missing response, e.g. the number of peptides behind
#' each protein intensity of a sample. Counts per sample, observed values only.
#'
#' @param data data.frame
#' @param response character - intensity column name
#' @param hierarchy_keys_depth character vector - hierarchy columns at current depth
#' @param file_name character - file name column
#' @param nr_children_col character - nr_children column name
#' @param new_child character - output column name
#' @return The computed result.
#' @export
#' @family hierarchy counting
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
    group_by(across(all_of(c(hierarchy_keys_depth, file_name)))) |>
    summarize(!!new_child := sum(.data[[nr_children_col]], na.rm = TRUE), .groups = "drop")
  nr_children
}


#' Max nr_children across samples per hierarchy unit
#'
#' Aggregates nr_children per sample (via nr_obs_sample), then takes the max per hierarchy unit.
#' Counts experiment-wide, observed values only.
#'
#' @param data data.frame
#' @param response character - intensity column name
#' @param hierarchy_keys_depth character vector - hierarchy columns at current depth
#' @param file_name character - file name column
#' @param nr_children_col character - nr_children column name
#' @param name_nr_child character - output column name
#' @return The computed result.
#' @export
#' @family hierarchy counting
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
  # experiment-wide count = the largest per-sample count across all samples; counts are never NA
  per_sample |>
    dplyr::group_by(across(all_of(hierarchy_keys_depth))) |>
    dplyr::summarize(!!name_nr_child := max(.data[[per_sample_col]]), .groups = "drop")
}

#' Count distinct child features per hierarchy unit
#'
#' Counts the number of distinct entries of the hierarchy levels below
#' `hierarchy_keys_depth` (e.g. peptides per protein). Counts experiment-wide,
#' whether observed or not.
#'
#' @param data data.frame
#' @param hierarchy_keys character vector - all hierarchy column names
#' @param hierarchy_keys_depth character vector - hierarchy columns at current depth
#' @param name_nr_child character - output column name
#' @return The computed result.
#' @export
#' @family hierarchy counting
#' @examples
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(dd$data, dd$config)
#' xd <- nr_features_experiment(lfq$data_long(), lfq$hierarchy_keys(),
#'   lfq$relevant_hierarchy_keys())
#' stopifnot(min(xd$nr_features_exp) == 1)
#'
nr_features_experiment <- function(data, hierarchy_keys, hierarchy_keys_depth, name_nr_child = "nr_features_exp") {
  nr_B_in_A(data, hierarchy_keys_depth, setdiff(hierarchy_keys, hierarchy_keys_depth), name_nr_child)
}

# Summarize Intensities by Intensity or NAs ----
#' Rank precursors or peptides by intensity
#' @param pdata data.frame
#' @param response character - intensity column name
#' @param hierarchy_keys character vector - all hierarchy column names
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
    dplyr::group_by(across(all_of(hierarchy_keys))) |>
    dplyr::summarize(!!summary_column := mean(.data[[response]], na.rm = TRUE)) |>
    dplyr::arrange(.data[[hierarchy_keys[1]]]) |>
    dplyr::group_by(.data[[hierarchy_keys[1]]]) |>
    dplyr::mutate(!!rank_column := min_rank(desc(.data[[summary_column]])))
  pdata <- dplyr::inner_join(pdata, ranked, by = hierarchy_keys, relationship = "many-to-one")

  inform_prolfqua("Columns added : ", summary_column, " ", rank_column)
  pdata
}


#' Convergence controls of an aggregation method, merged over its defaults
#' @param method "medpolish" or "rlm"
#' @param control named list of user controls
#' @noRd
.aggregation_control <- function(method, control = list()) {
  defaults <- switch(method, medpolish = list(maxiter = 10, eps = 0.01), rlm = list(maxit = 20, acc = 1e-4))
  .merge_control(control, defaults)
}
