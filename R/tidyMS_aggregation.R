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
  log_y = FALSE,
  show.legend = FALSE
) {
  # missing intensities are expected; lines still break at them
  points <- geom_point(show.legend = show.legend, na.rm = TRUE)
  lines <- geom_line(show.legend = show.legend, na.rm = TRUE)
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
#' @param res data.frame
#' @param protein_name title of plot
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param show.legend logical, show legend in plot
#' @return a ggplot object with one line per feature across samples, coloured by peptide, faceted by factor
#'   levels and titled `protein_name`
#' @noRd
.plot_hierarchies_line <- function(res, protein_name, lfqdata, show.legend = FALSE) {
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
#' @noRd
.reestablish_condition <- function(data, estimates, subject_keys, sample_name, factor_keys, file_name) {
  annotation <- data |>
    dplyr::select(dplyr::all_of(c(sample_name, factor_keys, file_name))) |>
    dplyr::distinct()
  res <- dplyr::inner_join(estimates, annotation, by = sample_name, relationship = "many-to-one")
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
#'
#' Base R per subject (the subject loop calls it thousands of times): `MASS::rlm()` of the response on sample and
#' feature effects, feature effects with sum contrasts, so a sample effect is the mean over the features.
#' @param pdata long-format measurement table
#' @param response name of the numeric response column
#' @param feature child-feature identifier column(s); several columns are joined with `"_"`, as `tidyr::unite()`
#' @param samples name of the sample identifier column
#' @param control `maxit` and `acc` of \code{MASS::rlm}
#' @return one row per sample of `pdata`, in order of appearance: `samples`, `mean.<response>`, `lmrob` and
#'   `weights`; `NA` for samples without an observed response
#' @noRd
.rlm_estimate <- function(pdata, response, feature, samples = "samples", control = .aggregation_control("rlm")) {
  sample_id <- pdata[[samples]]
  observed <- !is.na(pdata[[response]]) & !is.na(sample_id)
  data <- data.frame(
    sample_id[observed],
    do.call(paste, c(unname(as.list(pdata[feature])), sep = "_"))[observed],
    stringsAsFactors = FALSE
  )
  names(data) <- c(samples, "feature")
  data <- droplevels(data)
  y <- pdata[[response]][observed]
  usamples <- unique(data[[samples]])
  sid <- match(data[[samples]], usamples)
  per_sample_mean <- function(x) vapply(split(x, sid), mean, numeric(1), USE.NAMES = FALSE)
  mean_y <- per_sample_mean(y)

  ## All branches hand off here so the output schema is identical: a single, fixed column order
  ## (samples, mean.<response>, lmrob, weights) with one row per sample present in the input, also for samples
  ## without an observed response, so downstream unnest()/.reestablish_condition() always bind the same columns.
  all_samples <- unique(sample_id)
  finalize <- function(lmrob, weights) {
    idx <- match(all_samples, usamples)
    tibble(
      !!samples := all_samples,
      !!paste0("mean.", response) := mean_y[idx],
      lmrob = lmrob[idx],
      weights = weights[idx]
    )
  }

  ## If there is only 1 peptide for all samples, or only one sample, use the mean response per sample directly
  ## (model.matrix breaks on factors with 1 level)
  if (length(unique(data$feature)) == 1L || length(usamples) == 1L) {
    return(finalize(mean_y, rep(1, length(usamples))))
  }

  ## sum contrast on peptide level so sample effect will be mean over all peptides instead of reference level
  formula <- as.formula(paste0("~ -1 + ", samples, " + feature"))
  X <- model.matrix(formula, data = data, contrasts.arg = list(feature = "contr.sum"))
  ## MASS::rlm breaks on aliased columns: keep the columns of full rank found by the pivoted QR of the
  ## least-squares fit. .lm.fit() reports coefficients in pivoted order, so select by pivot and rank; the sample
  ## columns keep X of rank at least one.
  fit <- .lm.fit(X, y)
  X <- X[, sort(fit$pivot[seq_len(fit$rank)]), drop = FALSE]
  fit <- MASS::rlm(X, y, maxit = control$maxit, acc = control$acc)
  lmrob <- unname(fit$coefficients[paste0(samples, usamples)])
  weights <- 1 / per_sample_mean(fit$residuals^2)
  if (any(is.infinite(weights) | is.na(weights) | weights > 10e6)) {
    weights <- rep(1, length(usamples))
  }
  if (sum(is.na(mean_y)) < sum(is.na(lmrob))) {
    lmrob <- mean_y
  }
  finalize(lmrob, weights)
}

#' Count measured children per parent and sample after aggregation
#' @param data long-format measurement data frame
#' @param aggregated_data protein-level estimates returned by aggregation
#' @param hierarchy_keys_depth identifier columns at the output depth
#' @param file_name column identifying the input measurement files
#' @param response name of the numeric response column
#' @param nr_children_col column tracking measured child-feature counts
#' @param confidence_col column with the identification confidence of the children
#' @param newconfig configuration of the aggregated output
#' @return list with `data`, the aggregated data with the per-sample child count, the identification confidence
#'   (maximum over the observed children) and both feature-level columns, and `config`
#' @noRd
.add_nr_children <- function(
  data,
  aggregated_data,
  hierarchy_keys_depth,
  file_name,
  response,
  nr_children_col,
  confidence_col,
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
  confidence <- .aggregate_confidence(
    data,
    response,
    c(hierarchy_keys_depth, file_name),
    confidence_col,
    newconfig$identification_confidence
  )
  res_nr_children <- dplyr::inner_join(
    res_nr_children,
    confidence,
    by = c(hierarchy_keys_depth, file_name),
    relationship = "one-to-one"
  )
  result <- inner_join(
    aggregated_data,
    res_nr_children,
    by = c(hierarchy_keys_depth, file_name),
    relationship = "many-to-one"
  )
  newconfig$nr_children <- new_child
  list(data = .add_feature_columns(result, hierarchy_keys_depth, newconfig), config = newconfig)
}

#' Identification confidence of aggregated features: the maximum over the observed children
#' @param data long-format child-level data
#' @param response name of the numeric response column
#' @param by columns identifying one aggregated feature in one sample
#' @param confidence_col identification confidence column of the children
#' @param new_col name of the output column
#' @return a tibble with the `by` columns and `new_col`, one row per aggregated feature and sample with an observed
#'   child
#' @noRd
.aggregate_confidence <- function(data, response, by, confidence_col, new_col) {
  observed <- data[!is.na(data[[response]]), , drop = FALSE]
  dplyr::summarise(
    observed,
    !!new_col := .max_or_na(.data[[confidence_col]]),
    .by = dplyr::all_of(by)
  )
}

#' Add the feature-level identification confidence and child count: the maximum over samples
#' @param data long-format data with the per-sample columns named in `config`
#' @param keys columns identifying one feature
#' @param config configuration naming the per-sample and feature-level columns
#' @noRd
.add_feature_columns <- function(data, keys, config) {
  dplyr::mutate(
    data,
    !!config$identification_confidence_feature := .max_or_na(.data[[config$identification_confidence]]),
    !!config$nr_children_feature := .max_or_na(.data[[config$nr_children]]),
    .by = dplyr::all_of(keys)
  )
}

#' Aggregates e.g. protein abundances from peptide abundances
#'
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param method "medpolish" (Tukey's median polish, column `medpolish`) or "rlm" (robust regression with
#'   `MASS::rlm`, column `lmrob`)
#' @param control convergence controls merged over the defaults. medpolish: `maxiter` (default 10) and `eps`
#'   (default 0.01; the polish stops when the sum of absolute residuals changes by less than `eps` times that sum,
#'   a relative criterion independent of the input scale); rlm: `maxit` (default 20) and `acc` (default 1e-4,
#'   relative change of the IRLS residuals).
#' @return returns list with data (data.frame), config (AnalysisConfiguration) and failures (tibble with the
#'   subject `id` and `warning` text of every warning raised while estimating one subject; the warnings are not
#'   re-raised)
#' @importFrom MASS rlm
#' @noRd
.estimate_intensity <- function(lfqdata, method = c("medpolish", "rlm"), control = list()) {
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
  nrc <- lfqdata$nr_children_col()
  feature <- base::setdiff(hkeys, hkeysd)
  estimate <- if (method == "medpolish") {
    function(d) .medpolish_estimate(d, response, feature, sname, control)
  } else {
    function(d) .rlm_estimate(d, response, feature, sname, control)
  }

  xnested <- data |>
    group_by(across(all_of(hkeysd))) |>
    nest()

  pb <- .make_progress(nrow(xnested))
  # per-subject warnings (e.g. rlm or medpolish not converging) are collected with the subject id and returned
  ids <- .key_ids(xnested, hkeysd)
  fits <- purrr::map(xnested$data, function(d) {
    pb$tick()
    .fit_capturing_warnings(estimate, d)
  })
  subject_warnings <- purrr::map(fits, "warnings")
  failures <- tibble::tibble(
    id = rep(ids, lengths(subject_warnings)),
    warning = as.character(unlist(subject_warnings))
  )

  xnested[[make_name]] <- purrr::map(fits, "fit")
  newconfig <- .make_reduced_hierarchy_config(
    config,
    work_intensity = make_name,
    hierarchy = config$hierarchy_keys_depth(names = FALSE)
  )

  unnested <- xnested |>
    dplyr::select(all_of(c(hkeysd, make_name))) |>
    tidyr::unnest(cols = dplyr::all_of(make_name)) |>
    dplyr::ungroup()
  unnested <- .reestablish_condition(data, unnested, hkeysd, sname, fkeys, fname)

  result <- .add_nr_children(data, unnested, hkeysd, fname, response, nrc, config$identification_confidence, newconfig)
  result$failures <- failures
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
#' @noRd
.plot_estimate <- function(lfqdata, lfqdata_agg, show.legend = FALSE) {
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
  p <- .plot_hierarchies_line(
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
#' @noRd
.aggregate_intensity_top_n <- function(ranked_data, lfqdata, .func, N = 3) {
  newcol <- make.names(paste0("srm_", .func(name = TRUE), "_", N))
  config <- lfqdata$get_config()$clone(deep = TRUE)

  top_intensities <-
    ranked_data |> dplyr::filter(.data$srm_meanIntRank <= .env$N)

  top_intensities <- top_intensities |>
    dplyr::group_by(across(all_of(c(
      lfqdata$relevant_hierarchy_keys(),
      lfqdata$sample_name(),
      lfqdata$file_name(),
      lfqdata$factor_keys()
    ))))
  sum_top_intensities <- top_intensities |>
    dplyr::summarize(
      !!newcol := .func(.data[[lfqdata$response()]]),
      .groups = "drop"
    )

  newconfig <- .make_reduced_hierarchy_config(
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
    config$identification_confidence,
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
#' @return a tibble with one row per hierarchy element and sample: `hierarchy_keys_depth`, `file_name` and
#'   `new_child`, the sum of `nr_children_col` over rows with an observed response.
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
  observed <- !is.na(data[[response]])
  data <- data[observed, unique(c(hierarchy_keys_depth, file_name, nr_children_col)), drop = FALSE]
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
#' @return a tibble with one row per hierarchy element: `hierarchy_keys_depth` and `name_nr_child`, the largest
#'   per-sample count of observed children.
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
#' @return a tibble with one row per hierarchy element: `hierarchy_keys_depth` and `name_nr_child`, the number
#'   of distinct child features, observed or not.
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
#' @noRd
.rank_peptide_by_intensity <- function(pdata, response, hierarchy_keys) {
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
