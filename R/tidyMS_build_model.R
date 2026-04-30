#' Likelihood ratio test
#' @family modelling
#' @return The computed result.
#' @export
#' @param complete_models table with models (see build model)
#' @param model_name name of model
#' @param complete_models_int reduced model
#' @param model_name_int name of reduced model
#' @param subject_id subject id typically Assession or protein_Id
#' @param path default NULL, set to a directory if you need to write diagnostic plots.
#' @examples
#' data_2Factor <- prolfqua::sim_lfq_data_2factor_config(
#'  Nprot = 200,
#'  with_missing = TRUE,
#'  weight_missing = 2)
#'
#' pMerged <- LFQData$new(data_2Factor$data, data_2Factor$config)
#'
#' pMerged$response()
#' pMerged$factors()
#'
#' formula_condition_and_Batches <-
#'   prolfqua::strategy_lm("abundance ~ Treatment + Background")
#' modCB <- prolfqua::build_model(
#'   pMerged,
#'   formula_condition_and_Batches)
#'
#' formula_condition <-
#'   prolfqua::strategy_lm("abundance ~ Treatment")
#' modC <- prolfqua::build_model(
#'   pMerged,
#'   formula_condition)
#'
#' tmp <- LR_test(modCB$model_df, "modCB", modC$model_df, "modB")
#' hist(tmp$likelihood_ratio_test.pValue)
#'
LR_test <- function(
  complete_models,
  model_name,
  complete_models_int,
  model_name_int,
  subject_id = "protein_Id",
  path = NULL
) {
  # Model Comparison
  reg <- dplyr::inner_join(
    dplyr::select(complete_models, !!sym(subject_id), "linear_model"),
    dplyr::select(complete_models_int, !!sym(subject_id), "linear_model"),
    by = subject_id
  )

  reg <- reg |>
    dplyr::mutate(
      modelComparisonLikelihoodRatioTest = map2(
        !!sym("linear_model.x"),
        !!sym("linear_model.y"),
        .likelihood_ratio_test
      )
    )
  likelihood_ratio_test_result <- reg |>
    dplyr::select(!!sym(subject_id), dplyr::all_of("modelComparisonLikelihoodRatioTest")) |>
    tidyr::unnest(cols = c("modelComparisonLikelihoodRatioTest"))
  likelihood_ratio_test_result <- likelihood_ratio_test_result |>
    dplyr::rename(likelihood_ratio_test.pValue = .data$modelComparisonLikelihoodRatioTest)

  if (!is.null(path)) {
    plot_lrt_diagnostics(likelihood_ratio_test_result, model_name, model_name_int, path)
  }

  return(likelihood_ratio_test_result)
}

#' Write LRT diagnostic plots to PDF
#'
#' Writes a histogram of p-values and an empirical CDF to a PDF file.
#' Called by \code{\link{LR_test}} when \code{path} is non-NULL.
#'
#' @param result tibble with a \code{likelihood_ratio_test.pValue} column
#' @param model_name name of the full model
#' @param model_name_int name of the reduced/interaction model
#' @param path directory to write the PDF into
#' @keywords internal
plot_lrt_diagnostics <- function(result, model_name, model_name_int, path) {
  file_name <- paste("hist_LRT_", model_name, "_", model_name_int, ".pdf", sep = "")
  file_name <- file.path(path, file_name)
  message("writing figure : ", file_name, "\n")
  pdf(file_name)
  par(mfrow = c(2, 1))
  hist(result$likelihood_ratio_test.pValue, breaks = 20)
  plot(ecdf(result$likelihood_ratio_test.pValue))
  abline(v = c(0.01, 0.05), col = c(3, 2))
  dev.off()
}


#' Build protein models from data
#'
#' @param data data - a data frame or LFQData object
#' @param model_strategy model strategy object (e.g. from strategy_lmer or strategy_lm)
#' @param subject_id grouping variable
#' @param model_name model name
#' @return
#' a object of class \code{\link{Model}}
#' @family modelling
#' @seealso \code{\link{model_analyse}}, \code{\link{strategy_lmer}} \code{\link{strategy_lm}}
#'
#' @export
#' @examples
#' D <- prolfqua::sim_lfq_data_peptide_config(Nprot = 20, weight_missing = 0.1)
#' D$data$abundance |> is.na() |> sum()
#' D <- prolfqua::sim_lfq_data_peptide_config(Nprot = 20, weight_missing = 0.1, seed =3)
#' D$data$abundance |> is.na() |> sum()
#' model_name <- "f_condtion_r_peptide"
#' formula_randomPeptide <-
#'   strategy_lmer("abundance  ~ group_ + (1 | peptide_Id) + (1 | sampleName)",
#'    model_name = model_name)
#'
#'
#' mod <- prolfqua::build_model(
#'  D$data,
#'  formula_randomPeptide,
#'  model_name = model_name,
#'  subject_id = D$config$hierarchy_keys_depth())
#' aovtable <- mod$get_anova()
#'
#' mod <- prolfqua::build_model(
#'  LFQData$new(D$data, D$config),
#'  formula_randomPeptide,
#'  model_name = model_name)
#' model_summary(mod)
#'
#'
build_model <- function(
  data,
  model_strategy,
  subject_id = if ("LFQData" %in% class(data)) {
    data$subject_id()
  } else {
    "protein_Id"
  },
  model_name = model_strategy$model_name
) {
  nested_data <- if ("LFQData" %in% class(data)) data$data_long() else data
  modelling_result <- model_analyse(nested_data, model_strategy, model_name = model_name, subject_id = subject_id)
  return(Model$new(
    model_df = modelling_result$model_df,
    model_strategy = model_strategy,
    model_name = modelling_result$model_name,
    subject_id = subject_id
  ))
}


#' Build protein models with LOD imputation for failed fits
#'
#' Fits per-protein models, then re-fits failed/singular proteins after
#' imputing missing values with the limit of detection (LOD) and clamping.
#' Covariance is borrowed from successful fits so that variance is not
#' underestimated by the constant imputation.
#'
#' @param lfqdata LFQData object (aggregated to protein level)
#' @param model_strategy model strategy object (e.g. from strategy_lm)
#' @param model_name model name (default appends "Imputed")
#' @param lod numeric limit of detection; if NULL, auto-computed from data
#' @param borrow_method "sigma" borrows scalar sigma and uses per-protein
#'   (X'X)^-1; "vcov" borrows element-wise median of full vcov matrices
#' @param df_method "observed" uses max(n_observed - p, 1);
#'   "borrowed" uses median df from successful fits
#' @return a object of class \code{\link{Model}}
#' @family modelling
#' @seealso \code{\link{build_model}}, \code{\link{impute_refit_singular}}
#' @export
#' @examples
#' istar <- sim_lfq_data_protein_config(Nprot = 30, weight_missing = 0.5)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$rename_response("transformedIntensity")
#' strat <- strategy_lm(paste(lfqdata$response(), "~ group_"))
#' mod <- build_model_impute(lfqdata, strat)
#'
build_model_impute <- function(
  lfqdata,
  model_strategy,
  model_name = paste0(model_strategy$model_name, "Imputed"),
  lod = NULL,
  borrow_method = c("sigma", "vcov"),
  df_method = c("observed", "borrowed")
) {
  borrow_method <- match.arg(borrow_method)
  df_method <- match.arg(df_method)
  subject_id <- lfqdata$subject_id()
  response <- lfqdata$response()

  modelling_result <- model_analyse(
    lfqdata$data_long(),
    model_strategy,
    model_name = model_name,
    subject_id = subject_id
  )

  if (is.null(lod)) {
    mh <- MissingHelpers$new(lfqdata$data_long(), lfqdata$get_config())
    lod <- mh$get_lod()
  }

  # Build sample template from annotation columns only (fileName, sampleName,
  # factors). These don't vary per protein, so distinct() gives one row per sample.
  # Value columns like nr_children are protein-specific and come from dat via join.
  annotation_cols <- intersect(
    lfqdata$get_config()$annotation_vars(),
    colnames(lfqdata$data_long())
  )
  sample_template <- lfqdata$data_long() |>
    dplyr::select(dplyr::all_of(annotation_cols)) |>
    dplyr::distinct()

  # Resolve nr_children column name (used as lm weights; must be filled for imputed rows)
  nr_children_col <- lfqdata$nr_children_col()
  if (!is.null(nr_children_col) && !nr_children_col %in% colnames(lfqdata$data_long())) {
    nr_children_col <- NULL
  }

  modelling_result$model_df <- impute_refit_singular(
    modelling_result$model_df,
    model_strategy,
    lod = lod,
    response = response,
    sample_template = sample_template,
    borrow_method = borrow_method,
    df_method = df_method,
    nr_children_col = nr_children_col
  )

  return(Model$new(
    model_df = modelling_result$model_df,
    model_strategy = model_strategy,
    model_name = modelling_result$model_name,
    subject_id = subject_id
  ))
}


#' Summarize modelling and error reporting
#' @param mod model table see \code{\link{build_model}}
#' @keywords internal
#' @family modelling
#' @export
#' @examples
#' D <- prolfqua::sim_lfq_data_peptide_config(Nprot = 20, weight_missing = 0.1)
#' formula_rp <- strategy_lmer("abundance ~ group_ + (1 | peptide_Id) + (1 | sampleName)")
#' mod <- prolfqua::build_model(
#'   LFQData$new(D$data, D$config), formula_rp)
#' res <- model_summary(mod)
#' stopifnot(is.list(res))
#' stopifnot(all(c("exists", "isSingular") %in% names(res)))
model_summary <- function(mod) {
  res <- list()
  res$exists <- table(mod$model_df$has_model_fit)
  res$isSingular <- table(mod$model_df$isSingular)
  return(res)
}


# Model fitting internals ----

.likelihood_ratio_test <- function(modelNO, model) {
  res <- tryCatch(anova(modelNO, model), error = function(x) NULL)
  if (!is.null(res)) {
    p_value_col <- grep("^Pr\\(", colnames(res), value = TRUE)
    if (length(p_value_col) == 0 || nrow(res) < 2) {
      return(NA_real_)
    }
    return(as.numeric(res[[p_value_col[1]]][2]))
  } else {
    return(NA_real_)
  }
}


# lm_imputed S3 class ----

#' Create an imputed lm wrapper
#'
#' Wraps an lm object with borrowed covariance information.
#' S3 generics vcov(), sigma(), and df.residual() dispatch
#' to the borrowed values instead of the original model's.
#'
#' @param model lm object fitted on imputed data
#' @param borrowed_vcov matrix, borrowed variance-covariance matrix
#' @param borrowed_sigma numeric, borrowed residual standard error
#' @param borrowed_df numeric, borrowed residual degrees of freedom
#' @param n_observed integer, number of non-imputed observations
#' @return lm_imputed object
#' @keywords internal
#' @family modelling
#' @examples
#' # Fit a normal lm, then wrap it with borrowed covariance
#' dat <- data.frame(group_ = rep(c("A", "B"), each = 4),
#'                   y = c(20.1, 20.5, 19.8, 20.3, 22.1, 22.4, 21.9, 22.2))
#' fit <- lm(y ~ group_, data = dat)
#'
#' # Wrap with borrowed variance (in practice these come from donor pool)
#' wrapped <- prolfqua:::new_lm_imputed(fit,
#'   borrowed_vcov = vcov(fit),
#'   borrowed_sigma = 0.8,
#'   borrowed_df = 6,
#'   n_observed = 5)
#'
#' # S3 dispatch returns borrowed values
#' stopifnot(inherits(wrapped, "lm_imputed"))
#' stopifnot(sigma(wrapped) == 0.8)
#' stopifnot(df.residual(wrapped) == 6)
#' # coefficients() still dispatches to underlying lm
#' stopifnot(identical(coefficients(wrapped), coefficients(fit)))
new_lm_imputed <- function(model, borrowed_vcov, borrowed_sigma, borrowed_df, n_observed) {
  attr(model, "borrowed_vcov") <- borrowed_vcov
  attr(model, "borrowed_sigma") <- borrowed_sigma
  attr(model, "borrowed_df") <- borrowed_df
  attr(model, "n_observed") <- n_observed
  class(model) <- c("lm_imputed", class(model))
  model
}

#' @export
vcov.lm_imputed <- function(object, ...) {
  attr(object, "borrowed_vcov")
}

#' @export
sigma.lm_imputed <- function(object, ...) {
  attr(object, "borrowed_sigma")
}

#' @export
df.residual.lm_imputed <- function(object, ...) {
  attr(object, "borrowed_df")
}


# Covariance borrowing helpers ----

#' Compute borrowed variance from successful model fits
#'
#' @param model_df tibble from model_analyse
#' @param method "sigma" borrows scalar sigma and uses per-protein (X'X)^-1,
#'   "vcov" borrows element-wise median of full vcov matrices
#' @return list with sigma, df, method, and optionally vcov
#' @keywords internal
#' @family modelling
#' @examples
#' mod <- sim_build_models_lm(model = "parallel3", weight_missing = 1)
#'
#' # Sigma method (default): returns median sigma and df from donors
#' borrowed_s <- prolfqua:::compute_borrowed_variance(
#'   mod$model_df, method = "sigma")
#' stopifnot(borrowed_s$method == "sigma")
#' stopifnot(is.numeric(borrowed_s$sigma) && borrowed_s$sigma > 0)
#' stopifnot(is.numeric(borrowed_s$df) && borrowed_s$df > 0)
#'
#' # Vcov method: element-wise median vcov from donors.
#' # Falls back to sigma if donor models have different coefficient counts.
#' mod_no_missing <- sim_build_models_lm(model = "parallel3",
#'   Nprot = 10, with_missing = FALSE)
#' borrowed_v <- prolfqua:::compute_borrowed_variance(
#'   mod_no_missing$model_df, method = "vcov")
#' stopifnot(borrowed_v$method == "vcov")
#' stopifnot(is.matrix(borrowed_v$vcov))
compute_borrowed_variance <- function(model_df, method = c("sigma", "vcov")) {
  method <- match.arg(method)
  good <- get_complete_model_fit(model_df)
  good <- good |> dplyr::filter(.data$isSingular == FALSE)

  if (nrow(good) == 0) {
    stop("No successful model fits available to borrow variance from.", call. = FALSE)
  }

  borrowed_sigma <- stats::median(good$sigma, na.rm = TRUE)
  borrowed_df <- stats::median(good$df.residual, na.rm = TRUE)

  if (method == "sigma") {
    return(list(sigma = borrowed_sigma, df = borrowed_df, method = "sigma"))
  }

  vcov_list <- lapply(good$linear_model, stats::vcov)
  ref_dim <- dim(vcov_list[[1]])
  ref_names <- dimnames(vcov_list[[1]])
  # Check all have same dimensions
  dims_ok <- all(vapply(vcov_list, function(v) identical(dim(v), ref_dim), logical(1)))
  if (!dims_ok) {
    warning("vcov dimensions differ across successful fits; falling back to sigma method.")
    return(list(sigma = borrowed_sigma, df = borrowed_df, method = "sigma"))
  }
  vcov_array <- array(unlist(vcov_list), dim = c(ref_dim, length(vcov_list)))
  borrowed_vcov <- apply(vcov_array, c(1, 2), stats::median)
  dimnames(borrowed_vcov) <- ref_names

  return(list(vcov = borrowed_vcov, sigma = borrowed_sigma, df = borrowed_df, method = "vcov"))
}


#' Impute and refit singular/failed models
#'
#' For proteins where the initial lm fit failed or produced NA coefficients,
#' impute missing values with LOD, clamp, refit, and attach borrowed covariance.
#'
#' @param model_df tibble from model_analyse
#' @param model_strategy strategy list from strategy_lm etc.
#' @param lod numeric, limit of detection value
#' @param response character, response column name in nested data
#' @param sample_template data.frame with all sample/group combinations
#'   (columns matching the nested data minus the response). Used to complete
#'   proteins that are entirely missing in one or more groups.
#' @param borrow_method "sigma" or "vcov"
#' @param df_method "observed" uses max(n_observed - p, 1),
#'   "borrowed" uses median df from successful fits
#' @param nr_children_col optional column name for nr_children (peptide counts);
#'   NA values in this column are filled with 1 for imputed rows so that
#'   weighted lm fits do not fail
#' @return modified model_df with imputed models replacing failed/singular ones
#' @keywords internal
#' @family modelling
impute_refit_singular <- function(
  model_df,
  model_strategy,
  lod,
  response,
  sample_template,
  borrow_method = c("sigma", "vcov"),
  df_method = c("observed", "borrowed"),
  nr_children_col = NULL
) {
  borrow_method <- match.arg(borrow_method)
  df_method <- match.arg(df_method)

  max_coef <- max(model_df$nr_coef, na.rm = TRUE)

  needs_impute <- (!model_df$has_model_fit) |
    (!is.na(model_df$isSingular) & model_df$isSingular) |
    (!is.na(model_df$nr_coef) & model_df$nr_coef < max_coef)

  if (!any(needs_impute)) {
    return(model_df)
  }

  borrowed <- compute_borrowed_variance(model_df, method = borrow_method)

  impute_idx <- which(needs_impute)
  results <- vector("list", length(impute_idx))
  for (j in seq_along(impute_idx)) {
    results[[j]] <- .impute_one_protein(
      model_df$data[[impute_idx[j]]],
      model_strategy,
      lod,
      response,
      sample_template,
      borrowed,
      df_method,
      nr_children_col
    )
  }

  # Bulk-assign to avoid per-element copy-on-modify (O(n²) → O(n))
  succeeded <- !vapply(results, is.null, logical(1))
  idx <- impute_idx[succeeded]
  results <- results[succeeded]

  if (length(idx) > 0) {
    model_df$linear_model[idx] <- lapply(results, `[[`, "linear_model")
    model_df$data[idx] <- lapply(results, `[[`, "data")
    model_df$has_model_fit[idx] <- TRUE
    model_df$isSingular[idx] <- FALSE
    model_df$sigma[idx] <- vapply(results, `[[`, numeric(1), "sigma")
    model_df$df.residual[idx] <- vapply(results, `[[`, numeric(1), "df.residual")
    model_df$nr_coef[idx] <- vapply(results, `[[`, numeric(1), "nr_coef")
    model_df$nr_coef_not_NA[idx] <- vapply(results, `[[`, numeric(1), "nr_coef_not_NA")
  }

  return(model_df)
}

# Impute and refit a single protein's model
#
# Takes one protein's nested data, imputes missing values with LOD,
# refits the model, and wraps it with borrowed covariance.
#
# @param dat data.frame for one protein (nested row)
# @param model_strategy strategy object with model_fun
# @param lod limit of detection value
# @param response response column name
# @param sample_template data.frame with all sample/group combinations
# @param borrowed list from compute_borrowed_variance()
# @param df_method "observed" or "borrowed"
# @return named list of updated fields, or NULL if refit fails
# @keywords internal
.impute_one_protein <- function(
  dat,
  model_strategy,
  lod,
  response,
  sample_template,
  borrowed,
  df_method,
  nr_children_col = NULL
) {
  n_observed <- sum(!is.na(dat[[response]]))

  # Complete data with all samples so missing groups get rows
  dat <- dplyr::left_join(sample_template, dat, by = intersect(colnames(sample_template), colnames(dat)))

  # Fill nr_children with 1 for imputed rows (no observed peptides in missing group)
  if (!is.null(nr_children_col) && nr_children_col %in% colnames(dat)) {
    dat[[nr_children_col]] <- ifelse(is.na(dat[[nr_children_col]]), 1, dat[[nr_children_col]])
  }

  # Impute NAs with LOD, clamp all values to max(value, LOD)
  dat[[response]] <- ifelse(is.na(dat[[response]]), lod, dat[[response]])
  dat[[response]] <- pmax(dat[[response]], lod)

  new_model <- model_strategy$model_fun(dat)
  if (is.character(new_model)) {
    return(NULL)
  }

  p <- length(stats::coefficients(new_model))

  # Degrees of freedom
  if (df_method == "observed") {
    imp_df <- max(n_observed - p, 1)
  } else {
    imp_df <- borrowed$df
  }

  # Borrowed covariance
  if (borrowed$method == "sigma") {
    cov_unscaled <- summary(new_model)$cov.unscaled
    imp_vcov <- borrowed$sigma^2 * cov_unscaled
  } else {
    imp_vcov <- borrowed$vcov
  }

  wrapped <- new_lm_imputed(
    new_model,
    borrowed_vcov = imp_vcov,
    borrowed_sigma = borrowed$sigma,
    borrowed_df = imp_df,
    n_observed = n_observed
  )

  list(
    linear_model = wrapped,
    data = dat,
    sigma = borrowed$sigma,
    df.residual = imp_df,
    nr_coef = p,
    nr_coef_not_NA = p
  )
}


# Fit the models to data ----

#' check if lm model is singular
#' @keywords internal
#' @family modelling
#' @export
#' @examples
#' fit <- stats::lm(Sepal.Length ~ Species, data = iris)
#' is_singular_lm(fit)
#'
is_singular_lm <- function(m) {
  has_na <- any(is.na(coefficients(m)))
  if (has_na) {
    return(TRUE)
  } else {
    if (df.residual(m) >= 2) {
      return(FALSE)
    }
    return(TRUE)
  }
}

#' retrieve complete models.
#' @keywords internal
#' @family modelling
#' @export
#' @examples
#' x <- sim_build_models_lmer(model = "factors", Nprot = 10)
#' cfits <- get_complete_model_fit(x$model_df)
#' stopifnot(nrow(cfits) == 6)
get_complete_model_fit <- function(complete_models) {
  complete_models <- complete_models |> dplyr::filter(.data$has_model_fit == TRUE)
  complete_models <- complete_models |>
    dplyr::filter(.data$nr_coef_not_NA == max(.data$nr_coef_not_NA)) |>
    dplyr::arrange(dplyr::desc(.data$nr_coef_not_NA))
  complete_models <- complete_models |> dplyr::filter(df.residual > 1)
  return(complete_models)
}

#' analyses lmer4 and lm models created using help function `strategy_lm` or `strategy_lmer`
#'
#' used in project p2901
#'
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#'
#' x <- sim_lfq_data_peptide_config()
#' formula_randomPeptide <-
#'   strategy_lmer("abundance  ~ group_ + (1 | peptide_Id)")
#' mr <- model_analyse( x$data,
#'  formula_randomPeptide,
#'  subject_id = x$config$hierarchy_keys_depth())
#' stopifnot(nrow(get_complete_model_fit(mr$model_df)) == 6)
#'
model_analyse <- function(
  pepIntensity,
  model_strategy,
  subject_id = "protein_Id",
  model_name = "Model"
) {
  nested_proteins <- pepIntensity |>
    dplyr::group_by(!!!syms(subject_id)) |>
    tidyr::nest()

  lmermodel <- "linear_model"

  pb <- progress::progress_bar$new(total = nrow(nested_proteins))
  model_proteins <- nested_proteins |>
    dplyr::mutate(!!lmermodel := purrr::map(data, model_strategy$model_fun, pb = pb))

  model_proteins <- model_proteins |>
    dplyr::mutate(
      !!"has_model_fit" := purrr::map_lgl(!!sym(lmermodel), function(x) {
        !is.character(x)
      })
    )

  count_coef <- function(x) {
    cc <- coefficients(x)
    if (inherits(cc, "numeric")) length(cc) else ncol(cc[[1]])
  }

  count_coef_not_NA <- function(x) {
    cc <- coefficients(x)
    if (inherits(cc, "numeric")) sum(!is.na(cc)) else ncol(cc[[1]])
  }

  model_proteins <- model_proteins |>
    dplyr::mutate(
      isSingular = purrr::map2_lgl(
        !!sym(lmermodel),
        .data$has_model_fit,
        function(m, ok) if (ok) model_strategy$isSingular(m) else NA
      ),
      df.residual = purrr::map2_dbl(
        !!sym(lmermodel),
        .data$has_model_fit,
        function(m, ok) if (ok) model_strategy$df_residual(m) else NA_real_
      ),
      sigma = purrr::map2_dbl(
        !!sym(lmermodel),
        .data$has_model_fit,
        function(m, ok) if (ok) model_strategy$sigma(m) else NA_real_
      ),
      nr_coef = purrr::map2_int(
        !!sym(lmermodel),
        .data$has_model_fit,
        function(m, ok) if (ok) count_coef(m) else NA_integer_
      ),
      nr_coef_not_NA = purrr::map2_int(
        !!sym(lmermodel),
        .data$has_model_fit,
        function(m, ok) if (ok) count_coef_not_NA(m) else NA_integer_
      )
    )

  return(list(model_df = model_proteins, model_name = model_name))
}


# visualize lmer modelling results ----

#' Plot prdictions
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#' m <- sim_make_model_lmer()
#' plot_lmer_peptide_predictions(m, intensity = "abundance")
#' m <- sim_make_model_lmer("interaction")
#' plot_lmer_peptide_predictions(m, intensity = "abundance")
plot_lmer_peptide_predictions <- function(m, intensity = "abundance") {
  data <- m@frame
  data$prediction <- predict(m)
  interaction_columns <- intersect(attributes(terms(m))$term.labels, colnames(data))
  data <- make_interaction_column(data, interaction_columns, sep = ":")
  gg <- ggplot(data, aes(x = .data$interaction, y = !!sym(intensity))) + geom_point()
  gg <- gg + geom_point(aes(x = .data$interaction, y = .data$prediction), color = 2) + facet_wrap(~peptide_Id)
  gg <- gg + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  return(gg)
}
