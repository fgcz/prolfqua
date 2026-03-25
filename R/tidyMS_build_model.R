#' Likelihood ratio test
#' @family modelling
#' @export
#' @param modelProteinF table with models (see build model)
#' @param modelName name of model
#' @param modelProteinF_Int reduced model
#' @param modelName_Int name of reduced model
#' @param subject_Id subject id typically Assession or protein_Id
#' @param path default NULL, set to a directory if you need to write diagnostic plots.
#' @examples
#' data_2Factor <- prolfqua::sim_lfq_data_2Factor_config(
#'  Nprot = 200,
#'  with_missing = TRUE,
#'  weight_missing = 2)
#'
#' pMerged <- LFQData$new(data_2Factor$data, data_2Factor$config)
#'
#' pMerged$config$get_response()
#' pMerged$factors()
#'
#' formula_condition_and_Batches <-
#'   prolfqua::strategy_lm("abundance ~ Treatment + Background")
#' modCB <- prolfqua::build_model(
#'   pMerged$data,
#'   formula_condition_and_Batches,
#'   subject_Id = pMerged$config$hierarchy_keys() )
#'
#' formula_condition <-
#'   prolfqua::strategy_lm("abundance ~ Treatment")
#' modC <- prolfqua::build_model(
#'   pMerged$data,
#'   formula_condition,
#'   subject_Id = pMerged$config$hierarchy_keys() )
#'
#' tmp <- LR_test(modCB$modelDF, "modCB", modC$modelDF, "modB")
#' hist(tmp$likelihood_ratio_test.pValue)
#'
LR_test <- function(
  modelProteinF,
  modelName,
  modelProteinF_Int,
  modelName_Int,
  subject_Id = "protein_Id",
  path = NULL
) {
  # Model Comparison
  reg <- dplyr::inner_join(
    dplyr::select(modelProteinF, !!sym(subject_Id), "linear_model"),
    dplyr::select(modelProteinF_Int, !!sym(subject_Id), "linear_model"),
    by = subject_Id
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
    dplyr::select(!!sym(subject_Id), dplyr::all_of("modelComparisonLikelihoodRatioTest")) |>
    tidyr::unnest(cols = c("modelComparisonLikelihoodRatioTest"))
  likelihood_ratio_test_result <- likelihood_ratio_test_result |>
    dplyr::rename(likelihood_ratio_test.pValue = .data$modelComparisonLikelihoodRatioTest)

  if (!is.null(path)) {
    plot_lrt_diagnostics(likelihood_ratio_test_result, modelName, modelName_Int, path)
  }

  return(likelihood_ratio_test_result)
}

#' Write LRT diagnostic plots to PDF
#'
#' Writes a histogram of p-values and an empirical CDF to a PDF file.
#' Called by \code{\link{LR_test}} when \code{path} is non-NULL.
#'
#' @param result tibble with a \code{likelihood_ratio_test.pValue} column
#' @param modelName name of the full model
#' @param modelName_Int name of the reduced/interaction model
#' @param path directory to write the PDF into
#' @keywords internal
plot_lrt_diagnostics <- function(result, modelName, modelName_Int, path) {
  fileName <- paste("hist_LRT_", modelName, "_", modelName_Int, ".pdf", sep = "")
  fileName <- file.path(path, fileName)
  message("writing figure : ", fileName, "\n")
  pdf(fileName)
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
#' @param subject_Id grouping variable
#' @param modelName model name
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
#' modelName <- "f_condtion_r_peptide"
#' formula_randomPeptide <-
#'   strategy_lmer("abundance  ~ group_ + (1 | peptide_Id) + (1 | sampleName)",
#'    model_name = modelName)
#'
#'
#' mod <- prolfqua::build_model(
#'  D$data,
#'  formula_randomPeptide,
#'  modelName = modelName,
#'  subject_Id = D$config$hierarchy_keys_depth())
#' aovtable <- mod$get_anova()
#'
#' mod <- prolfqua::build_model(
#'  LFQData$new(D$data, D$config),
#'  formula_randomPeptide,
#'  modelName = modelName)
#' model_summary(mod)
#'
#'
build_model <- function(
  data,
  model_strategy,
  subject_Id = if ("LFQData" %in% class(data)) {
    data$subject_Id()
  } else {
    "protein_Id"
  },
  modelName = model_strategy$model_name
) {
  dataX <- if ("LFQData" %in% class(data)) data$data else data
  modellingResult <- model_analyse(dataX, model_strategy, modelName = modelName, subject_Id = subject_Id)
  return(Model$new(
    modelDF = modellingResult$modelDF,
    model_strategy = model_strategy,
    modelName = modellingResult$modelName,
    subject_Id = subject_Id
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
#' @param modelName model name (default appends "Imputed")
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
#' strat <- strategy_lm(paste(lfqdata$config$get_response(), "~ group_"))
#' mod <- build_model_impute(lfqdata, strat)
#'
build_model_impute <- function(
  lfqdata,
  model_strategy,
  modelName = paste0(model_strategy$model_name, "Imputed"),
  lod = NULL,
  borrow_method = c("sigma", "vcov"),
  df_method = c("observed", "borrowed")
) {
  borrow_method <- match.arg(borrow_method)
  df_method <- match.arg(df_method)
  subject_Id <- lfqdata$subject_Id()
  response <- lfqdata$config$get_response()

  modellingResult <- model_analyse(
    lfqdata$data,
    model_strategy,
    modelName = modelName,
    subject_Id = subject_Id
  )

  if (is.null(lod)) {
    mh <- MissingHelpers$new(lfqdata$data, lfqdata$config)
    lod <- mh$get_LOD()
  }

  # Build sample template: all unique sample rows (excluding subject_Id and response)
  non_subject_cols <- setdiff(colnames(lfqdata$data), c(subject_Id, response))
  sample_template <- lfqdata$data |>
    dplyr::select(dplyr::all_of(non_subject_cols)) |>
    dplyr::distinct()

  modellingResult$modelDF <- impute_refit_singular(
    modellingResult$modelDF,
    model_strategy,
    lod = lod,
    response = response,
    sample_template = sample_template,
    borrow_method = borrow_method,
    df_method = df_method
  )

  return(Model$new(
    modelDF = modellingResult$modelDF,
    model_strategy = model_strategy,
    modelName = modellingResult$modelName,
    subject_Id = subject_Id
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
  res$exists <- table(mod$modelDF$has_model_fit)
  res$isSingular <- table(mod$modelDF$isSingular)
  return(res)
}


# Model fitting internals ----

.likelihood_ratio_test <- function(modelNO, model) {
  res <- tryCatch(anova(modelNO, model), error = function(x) NULL)
  if (!is.null(res)) {
    res <- suppressWarnings(broom::tidy(res))[2, "p.value"]
    return(as.numeric(res))
  } else {
    return(NA)
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
#' @param modelDF tibble from model_analyse
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
#'   mod$modelDF, method = "sigma")
#' stopifnot(borrowed_s$method == "sigma")
#' stopifnot(is.numeric(borrowed_s$sigma) && borrowed_s$sigma > 0)
#' stopifnot(is.numeric(borrowed_s$df) && borrowed_s$df > 0)
#'
#' # Vcov method: element-wise median vcov from donors.
#' # Falls back to sigma if donor models have different coefficient counts.
#' mod_no_missing <- sim_build_models_lm(model = "parallel3",
#'   Nprot = 10, with_missing = FALSE)
#' borrowed_v <- prolfqua:::compute_borrowed_variance(
#'   mod_no_missing$modelDF, method = "vcov")
#' stopifnot(borrowed_v$method == "vcov")
#' stopifnot(is.matrix(borrowed_v$vcov))
compute_borrowed_variance <- function(modelDF, method = c("sigma", "vcov")) {
  method <- match.arg(method)
  good <- get_complete_model_fit(modelDF)
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
#' @param modelDF tibble from model_analyse
#' @param model_strategy strategy list from strategy_lm etc.
#' @param lod numeric, limit of detection value
#' @param response character, response column name in nested data
#' @param sample_template data.frame with all sample/group combinations
#'   (columns matching the nested data minus the response). Used to complete
#'   proteins that are entirely missing in one or more groups.
#' @param borrow_method "sigma" or "vcov"
#' @param df_method "observed" uses max(n_observed - p, 1),
#'   "borrowed" uses median df from successful fits
#' @return modified modelDF with imputed models replacing failed/singular ones
#' @keywords internal
#' @family modelling
impute_refit_singular <- function(
  modelDF,
  model_strategy,
  lod,
  response,
  sample_template,
  borrow_method = c("sigma", "vcov"),
  df_method = c("observed", "borrowed")
) {
  borrow_method <- match.arg(borrow_method)
  df_method <- match.arg(df_method)

  max_coef <- max(modelDF$nr_coef, na.rm = TRUE)

  needs_impute <- (!modelDF$has_model_fit) |
    (!is.na(modelDF$isSingular) & modelDF$isSingular) |
    (!is.na(modelDF$nr_coef) & modelDF$nr_coef < max_coef)

  if (!any(needs_impute)) {
    return(modelDF)
  }

  borrowed <- compute_borrowed_variance(modelDF, method = borrow_method)

  impute_idx <- which(needs_impute)
  results <- vector("list", length(impute_idx))
  for (j in seq_along(impute_idx)) {
    results[[j]] <- .impute_one_protein(
      modelDF$data[[impute_idx[j]]],
      model_strategy,
      lod,
      response,
      sample_template,
      borrowed,
      df_method
    )
  }

  # Bulk-assign to avoid per-element copy-on-modify (O(n²) → O(n))
  succeeded <- !vapply(results, is.null, logical(1))
  idx <- impute_idx[succeeded]
  results <- results[succeeded]

  if (length(idx) > 0) {
    modelDF$linear_model[idx] <- lapply(results, `[[`, "linear_model")
    modelDF$data[idx] <- lapply(results, `[[`, "data")
    modelDF$has_model_fit[idx] <- TRUE
    modelDF$isSingular[idx] <- FALSE
    modelDF$sigma[idx] <- vapply(results, `[[`, numeric(1), "sigma")
    modelDF$df.residual[idx] <- vapply(results, `[[`, numeric(1), "df.residual")
    modelDF$nr_coef[idx] <- vapply(results, `[[`, numeric(1), "nr_coef")
    modelDF$nr_coef_not_NA[idx] <- vapply(results, `[[`, numeric(1), "nr_coef_not_NA")
  }

  return(modelDF)
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
.impute_one_protein <- function(dat, model_strategy, lod, response, sample_template, borrowed, df_method) {
  n_observed <- sum(!is.na(dat[[response]]))

  # Complete data with all samples so missing groups get rows
  dat <- dplyr::left_join(sample_template, dat, by = intersect(colnames(sample_template), colnames(dat)))

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
#'
isSingular_lm <- function(m) {
  anyNA <- any(is.na(coefficients(m)))
  if (anyNA) {
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
#' cfits <- get_complete_model_fit(x$modelDF)
#' stopifnot(nrow(cfits) == 6)
get_complete_model_fit <- function(modelProteinF) {
  modelProteinF <- modelProteinF |> dplyr::filter(.data$has_model_fit == TRUE)
  modelProteinF <- modelProteinF |>
    dplyr::filter(.data$nr_coef_not_NA == max(.data$nr_coef_not_NA)) |>
    dplyr::arrange(dplyr::desc(.data$nr_coef_not_NA))
  modelProteinF <- modelProteinF |> dplyr::filter(df.residual > 1)
  return(modelProteinF)
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
#'  subject_Id = x$config$hierarchy_keys_depth())
#' stopifnot(nrow(get_complete_model_fit(mr$modelDF)) == 6)
#'
model_analyse <- function(
  pepIntensity,
  model_strategy,
  subject_Id = "protein_Id",
  modelName = "Model"
) {
  nestProtein <- pepIntensity |>
    dplyr::group_by(!!!syms(subject_Id)) |>
    tidyr::nest()

  lmermodel <- "linear_model"

  pb <- progress::progress_bar$new(total = nrow(nestProtein))
  modelProtein <- nestProtein |>
    dplyr::mutate(!!lmermodel := purrr::map(data, model_strategy$model_fun, pb = pb))

  modelProtein <- modelProtein |>
    dplyr::mutate(
      !!"has_model_fit" := purrr::map_lgl(!!sym(lmermodel), function(x) {
        !is.character(x)
      })
    )

  modelProteinF <- modelProtein |>
    dplyr::filter(!!sym("has_model_fit") == TRUE)

  count_coef <- function(x) {
    cc <- coefficients(x)
    if (inherits(cc, "numeric")) {
      return(length(cc))
    } else {
      return(ncol(cc[[1]]))
    }
  }

  count_coef_not_NA <- function(x) {
    cc <- coefficients(x)
    if (inherits(cc, "numeric")) {
      return(sum(!is.na(cc)))
    } else {
      return(ncol(cc[[1]]))
    }
  }

  modelProteinF <- modelProteinF |>
    dplyr::mutate(
      isSingular = purrr::map_lgl(!!sym(lmermodel), model_strategy$isSingular),
      df.residual = purrr::map_dbl(!!sym(lmermodel), model_strategy$df_residual),
      sigma = purrr::map_dbl(!!sym(lmermodel), model_strategy$sigma),
      nr_coef = purrr::map_int(!!sym(lmermodel), count_coef),
      nr_coef_not_NA = purrr::map_int(!!sym(lmermodel), count_coef_not_NA)
    )

  modelProteinF <- modelProteinF |>
    dplyr::select(all_of(c(subject_Id, "isSingular", "df.residual", "sigma", "nr_coef", "nr_coef_not_NA")))
  modelProtein <- dplyr::left_join(modelProtein, modelProteinF)

  return(list(modelDF = modelProtein, modelName = modelName))
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
  interactionColumns <- intersect(attributes(terms(m))$term.labels, colnames(data))
  data <- make_interaction_column(data, interactionColumns, sep = ":")
  gg <- ggplot(data, aes(x = .data$interaction, y = !!sym(intensity))) + geom_point()
  gg <- gg + geom_point(aes(x = .data$interaction, y = .data$prediction), color = 2) + facet_wrap(~peptide_Id)
  gg <- gg + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  return(gg)
}
