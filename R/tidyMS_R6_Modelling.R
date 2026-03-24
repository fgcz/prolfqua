# Creating models from configuration ----

.ehandler <- function(e) {
  warning("WARN :", e)
  # return string here
  as.character(e)
}

#' Residual degrees of freedom for rlm objects
#'
#' \code{stats::df.residual} returns \code{NA} for \code{\link[MASS]{rlm}}
#' objects. This S3 method computes weighted residual df instead.
#'
#' @param object an \code{rlm} object
#' @param ... ignored
#' @return numeric scalar
#' @keywords internal
#' @export
df.residual.rlm <- function(object, ...) {
  sum(object$w) - object$rank
}

#' Residual scale estimate for rlm objects
#'
#' \code{stats::sigma} returns \code{NA} for \code{\link[MASS]{rlm}} objects.
#' This S3 method computes a weighted residual scale estimate instead.
#'
#' @param object an \code{rlm} object
#' @param ... ignored
#' @return numeric scalar
#' @keywords internal
#' @export
sigma.rlm <- function(object, ...) {
  sqrt(sqrt(sum(object$w * object$resid^2) / (sum(object$w) - object$rank)))
}

#' Linear mixed-effects model strategy (R6 class)
#'
#' Encapsulates everything needed to fit per-protein linear mixed-effects models
#' via \code{\link[lmerTest]{lmer}} and extract contrasts.
#'
#' @export
#' @family modelling
#' @examples
#' istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 10, with_missing = FALSE)
#' istar <- prolfqua::LFQData$new(istar$data, istar$config)
#' istar$data <- istar$data |> dplyr::group_by(protein_Id) |>
#'   dplyr::mutate(abundanceC = abundance - mean(abundance)) |> dplyr::ungroup()
#' strat <- StrategyLmer$new("abundanceC ~ group_ + (1|peptide_Id)",
#'   model_name = "random_example")
#' mod <- build_model(istar, strat)
#' sum(mod$modelDF$exists_lmer)
StrategyLmer <- R6::R6Class(
  "StrategyLmer",
  public = list(
    #' @field formula model formula
    formula = NULL,
    #' @field model_name name of model
    model_name = NULL,
    #' @field report_columns columns to report
    report_columns = NULL,
    #' @field is_mixed always TRUE for lmer
    is_mixed = TRUE,
    #' @field anova_df list with anova function and column names
    anova_df = NULL,

    #' @description Create a new StrategyLmer
    #' @param modelstr model formula string
    #' @param model_name name of model
    #' @param report_columns columns to report
    initialize = function(
      modelstr,
      model_name = "Model",
      report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value", "moderated.p.value.adjusted")
    ) {
      self$formula <- as.formula(modelstr)
      self$model_name <- model_name
      self$report_columns <- report_columns
      self$anova_df <- get_anova_df(test = "F")
    },

    #' @description Fit lmer to one protein's data
    #' @param x data.frame for one protein
    #' @param pb optional progress bar
    #' @param get_formula if TRUE, return formula instead of fitting
    model_fun = function(x, pb, get_formula = FALSE) {
      if (get_formula) {
        return(self$formula)
      }
      if (!missing(pb)) {
        pb$tick()
      }
      tryCatch(lmerTest::lmer(self$formula, data = x), error = .ehandler)
    },

    #' @description Check if model is singular
    #' @param model fitted model
    isSingular = function(model) lme4::isSingular(model),

    #' @description Compute contrasts from fitted model
    #' @param ... passed to \code{\link{my_contest}}
    contrast_fun = function(...) my_contest(...),

    #' @description Get residual degrees of freedom
    #' @param model fitted model
    df_residual = function(model) stats::df.residual(model),

    #' @description Get residual standard error
    #' @param model fitted model
    sigma = function(model) stats::sigma(model)
  )
)


#' Create linear mixed-effects model strategy
#'
#' Convenience wrapper that creates a \code{\link{StrategyLmer}} object.
#' @rdname strategy
#' @export
#' @param modelstr model formula
#' @param model_name name of model
#' @param report_columns columns to report
#' @family modelling
#' @return a \code{\link{StrategyLmer}} object
#' @examples
#' modelFunction <- strategy_lmer("abundanceC ~ group_ + (1|peptide_Id)",
#'   model_name = "random_example")
#' modelFunction$model_fun(get_formula = TRUE)
strategy_lmer <- function(
  modelstr,
  model_name = "Model",
  report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value", "moderated.p.value.adjusted")
) {
  StrategyLmer$new(modelstr, model_name, report_columns)
}

#' Linear model strategy (R6 class)
#'
#' Encapsulates everything needed to fit per-protein linear models and extract
#' contrasts: the formula, model fitting function, singularity check, contrast
#' computation, ANOVA, and residual statistics.
#'
#' @export
#' @family modelling
#' @examples
#' strat <- StrategyLM$new("Intensity ~ condition", model_name = "parallel design")
#' strat$model_fun(get_formula = TRUE)
#' strat$weights
StrategyLM <- R6::R6Class(
  "StrategyLM",
  public = list(
    #' @field formula model formula
    formula = NULL,
    #' @field model_name name of model
    model_name = NULL,
    #' @field report_columns columns to report
    report_columns = NULL,
    #' @field weights optional character string naming a column in the data
    #'   containing per-observation weights, passed to \code{\link[stats]{lm}}.
    weights = NULL,
    #' @field is_mixed always FALSE for lm
    is_mixed = FALSE,
    #' @field anova_df list with anova function and column names
    anova_df = NULL,

    #' @description Create a new StrategyLM
    #' @param modelstr model formula string
    #' @param model_name name of model
    #' @param report_columns columns to report
    #' @param weights optional character string naming a column in the data
    #'   containing per-observation weights
    initialize = function(
      modelstr,
      model_name = "Model",
      report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value", "moderated.p.value.adjusted"),
      weights = NULL
    ) {
      self$formula <- as.formula(modelstr)
      self$model_name <- model_name
      self$report_columns <- report_columns
      self$weights <- weights
      self$anova_df <- get_anova_df(test = "F")
    },

    #' @description Fit lm to one protein's data
    #' @param x data.frame for one protein
    #' @param pb optional progress bar
    #' @param get_formula if TRUE, return formula instead of fitting
    model_fun = function(x, pb, get_formula = FALSE) {
      if (get_formula) {
        return(self$formula)
      }
      if (!missing(pb)) {
        pb$tick()
      }
      if (!is.null(self$weights)) {
        call <- bquote(lm(.(self$formula), data = x, weights = .(as.name(self$weights))))
        tryCatch(eval(call), error = .ehandler)
      } else {
        tryCatch(lm(self$formula, data = x), error = .ehandler)
      }
    },

    #' @description Check if model is singular
    #' @param model fitted model
    isSingular = function(model) isSingular_lm(model),

    #' @description Compute contrasts from fitted model
    #' @param ... passed to \code{\link{my_contrast_V2}}
    contrast_fun = function(...) my_contrast_V2(...),

    #' @description Get residual degrees of freedom
    #' @param model fitted model
    df_residual = function(model) stats::df.residual(model),

    #' @description Get residual standard error
    #' @param model fitted model
    sigma = function(model) stats::sigma(model)
  )
)


#' Create linear model strategy
#'
#' Convenience wrapper that creates a \code{\link{StrategyLM}} object.
#' @rdname strategy
#' @export
#' @param modelstr model formula
#' @param model_name name of model
#' @param report_columns columns to report
#' @param weights optional character string naming a column in the data
#'   containing per-observation weights, passed to \code{\link[stats]{lm}}.
#'   Default \code{NULL} (unweighted).
#' @family modelling
#' @return a \code{\link{StrategyLM}} object
#' @examples
#' tmp <- strategy_lm("Intensity ~ condition", model_name = "parallel design")
#' tmp$model_fun(get_formula = TRUE)
#' tmp$weights
strategy_lm <- function(
  modelstr,
  model_name = "Model",
  report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value", "moderated.p.value.adjusted"),
  weights = NULL
) {
  StrategyLM$new(modelstr, model_name, report_columns, weights)
}


#' Robust linear model strategy (R6 class)
#'
#' Encapsulates everything needed to fit per-protein robust linear models
#' via \code{\link[MASS]{rlm}} and extract contrasts.
#'
#' @export
#' @family modelling
#' @examples
#' strat <- StrategyRLM$new("Intensity ~ condition", model_name = "parallel design")
#' strat$model_fun(get_formula = TRUE)
StrategyRLM <- R6::R6Class(
  "StrategyRLM",
  public = list(
    #' @field formula model formula
    formula = NULL,
    #' @field model_name name of model
    model_name = NULL,
    #' @field report_columns columns to report
    report_columns = NULL,
    #' @field is_mixed always FALSE for rlm
    is_mixed = FALSE,
    #' @field anova_df list with anova function and column names
    anova_df = NULL,

    #' @description Create a new StrategyRLM
    #' @param modelstr model formula string
    #' @param model_name name of model
    #' @param report_columns columns to report
    initialize = function(
      modelstr,
      model_name = "Model",
      report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value", "moderated.p.value.adjusted")
    ) {
      self$formula <- as.formula(modelstr)
      self$model_name <- model_name
      self$report_columns <- report_columns
      self$anova_df <- get_anova_df(test = "F")
    },

    #' @description Fit rlm to one protein's data
    #' @param x data.frame for one protein
    #' @param pb optional progress bar
    #' @param get_formula if TRUE, return formula instead of fitting
    model_fun = function(x, pb, get_formula = FALSE) {
      if (get_formula) {
        return(self$formula)
      }
      if (!missing(pb)) {
        pb$tick()
      }
      tryCatch(MASS::rlm(self$formula, data = x, method = "M"), error = .ehandler)
    },

    #' @description Check if model is singular (NA coefficients or df < 2)
    #' @param model fitted model
    isSingular = function(model) {
      if (any(is.na(coefficients(model)))) {
        return(TRUE)
      }
      df <- df.residual(model)
      if (is.na(df) || df < 2) {
        return(TRUE)
      }
      FALSE
    },

    #' @description Compute contrasts from fitted model
    #' @param ... passed to \code{\link{my_contrast_V2}}
    contrast_fun = function(...) my_contrast_V2(...),

    #' @description Get residual degrees of freedom
    #' @param model fitted model
    df_residual = function(model) stats::df.residual(model),

    #' @description Get residual standard error
    #' @param model fitted model
    sigma = function(model) stats::sigma(model)
  )
)


#' Create robust linear model strategy
#'
#' Convenience wrapper that creates a \code{\link{StrategyRLM}} object.
#' @rdname strategy
#' @export
#' @param modelstr model formula
#' @param model_name name of model
#' @param report_columns columns to report
#' @family modelling
#' @return a \code{\link{StrategyRLM}} object
#' @examples
#' tmp <- strategy_rlm("Intensity ~ condition", model_name = "parallel design")
#' tmp$model_fun(get_formula = TRUE)
strategy_rlm <- function(
  modelstr,
  model_name = "Model",
  report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value", "moderated.p.value.adjusted")
) {
  StrategyRLM$new(modelstr, model_name, report_columns)
}


#' R6 class for extracting ANOVA results as a data frame
#'
#' Wraps the ANOVA extraction logic and associated column names.
#'
#' @export
#' @family modelling
AnovaExtractor <- R6::R6Class(
  "AnovaExtractor",
  public = list(
    #' @field test statistical test type ("F", "Chisq", etc.)
    test = "F",
    #' @field col_pval p-value column name in ANOVA output
    col_pval = character(),
    #' @field col_fdr FDR column name
    col_fdr = character(),
    #' @description Create a new AnovaExtractor
    #' @param test statistical test type
    initialize = function(test = "F") {
      self$test <- test
      self$col_pval <- paste0("Pr..", substr(test, 1, 3), ".")
      self$col_fdr <- paste0("FDR.Pr..", substr(test, 1, 3), ".")
    },
    #' @description Extract ANOVA table from a fitted model
    #' @param x a fitted model
    #' @return data.frame with factor column and ANOVA statistics
    fun = function(x) {
      x <- anova(x, test = self$test)
      colnames(x) <- make.names(colnames(x))
      x <- data.frame(factor = rownames(x), x)
      return(x)
    }
  )
)

#' anova returning dataframe
#' @keywords internal
#' @family modelling
#' @export
#' @examples
#' x <- get_anova_df(test = "F")
#' x <- get_anova_df(test = "Chisq")
get_anova_df <- function(test = "F") {
  AnovaExtractor$new(test = test)
}


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

  max_coef <- max(modelDF$nrcoef, na.rm = TRUE)

  needs_impute <- (!modelDF$exists_lmer) |
    (!is.na(modelDF$isSingular) & modelDF$isSingular) |
    (!is.na(modelDF$nrcoef) & modelDF$nrcoef < max_coef)

  if (!any(needs_impute)) {
    return(modelDF)
  }

  borrowed <- compute_borrowed_variance(modelDF, method = borrow_method)

  for (i in which(needs_impute)) {
    dat <- modelDF$data[[i]]
    n_observed <- sum(!is.na(dat[[response]]))

    # Complete data with all samples so missing groups get rows
    dat <- dplyr::left_join(sample_template, dat, by = intersect(colnames(sample_template), colnames(dat)))

    # Impute NAs with LOD, clamp all values to max(value, LOD)
    dat[[response]] <- ifelse(is.na(dat[[response]]), lod, dat[[response]])
    dat[[response]] <- pmax(dat[[response]], lod)

    new_model <- model_strategy$model_fun(dat)
    if (is.character(new_model)) {
      next
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

    modelDF$linear_model[[i]] <- wrapped
    modelDF$data[[i]] <- dat
    modelDF$exists_lmer[[i]] <- TRUE
    modelDF$isSingular[[i]] <- FALSE
    modelDF$sigma[[i]] <- borrowed$sigma
    modelDF$df.residual[[i]] <- imp_df
    modelDF$nrcoef[[i]] <- p
    modelDF$nrcoeff_not_NA[[i]] <- p
  }

  return(modelDF)
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
  modelProteinF <- modelProteinF |> dplyr::filter(.data$exists_lmer == TRUE)
  modelProteinF <- modelProteinF |>
    dplyr::filter(.data$nrcoeff_not_NA == max(.data$nrcoeff_not_NA)) |>
    dplyr::arrange(dplyr::desc(.data$nrcoeff_not_NA))
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
      !!"exists_lmer" := purrr::map_lgl(!!sym(lmermodel), function(x) {
        !is.character(x)
      })
    )

  modelProteinF <- modelProtein |>
    dplyr::filter(!!sym("exists_lmer") == TRUE)
  modelProteinF <- modelProteinF |>
    dplyr::mutate(!!"isSingular" := purrr::map_lgl(!!sym(lmermodel), model_strategy$isSingular))
  modelProteinF <- modelProteinF |>
    dplyr::mutate(!!"df.residual" := purrr::map_dbl(!!sym(lmermodel), model_strategy$df_residual))
  modelProteinF <- modelProteinF |>
    dplyr::mutate(!!"sigma" := purrr::map_dbl(!!sym(lmermodel), model_strategy$sigma))

  nrcoeff <- function(x) {
    cc <- coefficients(x)
    if (inherits(cc, "numeric")) {
      return(length(cc))
    } else {
      return(ncol(cc[[1]]))
    }
  }

  nrcoeff_not_NA <- function(x) {
    cc <- coefficients(x)
    if (inherits(cc, "numeric")) {
      return(sum(!is.na(cc)))
    } else {
      return(ncol(cc[[1]]))
    }
  }

  modelProteinF <- modelProteinF |> dplyr::mutate(nrcoef = purrr::map_int(!!sym(lmermodel), nrcoeff))
  modelProteinF <- modelProteinF |> dplyr::mutate(nrcoeff_not_NA = purrr::map_int(!!sym(lmermodel), nrcoeff_not_NA))

  modelProteinF <- modelProteinF |>
    dplyr::select(all_of(c(subject_Id, "isSingular", "df.residual", "sigma", "nrcoef", "nrcoeff_not_NA")))
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

# Generate linear functions -----

.model_coeff_matrix <- function(m) {
  data <- NULL
  if ("lmerModLmerTest" %in% class(m)) {
    data <- m@frame
    coeffs <- coefficients(summary(m))[, "Estimate"]
  } else {
    # for "lmerModLmerTest"
    data <- m$model
    coeffs <- coef(m)
  }
  interactionColumns <- intersect(attributes(terms(m))$term.labels, colnames(data))
  data <- make_interaction_column(data, interactionColumns, sep = ":")

  inter <- unique(data$interaction)
  mm <- matrix(0, nrow = length(inter), ncol = length(coeffs))
  rownames(mm) <- inter
  colnames(mm) <- names(coeffs)
  mm[, 1] <- 1
  coefi <- coeffs[-1]
  for (i in seq_along(coefi)) {
    # the grep is needed to extract coefficients of interaction terms belonging to a factor
    # I am using wor boundaries "\\b" to allow for factor levels that are substrings.
    positionIDX <- grep(paste0("\\b", names(coefi)[i], "\\b"), inter)
    mm[positionIDX, i + 1] <- 1
  }
  return(list(mm = mm, coeffs = coeffs))
}


.get_match_idx <- function(mm, factor_level) {
  ddd <- names_to_matrix(rownames(mm), split = ":")
  xd <- apply(
    ddd,
    2,
    function(x, factor_level) {
      x %in% factor_level
    },
    factor_level
  )
  idx <- which(apply(xd, 1, sum) > 0)
  return(idx)
}

.coeff_weights_factor_levels <- function(mm) {
  getCoeffs <- function(factor_level, mm) {
    idx <- .get_match_idx(mm, factor_level)
    x <- as.list(apply(mm[idx, , drop = FALSE], 2, mean))
    x <- tibble::as_tibble(x)
    tibble::add_column(x, "factor_level" = factor_level, .before = 1)
  }
  factor_levels <- unique(unlist(stringi::stri_split_fixed(rownames(mm), ":")))
  xx <- purrr::map_df(factor_levels, getCoeffs, mm)
  return(xx)
}

#' get linfct from model
#' @param m linear model
#' @family modelling
#' @export
#' @keywords internal
#' @examples
#'
#' m <- sim_make_model_lm("factors")
#' linfct <- linfct_from_model(m, as_list = TRUE)
#'
#' linfct$linfct_factors
#' linfct$linfct_interactions
#' lf <- matrix(
#' c(1, 1, 1, 1, 0.5, 0.5, 0, 1, 0, 1, 0.5, 0.5),
#' nrow = 4,
#' byrow = FALSE,
#' dimnames = list(c("BackgroundX", "BackgroundZ", "TreatmentA", "TreatmentB"),
#'                 c("(Intercept)", "TreatmentB", "BackgroundZ"))
#' )
#' stopifnot(lf == linfct$linfct_factors)
#' m <- sim_make_model_lm("interaction")
#' linfct <- linfct_from_model(m)
#'
#' m <- lm(Petal.Width ~ Species, data = iris)
#' linfct_from_model(m)
#' xx <- data.frame( Y = 1:10 , Condition = c(rep("a",5), rep("b",5)) )
#' m <- lm(Y ~ Condition, data = xx)
#' linfct_from_model(m)
#' xx <- data.frame( Y = 1:10 , Condition = c(rep("a",5), rep("b.b",5)) )
#' m <- lm(Y ~ Condition, data = xx)
#' linfct_from_model(m)
#' xx <- data.frame( Y = 1:10 , Condition = c(rep("a",5), rep("ab",5)) )
#' m <- lm(Y ~ Condition, data = xx)
#' linfct_from_model(m)
#'
linfct_from_model <- function(m, as_list = TRUE) {
  cm <- .model_coeff_matrix(m)
  cm_mm <- cm$mm[order(rownames(cm$mm)), ]

  l_factors <- .coeff_weights_factor_levels(cm_mm)
  linfct_factors <- l_factors |>
    dplyr::select(-dplyr::all_of("factor_level")) |>
    data.matrix()

  rownames(linfct_factors) <- l_factors$factor_level
  linfct_factors <- linfct_factors[order(rownames(linfct_factors)), ]
  res <- list(linfct_factors = linfct_factors, linfct_interactions = cm_mm)

  if (as_list) {
    return(res)
  } else {
    do.call(rbind, res)
  }
}


#' linfct_matrix_contrasts
#' @export
#' @param linfct linear functions as created by linfct_from_model
#' @param contrasts named character vector of contrasts to determine linear functions for
#' @param p.message print messages default FALSE
#'
#' @family modelling
#' @keywords internal
#' @examples
#'
#' m <- sim_make_model_lm( "factors")
#' Contr <- c("TreatmentA_vs_B" = "TreatmentA - TreatmentB",
#'     "BackgroundX_vs_Z" = "BackgroundX - BackgroundZ",
#'     "IntoflintoA" = "`TreatmentA:BackgroundX` - `TreatmentA:BackgroundZ`",
#'     "IntoflintoB" = "`TreatmentB:BackgroundX` - `TreatmentB:BackgroundZ`",
#'     "IntoflintoX" = "`TreatmentA:BackgroundX` - `TreatmentB:BackgroundX`",
#'     "IntoflintoZ" = "`TreatmentA:BackgroundZ` - `TreatmentB:BackgroundZ`",
#'     "interactXZ" = "IntoflintoX - IntoflintoZ",
#'     "interactAB" = "IntoflintoA - IntoflintoB"
#'      )
#' linfct <- linfct_from_model(m, as_list = FALSE)
#' x <- linfct_matrix_contrasts(linfct, Contr )
#' stopifnot(sum(x["interactXZ",]) == 0 )
#' stopifnot(sum(x["interactAB",]) == 0 )
#'
#' m <- sim_make_model_lm( "interaction")
#' linfct <- linfct_from_model(m, as_list = FALSE)
#' x <- linfct_matrix_contrasts(linfct, Contr )
#' stopifnot(sum(x["interactXZ",]) ==1 )
#' stopifnot(sum(x["interactAB",]) ==1 )
#'
linfct_matrix_contrasts <- function(linfct, contrasts, p.message = FALSE) {
  linfct <- t(linfct)
  df <- tibble::as_tibble(linfct, rownames = "interaction")
  make_contrasts <- function(data, contrasts) {
    cnams <- base::setdiff(colnames(data), "interaction")
    failures <- list()
    for (i in seq_along(contrasts)) {
      if (p.message) {
        message(names(contrasts)[i], "=", contrasts[i], "\n")
      }
      contrast_name <- names(contrasts)[i]
      if (is.null(contrast_name) || !nzchar(contrast_name)) {
        contrast_name <- paste0("contrast_", i)
      }
      err <- tryCatch(
        {
          data <- dplyr::mutate(data, !!names(contrasts)[i] := !!rlang::parse_expr(contrasts[i]))
          NULL
        },
        error = function(e) {
          e
        }
      )
      if (inherits(err, "error")) {
        failures[[length(failures) + 1]] <- list(
          contrast = contrast_name,
          message = conditionMessage(err)
        )
      }
    }
    res <- data |> dplyr::select(-all_of(cnams))
    if (length(failures) > 0) {
      failure_df <- dplyr::bind_rows(failures)
      failure_names <- paste(failure_df$contrast, collapse = ", ")
      failure_messages <- unique(failure_df$message)
      failure_summary <- paste(utils::head(failure_messages, 3), collapse = "; ")
      warning(
        paste0(
          "linfct_matrix_contrasts: computed ",
          ncol(res) - 1,
          "/",
          length(contrasts),
          " contrasts; failed ",
          nrow(failure_df),
          ": ",
          failure_names,
          ". ",
          failure_summary
        ),
        call. = FALSE
      )
    }
    return(res)
  }

  res <- make_contrasts(df, contrasts)
  res <- tibble::column_to_rownames(res, "interaction")
  res <- t(res)
  return(res)
}


#' create all possible contrasts
#' @export
#' @keywords internal
#' @family modelling
#' @examples
#' m <- sim_make_model_lm( "interaction")
#' linfct <- linfct_from_model(m)
#' xl <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_factors)
#' xx <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_interactions)
#' m <- sim_make_model_lm( "factor")
#' linfct <- linfct_from_model(m)
#' xl <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_factors)
#' xx <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_interactions)
#' m <- sim_make_model_lm( "parallel2")
#' linfct <- linfct_from_model(m)
#' xl <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_factors)
#' stopifnot(all(xl == c(0,-1)))
#' xx <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_interactions)
#' stopifnot(all(xx == c(0,-1)))
#'
#' m <- sim_make_model_lm( "parallel3")
#' linfct <- linfct_from_model(m)
#' xl <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_factors)
#' stopifnot(all(xl == c(0,0,0,-1,0,1,0,-1,-1)))
#' xx <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_interactions)
#' stopifnot(all(xl == c(0,0,0,-1,0,1,0,-1,-1)))
linfct_all_possible_contrasts <- function(lin_int) {
  combs <- combn(nrow(lin_int), 2)
  names <- rownames(lin_int)
  newnames <- rep("", ncol(combs))
  new_lin_fct <- matrix(NA, nrow = ncol(combs), ncol = ncol(lin_int))
  for (i in seq_len(ncol(combs))) {
    newnames[i] <- paste(names[combs[, i]], collapse = " - ")
    new_lin_fct[i, ] <- lin_int[combs[1, i], ] - lin_int[combs[2, i], ]
  }
  rownames(new_lin_fct) <- newnames
  colnames(new_lin_fct) <- colnames(lin_int)
  return(new_lin_fct)
}
#' create contrasts between factor levels
#'
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#'
#' m <- sim_make_model_lm( "interaction")
#' xl <- linfct_factors_contrasts(m)
#' m <- lm(Petal.Width ~ Species, data = iris)
#' linfct_factors_contrasts(m)
linfct_factors_contrasts <- function(m) {
  ffac <- attributes(terms(m))$term.labels
  ffac <- ffac[!grepl(":", ffac)] # remove interactions
  linfct_factors <- linfct_from_model(m)$linfct_factors

  factorDepths <- rownames(linfct_factors)
  res <- vector(length(ffac), mode = "list")
  for (i in seq_along(ffac)) {
    fac <- ffac[i]
    idx <- grep(fac, factorDepths)
    linfct_m <- linfct_factors[idx, ]
    res[[i]] <- linfct_all_possible_contrasts(linfct_m)
  }
  res <- do.call(rbind, res)
  return(res)
}

# Computing contrasts helpers -----

#' compute contrasts for full models
#' @param m linear model generated using lm
#' @param linfct linear function
#' @param strategy optional strategy for df and sigma computation
#' @param coef coefficients vector, default from model
#' @param Sigma.hat variance-covariance matrix, default from model
#' @param confint which confidence interval to determine
#'
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#'
#' m <-  sim_make_model_lm( "factors")
#' linfct <- linfct_from_model(m)$linfct_factors
#' my_contrast(m, linfct, confint = 0.95)
#' my_contrast(m, linfct, confint = 0.99)
#'
my_contrast <- function(m, linfct, strategy = NULL, coef = coefficients(m), Sigma.hat = vcov(m), confint = 0.95) {
  if (is.null(strategy)) {
    df <- df.residual(m)
    sigma <- sigma(m)
  } else {
    df <- strategy$df_residual(m)
    sigma <- strategy$sigma(m)
  }

  estimate <- linfct %*% t(t(coef))

  if (df > 0) {
    std.error <- sqrt(diag(linfct %*% Sigma.hat %*% t(linfct)))
    statistic <- estimate / std.error

    p.value <- pt(abs(statistic), df = df, lower.tail = FALSE) * 2
    prqt <- -qt((1 - confint) / 2, df = df)
    conf.low <- estimate - prqt * std.error
    conf.high <- estimate + prqt * std.error
  } else {
    std.error <- NA
    statistic <- NA
    p.value <- NA
    conf.low <- NA
    conf.high <- NA
  }

  res <- data.frame(
    lhs = rownames(linfct),
    sigma = sigma,
    df = df,
    estimate = estimate,
    std.error = std.error,
    statistic = statistic,
    p.value = p.value,
    conf.low = conf.low,
    conf.high = conf.high,
    stringsAsFactors = FALSE
  )
  return(res)
}

#' handles incomplete models
#'
#' only keeps non NA coefficients.
#'
#' @param m linear model generated using lm
#' @param linfct linear function
#' @param confint confidence interval default 0.95
#'
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#' m <- sim_make_model_lm( "factors")
#' linfct <- linfct_from_model(m)$linfct_factors
#' my_contrast_V2(m, linfct, confint = 0.95)
#' my_contrast_V2(m, linfct, confint = 0.99)
#'
my_contrast_V2 <- function(m, linfct, confint = 0.95) {
  Sigma.hat <- vcov(m)

  coef <- na.omit(coefficients(m))

  res <- vector(nrow(linfct), mode = "list")
  for (i in seq_len(nrow(linfct))) {
    linfct_v <- linfct[i, , drop = FALSE]
    idx <- which(linfct_v != 0)
    nam <- colnames(linfct_v)[idx]

    if (all(nam %in% names(coef))) {
      linfct_v_red <- linfct_v[, nam, drop = FALSE]
      Sigma.hat_red <- Sigma.hat[nam, nam, drop = FALSE]
      coef_red <- coef[nam]
      stopifnot(all.equal(colnames(linfct_v_red), colnames(Sigma.hat_red)))
      stopifnot(all.equal(colnames(linfct_v_red), names(coef_red)))
      res[[i]] <- my_contrast(m, linfct_v_red, coef = coef_red, Sigma.hat = Sigma.hat_red, confint = confint)
    } else {
      res[[i]] <- data.frame(
        lhs = rownames(linfct_v),
        sigma = sigma(m),
        df = df.residual(m),
        estimate = NA,
        std.error = NA,
        statistic = NA,
        p.value = NA,
        conf.low = NA,
        conf.high = NA,
        stringsAsFactors = FALSE
      )
    }
  }
  res <- dplyr::bind_rows(res)
  return(res)
}

#' applies contrast computation using lmerTest::contest function
#' @param model mixed effects model
#' @param linfct linear function
#' @param ddf method to determine denominator degrees of freedom
#' @family modelling
#' @export
#' @keywords internal
#' @examples
#'
#' mb <- sim_make_model_lmer("interaction")
#' summary(mb)
#'
#' linfct <- linfct_from_model(mb)
#' names(linfct)
#' my_contest(mb, linfct$linfct_factors)
#' my_contest(mb, linfct$linfct_interactions)
#' length(mb@beta)
#' lmerTest::contest(mb, c( 0 ,1 , 0 , 0),joint = FALSE)
#' summary(mb)
#'
my_contest <- function(model, linfct, ddf = c("Satterthwaite", "Kenward-Roger")) {
  ddf <- match.arg(ddf)
  if (length(lme4::fixef(model)) != ncol(linfct)) {
    warning("Model is rank deficient!")
    return(NA) # catch rank defficient
  } else {
    res <- lmerTest::contest(model, linfct, joint = FALSE, confint = TRUE, ddf = ddf)
  }
  res <- tibble::as_tibble(res, rownames = "lhs")
  res$sigma <- sigma(model)
  res <- res |>
    dplyr::rename(
      estimate = "Estimate",
      std.error = "Std. Error",
      statistic = "t value",
      p.value = "Pr(>|t|)",
      conf.low = "lower",
      conf.high = "upper"
    )
  return(res)
}

# computing contrast ----
#' pivot model contrasts matrix to wide format produced by `contrasts_linfct` and ...
#' @export
#' @family modelling
#' @param modelWithInteractionsContrasts data.frame with contrast results in long format
#' @param subject_Id column name(s) identifying subjects (e.g. protein_Id)
#' @param columns character vector of value columns to pivot wide
#' @param contrast column name containing contrast labels
#' @examples
#'
#' # this function is used by the contrast classes to implement the to wide method
#'
pivot_model_contrasts_2_Wide <- function(
  modelWithInteractionsContrasts,
  subject_Id = "protein_Id",
  columns = c("estimate", "p.value", "p.value.adjusted"),
  contrast = "lhs"
) {
  m_spread <- function(longContrasts, subject_Id, column, contrast) {
    res <- longContrasts |>
      dplyr::select(all_of(c(subject_Id, contrast, column)))
    res <- res |> dplyr::mutate(!!contrast := paste0(column, ".", !!sym(contrast)))
    res <- res |>
      tidyr::pivot_wider(
        names_from = dplyr::all_of(contrast),
        values_from = dplyr::all_of(column)
      )
    return(res)
  }
  res <- list()
  for (column in columns) {
    res[[column]] <- m_spread(modelWithInteractionsContrasts, subject_Id, column, contrast)
  }
  res <- res |> reduce(left_join, by = c(subject_Id))
  return(res)
}
#' compute group averages
#'
#' used in p2621, p2109
#'
#' @export
#' @keywords internal
#' @examples
#' modelSummary_A <- sim_build_models_lm()
#' m <- get_complete_model_fit(modelSummary_A$modelDF)
#'
#' factor_contrasts <- linfct_factors_contrasts( m$linear_model[[1]])
#'
#' factor_levelContrasts <- contrasts_linfct( m,
#'         factor_contrasts,
#'         subject_Id = "protein_Id",
#'         contrastfun = prolfqua::my_contrast_V2)
contrasts_linfct <- function(models, linfct, subject_Id = "protein_Id", contrastfun = prolfqua::my_contest) {
  message("contrasts_linfct")
  modelcol <- "linear_model"

  interaction_models <- vector(mode = "list", length = nrow(models))

  if ("matrix" %in% class(linfct)) {
    pb <- progress::progress_bar$new(total = length(models[[modelcol]]))
    for (i in seq_along(models[[modelcol]])) {
      interaction_models[[i]] <- contrastfun(models[[modelcol]][[i]], linfct = linfct)

      pb$tick()
    }
    interaction_model_matrix <- models
    interaction_model_matrix$contrast <- interaction_models
  } else if (("list" %in% class(linfct)) && (length(linfct) == nrow(models))) {
    pb <- progress::progress_bar$new(total = length(models[[modelcol]]))

    for (i in seq_along(models[[modelcol]])) {
      interaction_models[[i]] <- contrastfun(models[[modelcol]][[i]], linfct = linfct[[i]])
      pb$tick()
    }
    interaction_model_matrix <- models
    interaction_model_matrix$contrast <- interaction_models
  } else {
    stop("linct must be either a matrix or a list of length == nrow models")
  }

  mclass <- function(x) {
    class(x)[1]
  }

  interaction_model_matrix <- interaction_model_matrix |>
    dplyr::mutate(classC = purrr::map_chr(.data$contrast, mclass))

  failed_mask <- interaction_model_matrix$classC == "logical"
  n_failed <- sum(failed_mask)
  if (n_failed > 0) {
    failed_ids <- interaction_model_matrix[[subject_Id[1]]][failed_mask]
    warning(
      "contrasts_linfct: dropped ",
      n_failed,
      " of ",
      nrow(interaction_model_matrix),
      " proteins with failed contrasts: ",
      paste(failed_ids, collapse = ", "),
      call. = FALSE
    )
  }

  interaction_model_matrix <- interaction_model_matrix |>
    dplyr::filter(.data$classC != "logical")

  contrasts <- interaction_model_matrix |>
    dplyr::select(all_of(c(subject_Id, "contrast"))) |>
    tidyr::unnest(cols = c("contrast"))

  # take sigma and df from somewhere else.
  modelInfos <- models |>
    dplyr::select(all_of(c(subject_Id, "isSingular", "sigma.model" = "sigma", "df.residual.model" = "df.residual"))) |>

    dplyr::distinct()
  contrasts <- dplyr::inner_join(contrasts, modelInfos, by = subject_Id)
  return(ungroup(contrasts))
}


# LIMMA ----

#' Moderate p-values - limma approach
#' @export
#' @family modelling
#' @keywords internal
#'
moderated_p_limma <- function(mm, df = "df", estimate = "diff", robust = FALSE, confint = 0.95) {
  sv <- prolfqua::squeezeVarRob(mm$sigma^2, df = mm[[df]], robust = robust)

  # pior degrees of freedom are Inf
  if (all(is.infinite(sv$df.prior))) {
    sv$df.prior <- mean(mm[[df]]) * nrow(mm) / 10
  }

  sv <- tibble::as_tibble(sv)
  sv <- setNames(sv, paste0("moderated.", names(sv)))
  mm <- dplyr::bind_cols(mm, sv)
  mm <- mm |> dplyr::mutate(moderated.statistic = .data$statistic * .data$sigma / sqrt(.data$moderated.var.post))
  mm <- mm |> dplyr::mutate(moderated.df.total = !!sym(df) + .data$moderated.df.prior)
  mm <- mm |>
    dplyr::mutate(
      moderated.p.value = 2 * pt(abs(.data$moderated.statistic), df = .data$moderated.df.total, lower.tail = FALSE)
    )

  prqt <- -qt((1 - confint) / 2, df = mm$moderated.df.total)
  mm$moderated.conf.low <- mm[[estimate]] - prqt * sqrt(mm$moderated.var.post)
  mm$moderated.conf.high <- mm[[estimate]] + prqt * sqrt(mm$moderated.var.post)
  mm <- dplyr::ungroup(mm)

  return(mm)
}

#' Moderate p-value for long table
#' @param mm result of \code{\link{contrasts_linfct}}
#' @param group_by_col colnames with contrast description - default 'lhs'
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#'
#' mod <- sim_build_models_lm()
#' m <- get_complete_model_fit(mod$modelDF)
#' factor_contrasts <- linfct_factors_contrasts(m$linear_model[[1]])
#' factor_levelContrasts <- contrasts_linfct(
#'   mod$modelDF,
#'   factor_contrasts,
#'   subject_Id = "protein_Id",
#'   contrastfun = my_contrast_V2)
#'
#' mmm <- moderated_p_limma_long(factor_levelContrasts, group_by_col = "lhs")
#'
moderated_p_limma_long <- function(mm, group_by_col = "lhs", estimate = "estimate", robust = FALSE) {
  dfg <- mm |>
    dplyr::group_by(across(all_of(group_by_col))) |>
    dplyr::group_split()
  xx <- purrr::map_df(dfg, moderated_p_limma, estimate = estimate, robust = robust)
  return(xx)
}


#' adjust columns
#'
#' @export
#' @param mm data.frame with p-values to adjust
#' @param column name of column containing p-values
#' @param group_by_col column(s) to group by before adjusting (e.g. contrast), or NULL for no grouping
#' @param newname name of the new column with adjusted p-values
#' @examples
#'
#' bb <- c(runif(1000), rexp(1500,rate=5))
#' length(bb)
#' bb <- bb[bb < 1]
#' length(bb)
#' bb <- bb[1:2000]
#' hist(bb)
#' data <- data.frame(contrast = rep(LETTERS[1:5],400), p.value = bb)
#'
#' dataX <- adjust_p_values(data)
#' Adata <- dataX |> dplyr::filter(contrast == "A")
#' stopifnot(all.equal(Adata$FDR, p.adjust(Adata$p.value, method="BH")))
#' data2 <- adjust_p_values(data, group_by_col = NULL)
#' stopifnot(all.equal(data2$FDR, p.adjust(data2$p.value, method="BH")))
#'
#'
adjust_p_values <- function(
  mm,
  column = "p.value",
  group_by_col = "contrast",
  newname = "FDR"
) {
  dfg <- mm |>
    dplyr::group_by(across(all_of(group_by_col)))
  xx <- dplyr::mutate(dfg, !!newname := p.adjust(!!sym(column), method = "BH"))
  return(xx)
}


# HELPER ----

# ROPECA ----

#' p-value of protein from p.value of the median fold change peptide.
#' @param max.n limit number of peptides per protein.
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#' plot(get_p_values_pbeta(0.1,1:10,10), ylim=c(0,0.1))
#' plot(get_p_values_pbeta(0.1,1:10,3), ylim=c(0,0.1))
#' plot(get_p_values_pbeta(0.3,1:30, 3), ylim=c(0,0.1))
#' abline(h=.05,col = 2)
#' plot(seq(0.0,1.0,length=30),get_p_values_pbeta(seq(0.0,1.0,length=30),rep(10,30)))
#' abline(0,1)
#' plot(seq(0.0,1.0,length=30),get_p_values_pbeta(seq(0.0,1.0,length=30),rep(10,30),3))
#' abline(0,1)
#' testthat::expect_equal(get_p_values_pbeta(0.3,10, 3),0.216, tolerance = 1e-4)
#' testthat::expect_equal(get_p_values_pbeta(0,10, 3),0, tolerance = 1e-4)
#' testthat::expect_equal(get_p_values_pbeta(1,10, 3),1, tolerance = 1e-4)
#' testthat::expect_equal(get_p_values_pbeta(1,10, 3),get_p_values_pbeta(1,3, 10), tolerance = 1e-4)
#'
get_p_values_pbeta <- function(median.p.value, n.obs, max.n = 10) {
  n.obs <- pmin(n.obs, max.n)

  shape1 <- (n.obs / 2 + 0.5)
  shape2 <- (n.obs - (n.obs / 2 + 0.5) + 1)

  stopifnot(shape1 == shape2)
  res.p.value <- pbeta(median.p.value, shape1 = shape1, shape2 = shape2)
  return(res.p.value)
}


#' compute protein level fold changes and p.values (using beta distribution)
#' takes p-value of the scaled p-value
#'
#' @param contrasts_data data frame
#' @param contrast name of column with contrast identifier
#' @param subject_Id name of column with typically protein Id
#' @param estimate name of column with effect size estimate
#' @param statistic statistic name of column with statistic (typically t-statistics)
#' @param p.value name of column with moderated.p.value
#' @param max.n used to limit the number of peptides in probablity computation.
#' @export
#' @family modelling
#' @keywords internal
#' @return data.frame with columns
#'
#'
#' @examples
#'
#' set.seed(10)
#' nrPep <- 10000
#' nrProtein <- 800
#' p.value <- runif(nrPep)
#' estimate <- rnorm(nrPep)
#' avgAbd <- runif(nrPep)
#' protein_Id <- sample(1:800, size = nrPep,
#'   replace = TRUE, prob = dexp(seq(0,5,length = 800)))
#'
#' plot(table(table(protein_Id)))
#'
#' testdata <- data.frame(contrast = "contrast1",
#'   protein_Id = protein_Id,
#'   estimate = estimate,
#'   pseudo_estimate = estimate,
#'   p.value = p.value,
#'   avgAbd = avgAbd )
#'
#' xx30 <- summary_ROPECA_median_p.scaled(testdata,
#'                                     subject_Id = "protein_Id",
#'                                     estimate = "estimate",
#'                                     p.value = "p.value",
#'                                     max.n = 30)
#'
#' xx2 <- summary_ROPECA_median_p.scaled(testdata,
#'                                     subject_Id = "protein_Id",
#'                                     estimate = "estimate",
#'                                     p.value = "p.value",
#'                                     max.n = 1)
#'
#' testthat::expect_equal(mad(xx2$estimate, na.rm = TRUE),0.384409, tolerance = 1e-4)
#' testthat::expect_equal(median(xx2$estimate), -0.006874857, tolerance = 1e-4)
#' testthat::expect_equal(xx2$beta.based.significance[1],0.819, tolerance = 1e-3)
#' testthat::expect_equal(xx2$beta.based.significance[2],0.9234362,tolerance = 1e-3)
#'
#' # Uniform distribution
#' hist(testdata$p.value)
#' hist(xx30$median.p.scaled, breaks = 20)
#' hist(xx2$median.p.scaled, breaks = 20)
#' # shows that beta.based.significance has NO uniform distribution
#' # although H0 is true for all cases.
#'
#' hist(xx30$beta.based.significance, breaks = 20)
#' hist(xx2$beta.based.significance, breaks = 20)
#'
#' hist(xx2$median.p.value, breaks = 20)
#' hist(xx2$beta.based.significance, breaks = 20)
#' hist(estimate)
#'
summary_ROPECA_median_p.scaled <- function(
  contrasts_data,
  contrast = "contrast",
  subject_Id = "protein_Id",
  estimate = "diff",
  statistic = "statistic",
  p.value = "moderated.p.value",
  max.n = 10
) {
  nrpepsPerProt <- contrasts_data |>
    group_by(across(all_of(c(subject_Id, contrast)))) |>
    dplyr::summarize(n = dplyr::n())

  contrasts_data <- contrasts_data |>
    dplyr::mutate(
      scaled.p = ifelse(!!sym(estimate) > 0, 1 - !!sym(p.value), !!sym(p.value) - 1)
    )

  summarized.protein <- contrasts_data |>
    group_by(across(all_of(c(subject_Id, contrast)))) |>
    dplyr::summarize(
      n_not_na = n(),
      mad.estimate = mad(!!sym(estimate), na.rm = TRUE),
      estimate = median(!!sym(estimate), na.rm = TRUE),
      statistic = median(!!sym(statistic), na.rm = TRUE),
      median.p.scaled = median(.data$scaled.p, na.rm = TRUE),
      avgAbd = median(.data$avgAbd, na.rm = TRUE)
    )

  summarized.protein <- summarized.protein |>
    dplyr::mutate(median.p.value = 1 - abs(.data$median.p.scaled))

  summarized.protein <- summarized.protein |>
    dplyr::mutate(beta.based.significance = get_p_values_pbeta(.data$median.p.value, .data$n_not_na, max.n = max.n))
  summarized.protein <- summarized.protein |>
    dplyr::mutate(n.beta = pmin(.data$n_not_na, max.n))

  summarized.protein <- dplyr::inner_join(nrpepsPerProt, summarized.protein, by = c(subject_Id, contrast))

  summarized.protein$isSingular <- FALSE
  # scale it back here.
  return(ungroup(summarized.protein))
}


#' Fishers exact test on a datframe
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#' Nprot <- 1000
#' condA <- 8
#' condB <- 8
#' observedA <- sample(0:8, Nprot, replace = TRUE)
#' observedB <- sample(0:8, Nprot, replace = TRUE)
#' xb <- data.frame(observedA = observedA, observedB = observedB)
#'
#' xb$samplesA <- condA
#' xb$samplesB <- condB
#' proteinID <- unique(stringi::stri_rand_strings(Nprot + 20,5))[1:Nprot]
#' xb$proteinID <- proteinID
#' xlater <- xb
#' res <- contrasts_fisher_exact(xlater)
#'
contrasts_fisher_exact <- function(
  xb,
  observedA = "observedA",
  observedB = "observedB",
  samplesA = "samplesA",
  samplesB = "samplesB"
) {
  relativeRisk <- function(observedA, observedB, samplesA, samplesB) {
    rr <- (observedA / (observedA + observedB)) / (samplesA / (samplesA + samplesB))
    return(rr)
  }
  odsRatio <- function(observedA, observedB, samplesA, samplesB) {
    rr <- (observedA / observedB) / (samplesA / samplesB)
    return(rr)
  }
  apply_fischer <- function(proteinID, observedA, observedB, samplesA, samplesB) {
    mat <- matrix(c(observedA, samplesA - observedA, observedB, samplesB - observedB), nrow = 2)
    fisher_result <- fisher.test(mat)
    return(data.frame(
      proteinID = proteinID,
      p_value = fisher_result$p.value,
      OdsRatio = (fisher_result$estimate),
      conf.lower = (fisher_result$conf.int[1]),
      conf.higher = (fisher_result$conf.int[2])
    ))
  }

  xb$OdsRatioM <- odsRatio(
    observedA = xb[["observedA"]],
    observedB = xb[["observedB"]],
    samplesA = xb[["samplesA"]],
    samplesB = xb[["samplesB"]]
  )
  xb$relativeRiskM <- relativeRisk(
    observedA = xb[["observedA"]],
    observedB = xb[["observedB"]],
    samplesA = xb[["samplesA"]],
    samplesB = xb[["samplesB"]]
  )

  res <- vector(mode = "list", length(nrow(xb)))

  for (i in seq_len(nrow(xb))) {
    res[[i]] <- apply_fischer(
      xb[["proteinID"]][i],
      xb[["observedA"]][i],
      xb[["observedB"]][i],
      xb[["samplesA"]][i],
      xb[["samplesB"]][i]
    )
  }

  result <- dplyr::bind_rows(res)
  xx <- dplyr::inner_join(xb, result, by = c("proteinID" = "proteinID"))
  return(xx)
}
