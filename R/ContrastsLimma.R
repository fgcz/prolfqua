# StrategyLimma -----

#' R6 class for limma modelling strategy
#'
#' Analogous to \code{\link{strategy_lm}} but for limma's matrix-based pipeline.
#' Consumed by \code{\link{build_model_limma}}.
#'
#' @export
#' @family modelling
StrategyLimma <- R6::R6Class(
  "StrategyLimma",
  public = list(
    #' @field formula model formula
    formula = NULL,
    #' @field model_name name of model
    model_name = character(),
    #' @field trend logical, passed to \code{\link[limma]{eBayes}}
    trend = FALSE,
    #' @field robust logical, passed to \code{\link[limma]{eBayes}}
    robust = FALSE,
    #' @field weights either a character string (column name) or a numeric matrix
    weights = NULL,
    #' @description Create a new StrategyLimma
    #' @param modelstr model formula as string (e.g. "abundance ~ group_")
    #' @param model_name name of model
    #' @param trend logical, passed to \code{\link[limma]{eBayes}}
    #' @param robust logical, passed to \code{\link[limma]{eBayes}}
    #' @param weights either a character string (column name in annotation for
    #'   per-sample weights) or a numeric matrix (proteins x samples) passed to
    #'   \code{\link[limma]{lmFit}}. Default \code{NULL} (no weights).
    initialize = function(modelstr, model_name = "limma", trend = FALSE, robust = FALSE, weights = NULL) {
      self$formula <- as.formula(modelstr)
      self$model_name <- model_name
      self$trend <- trend
      self$robust <- robust
      self$weights <- weights
    }
  )
)

#' Create limma modelling strategy
#'
#' Wrapper that returns a \code{\link{StrategyLimma}} R6 object.
#'
#' @param modelstr model formula as string (e.g. "abundance ~ group_")
#' @param model_name name of model
#' @param trend logical, passed to \code{\link[limma]{eBayes}}
#' @param robust logical, passed to \code{\link[limma]{eBayes}}
#' @param weights either a character string (column name in annotation for
#'   per-sample weights) or a numeric matrix (proteins x samples) passed to
#'   \code{\link[limma]{lmFit}}. Default \code{NULL} (no weights).
#' @export
#' @family modelling
#' @examples
#' strat <- strategy_limma("abundance ~ group_")
#' strat$formula
#' strat$model_name
strategy_limma <- function(modelstr, model_name = "limma", trend = FALSE, robust = FALSE, weights = NULL) {
  StrategyLimma$new(modelstr, model_name = model_name, trend = trend, robust = robust, weights = weights)
}


# compute_borrowed_variance_limma -----

#' Compute borrowed variance from successful limma fits
#'
#' Extracts median sigma and df from proteins that fitted successfully
#' (no NA coefficients). Used by \code{\link{build_model_limma_impute}} to
#' replace the artificially low variance of LOD-imputed proteins.
#'
#' @param fit MArrayLM object from \code{\link[limma]{lmFit}}
#' @return list with \code{sigma} (median residual SD) and \code{df} (median
#'   residual df) from successful proteins
#' @keywords internal
#' @family modelling
compute_borrowed_variance_limma <- function(fit) {
  good <- which(
    rowSums(is.na(fit$coefficients)) == 0 &
      is.finite(fit$sigma) &
      fit$sigma > 0 &
      is.finite(fit$df.residual) &
      fit$df.residual > 1
  )
  if (length(good) == 0) {
    stop("No successful limma fits available to borrow variance from.", call. = FALSE)
  }
  list(
    sigma = stats::median(fit$sigma[good], na.rm = TRUE),
    df = stats::median(fit$df.residual[good], na.rm = TRUE)
  )
}


# .lfqdata_to_elist -----

#' Convert LFQData to limma EList with design matrix and metadata
#'
#' Shared preamble for all \code{build_model_limma*} and \code{build_model_limpa}
#' functions. Pivots LFQData to wide format, builds the design matrix from the
#' formula, resolves the subject_Id / isotopeLabel, and creates a dummy lm for
#' linfct extraction.
#'
#' @param lfqdata an \code{\link{LFQData}} object
#' @param formula a formula (with response and RHS)
#' @return a list with components:
#' \describe{
#'   \item{elist}{limma EList with \code{$E} = expression matrix}
#'   \item{expr_matrix}{the expression matrix (same as \code{elist$E})}
#'   \item{design}{the design matrix}
#'   \item{annotation}{sample-level annotation data.frame}
#'   \item{subject_Id}{character vector of hierarchy keys (possibly including isotopeLabel)}
#'   \item{rowdata}{data.frame with one row per feature, columns = subject_Id}
#'   \item{rhs_formula}{the RHS-only formula}
#'   \item{dummy_model}{a dummy \code{lm} fitted on one complete row}
#' }
#' @importFrom methods new
#' @keywords internal
.lfqdata_to_elist <- function(lfqdata, formula) {
  wide <- lfqdata$to_wide(as.matrix = TRUE)
  expr_matrix <- wide$data
  annotation <- wide$annotation
  subject_Id <- lfqdata$config$hierarchy_keys()
  rowdata <- wide$rowdata |> dplyr::select(dplyr::all_of(subject_Id))
  if (anyDuplicated(rowdata) && !is.null(lfqdata$config$isotopeLabel)) {
    rowdata <- wide$rowdata |>
      dplyr::select(dplyr::all_of(unique(c(subject_Id, lfqdata$config$isotopeLabel))))
    subject_Id <- colnames(rowdata)
  }

  rhs_formula <- formula(delete.response(terms(formula)))
  design <- model.matrix(rhs_formula, data = annotation)

  elist <- new("EList", list(E = expr_matrix))

  # Dummy model for linfct extraction: fit lm on one complete row
  complete_rows <- which(rowSums(is.na(expr_matrix)) == 0)
  if (length(complete_rows) == 0) {
    complete_rows <- which.min(rowSums(is.na(expr_matrix)))
  }
  idx <- complete_rows[1]
  dummy_data <- annotation
  dummy_data$.response <- as.numeric(expr_matrix[idx, ])
  dummy_formula <- update(rhs_formula, .response ~ .)
  dummy_model <- lm(dummy_formula, data = dummy_data)

  list(
    elist = elist,
    expr_matrix = expr_matrix,
    design = design,
    annotation = annotation,
    subject_Id = subject_Id,
    rowdata = rowdata,
    rhs_formula = rhs_formula,
    dummy_model = dummy_model
  )
}


# .resolve_weights -----

#' Resolve strategy weights to a matrix or vector for limma::lmFit
#'
#' Shared weight resolution logic for \code{build_model_limma*} functions.
#' Handles character (column name in annotation or LFQData), matrix, or NULL.
#'
#' @param lfqdata an \code{\link{LFQData}} object
#' @param strategy a strategy object with a \code{weights} field
#' @param annotation the sample-level annotation data.frame
#' @return a weight matrix, vector, or NULL
#' @keywords internal
.resolve_weights <- function(lfqdata, strategy, annotation) {
  if (is.null(strategy$weights)) {
    return(NULL)
  }

  if (is.matrix(strategy$weights)) {
    return(strategy$weights)
  }

  if (is.character(strategy$weights) && length(strategy$weights) == 1) {
    wcol <- strategy$weights
    if (wcol %in% colnames(annotation)) {
      return(annotation[[wcol]])
    }
    if (wcol %in% colnames(lfqdata$data)) {
      if (wcol %in% lfqdata$config$value_vars()) {
        wt_wide <- lfqdata$to_wide(as.matrix = TRUE, value = wcol)
        return(wt_wide$data)
      } else {
        fname_col <- lfqdata$config$table$fileName
        wt_df <- unique(lfqdata$data[, c(fname_col, wcol)])
        wt_df <- wt_df[match(annotation[[fname_col]], wt_df[[fname_col]]), ]
        return(wt_df[[wcol]])
      }
    }
  }

  NULL
}


# build_model_limma -----

#' Build limma model from LFQData
#'
#' Analogous to \code{\link{build_model}} but uses limma's matrix-based pipeline.
#' Takes an LFQData object and a strategy from \code{\link{strategy_limma}},
#' pivots data to wide format, fits with \code{\link[limma]{lmFit}}, and returns
#' a \code{\link{ModelLimma}} object.
#'
#' @param lfqdata an \code{\link{LFQData}} object
#' @param strategy output of \code{\link{strategy_limma}}
#' @param modelName name of model (default from strategy)
#' @return a \code{\link{ModelLimma}} object
#' @export
#' @family modelling
#' @examples
#'
#' istar <- sim_lfq_data_protein_config(Nprot = 50)
#' lProt <- LFQData$new(istar$data, istar$config)
#' lProt$rename_response("transformedIntensity")
#'
#' strat <- strategy_limma("transformedIntensity ~ group_")
#' mod_limma <- build_model_limma(lProt, strat)
#' mod_limma$get_coefficients()
#' mod_limma$get_anova()
#'
build_model_limma <- function(lfqdata, strategy, modelName = strategy$model_name) {
  setup <- .lfqdata_to_elist(lfqdata, strategy$formula)
  wt <- .resolve_weights(lfqdata, strategy, setup$annotation)
  fit <- limma::lmFit(setup$expr_matrix, setup$design, weights = wt)

  ModelLimma$new(
    fit = fit,
    design = setup$design,
    formula = strategy$formula,
    subject_Id = setup$subject_Id,
    modelName = modelName,
    rowdata = setup$rowdata,
    trend = strategy$trend,
    robust = strategy$robust,
    dummy_model = setup$dummy_model,
    p.adjust = prolfqua::adjust_p_values
  )
}


# build_model_limma_impute -----

#' Build limma model with LOD imputation for failed proteins
#'
#' Analogous to \code{\link{build_model_impute}} but for limma's matrix-based
#' pipeline. Fits all proteins with \code{\link[limma]{lmFit}}, identifies
#' proteins with NA coefficients (typically from entire missing groups), imputes
#' their missing values with the limit of detection (LOD), refits, and replaces
#' the variance with a borrowed estimate from successful proteins.
#'
#' The LOD imputation gives plausible coefficients (fold change direction),
#' while the borrowed sigma and corrected degrees of freedom ensure that
#' inference is not artificially precise from the constant imputation.
#'
#' @param lfqdata an \code{\link{LFQData}} object (aggregated to protein level)
#' @param strategy output of \code{\link{strategy_limma}}
#' @param modelName name of model (default: strategy name + "Imputed")
#' @param lod numeric limit of detection; if NULL, auto-computed from data
#'   via \code{\link{MissingHelpers}}
#' @param df_method how to set degrees of freedom for imputed proteins:
#'   \code{"observed"} (default) uses \code{max(n_observed - p, 1)} where
#'   \code{n_observed} counts only non-missing values;
#'   \code{"borrowed"} uses the median df from successful fits
#' @return a \code{\link{ModelLimma}} object with a hybrid fit
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_protein_config(Nprot = 50, weight_missing = 0.5)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$rename_response("transformedIntensity")
#'
#' strat <- strategy_limma("transformedIntensity ~ group_")
#' mod <- build_model_limma_impute(lfqdata, strat)
#' mod$get_coefficients()
#'
build_model_limma_impute <- function(
  lfqdata,
  strategy,
  modelName = paste0(strategy$model_name, "Imputed"),
  lod = NULL,
  df_method = c("observed", "borrowed")
) {
  df_method <- match.arg(df_method)

  setup <- .lfqdata_to_elist(lfqdata, strategy$formula)
  expr_matrix <- setup$expr_matrix
  design <- setup$design
  p <- ncol(design)
  wt <- .resolve_weights(lfqdata, strategy, setup$annotation)

  # Step 1: Fit on original data (with NAs)
  fit_na <- limma::lmFit(expr_matrix, design, weights = wt)

  # Step 2: Identify failed proteins (any NA coefficient)
  failed <- which(rowSums(is.na(fit_na$coefficients)) > 0)

  if (length(failed) == 0) {
    return(ModelLimma$new(
      fit = fit_na,
      design = design,
      formula = strategy$formula,
      subject_Id = setup$subject_Id,
      modelName = modelName,
      rowdata = setup$rowdata,
      trend = strategy$trend,
      robust = strategy$robust,
      dummy_model = setup$dummy_model,
      p.adjust = prolfqua::adjust_p_values
    ))
  }

  # Step 3: Compute LOD if not provided
  if (is.null(lod)) {
    mh <- MissingHelpers$new(lfqdata$data, lfqdata$config)
    lod <- mh$get_LOD()
  }

  # Step 4: Compute borrowed variance from successful proteins
  borrowed <- compute_borrowed_variance_limma(fit_na)

  # Step 5: Impute expression matrix with LOD
  na_mask <- is.na(expr_matrix)
  expr_imputed <- expr_matrix
  expr_imputed[na_mask] <- lod
  expr_imputed <- pmax(expr_imputed, lod)

  # Also impute weight matrix: NA weights for imputed positions get weight = 1
  wt_imputed <- wt
  if (is.matrix(wt_imputed)) {
    wt_imputed[na_mask] <- 1
  }

  # Step 6: Fit on imputed data
  fit_lod <- limma::lmFit(expr_imputed, design, weights = wt_imputed)

  # Step 7: Create hybrid fit — replace failed rows in fit_na
  fit_na$coefficients[failed, ] <- fit_lod$coefficients[failed, ]
  fit_na$stdev.unscaled[failed, ] <- fit_lod$stdev.unscaled[failed, ]
  fit_na$sigma[failed] <- borrowed$sigma
  fit_na$Amean[failed] <- fit_lod$Amean[failed]
  # CRITICAL: df correction — never use fit_lod$df.residual (counts imputed as real)
  n_observed <- rowSums(!is.na(expr_matrix))
  if (df_method == "observed") {
    fit_na$df.residual[failed] <- pmax(n_observed[failed] - p, 1)
  } else {
    fit_na$df.residual[failed] <- borrowed$df
  }

  # Dummy model from imputed data (all rows complete after imputation)
  dummy_data <- setup$annotation
  dummy_data$.response <- as.numeric(expr_imputed[1, ])
  dummy_formula <- update(setup$rhs_formula, .response ~ .)
  dummy_model <- lm(dummy_formula, data = dummy_data)

  ModelLimma$new(
    fit = fit_na,
    design = design,
    formula = strategy$formula,
    subject_Id = setup$subject_Id,
    modelName = modelName,
    rowdata = setup$rowdata,
    trend = strategy$trend,
    robust = strategy$robust,
    dummy_model = dummy_model,
    p.adjust = prolfqua::adjust_p_values
  )
}


# build_model_limma_voom -----

#' Build limma model with vooma precision weights (proteomics)
#'
#' Estimates observation-level precision weights from a mean-variance trend
#' (vooma) and fits a weighted least squares model via \code{\link[limma]{lmFit}}.
#' For proteomics data that is already log2-transformed.
#'
#' When \code{strategy$weights} is set (e.g. to \code{nr_children}), these
#' external weights enter the preliminary fit (so the trend estimation accounts
#' for measurement precision) and are multiplied element-wise with the vooma
#' weights for the final fit. See the voom integration notes in
#' \code{TODO/TODO_limma_voom_integration.md} for the mathematical basis.
#'
#' @param lfqdata an \code{\link{LFQData}} object (aggregated to protein level)
#' @param strategy output of \code{\link{strategy_limma}}
#' @param modelName name of model (default from strategy)
#' @param span lowess smoother span for mean-variance trend (default 0.5)
#' @param plot logical; if TRUE, plot the mean-variance trend
#' @return a \code{\link{ModelLimma}} object
#' @export
#' @family modelling
#' @examples
#'
#' istar <- sim_lfq_data_protein_config(Nprot = 50)
#' lProt <- LFQData$new(istar$data, istar$config)
#' lProt$rename_response("transformedIntensity")
#'
#' strat <- strategy_limma("transformedIntensity ~ group_")
#' mod <- build_model_limma_voom(lProt, strat)
#' mod$get_coefficients()
#'
build_model_limma_voom <- function(
  lfqdata,
  strategy,
  modelName = strategy$model_name,
  span = 0.5,
  plot = FALSE
) {
  setup <- .lfqdata_to_elist(lfqdata, strategy$formula)
  expr_matrix <- setup$expr_matrix
  design <- setup$design
  ext_wt <- .resolve_weights(lfqdata, strategy, setup$annotation)

  # Step 1: Preliminary fit (with external weights if available)
  fit0 <- limma::lmFit(expr_matrix, design, weights = ext_wt)

  # Step 2: Lowess trend on (mean expression, sqrt(sigma))
  sx <- fit0$Amean
  sy <- sqrt(fit0$sigma)
  l <- stats::lowess(sx, sy, f = span)
  f <- stats::approxfun(l, rule = 2, ties = list("ordered", mean))

  # Step 3: Observation-level vooma weights
  fitted_vals <- fit0$coefficients %*% t(fit0$design)
  w_voom <- 1 / f(fitted_vals)^4
  dim(w_voom) <- dim(fitted_vals)

  # Step 4: Combine with external weights
  if (!is.null(ext_wt)) {
    if (is.matrix(ext_wt)) {
      w_combined <- w_voom * ext_wt
    } else {
      w_combined <- w_voom * rep(ext_wt, each = nrow(w_voom))
    }
  } else {
    w_combined <- w_voom
  }

  # Optional diagnostic plot
  if (plot) {
    graphics::plot(
      sx,
      sy,
      xlab = "Average log2 expression",
      ylab = expression(sqrt(sigma)),
      main = "vooma: Mean-variance trend",
      pch = 16,
      cex = 0.3
    )
    graphics::lines(l, col = "red", lwd = 2)
  }

  # Step 5: Final WLS fit with combined weights
  fit <- limma::lmFit(expr_matrix, design, weights = w_combined)

  ModelLimma$new(
    fit = fit,
    design = design,
    formula = strategy$formula,
    subject_Id = setup$subject_Id,
    modelName = modelName,
    rowdata = setup$rowdata,
    trend = strategy$trend,
    robust = strategy$robust,
    dummy_model = setup$dummy_model,
    p.adjust = prolfqua::adjust_p_values
  )
}


# build_model_limma_voom_impute -----

#' Build limma-voom model with LOD imputation for failed proteins
#'
#' Combines vooma precision weights with LOD imputation for proteins that have
#' entire missing groups (NA coefficients). Mirrors
#' \code{\link{build_model_limma_impute}} but uses vooma weights.
#'
#' @param lfqdata an \code{\link{LFQData}} object (aggregated to protein level)
#' @param strategy output of \code{\link{strategy_limma}}
#' @param modelName name of model (default: strategy name + "Imputed")
#' @param lod numeric limit of detection; if NULL, auto-computed from data
#' @param df_method how to set degrees of freedom for imputed proteins:
#'   \code{"observed"} (default) uses \code{max(n_observed - p, 1)};
#'   \code{"borrowed"} uses the median df from successful fits
#' @param span lowess smoother span for vooma trend (default 0.5)
#' @param plot logical; if TRUE, plot the mean-variance trend
#' @return a \code{\link{ModelLimma}} object with a hybrid fit
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_protein_config(Nprot = 50, weight_missing = 0.5)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$rename_response("transformedIntensity")
#'
#' strat <- strategy_limma("transformedIntensity ~ group_")
#' mod <- build_model_limma_voom_impute(lfqdata, strat)
#' mod$get_coefficients()
#'
build_model_limma_voom_impute <- function(
  lfqdata,
  strategy,
  modelName = paste0(strategy$model_name, "Imputed"),
  lod = NULL,
  df_method = c("observed", "borrowed"),
  span = 0.5,
  plot = FALSE
) {
  df_method <- match.arg(df_method)

  setup <- .lfqdata_to_elist(lfqdata, strategy$formula)
  expr_matrix <- setup$expr_matrix
  design <- setup$design
  p <- ncol(design)
  ext_wt <- .resolve_weights(lfqdata, strategy, setup$annotation)

  # Helper: compute vooma weights and final fit from an expression matrix
  .voom_fit <- function(expr, wt_ext) {
    fit0 <- limma::lmFit(expr, design, weights = wt_ext)
    sx <- fit0$Amean
    sy <- sqrt(fit0$sigma)
    l <- stats::lowess(sx, sy, f = span)
    f_trend <- stats::approxfun(l, rule = 2, ties = list("ordered", mean))
    fitted_vals <- fit0$coefficients %*% t(fit0$design)
    w_voom <- 1 / f_trend(fitted_vals)^4
    dim(w_voom) <- dim(fitted_vals)
    if (!is.null(wt_ext)) {
      if (is.matrix(wt_ext)) {
        w_combined <- w_voom * wt_ext
      } else {
        w_combined <- w_voom * rep(wt_ext, each = nrow(w_voom))
      }
    } else {
      w_combined <- w_voom
    }
    if (plot) {
      graphics::plot(
        sx,
        sy,
        xlab = "Average log2 expression",
        ylab = expression(sqrt(sigma)),
        main = "vooma: Mean-variance trend",
        pch = 16,
        cex = 0.3
      )
      graphics::lines(l, col = "red", lwd = 2)
    }
    limma::lmFit(expr, design, weights = w_combined)
  }

  # Step 1: Fit on original data with vooma weights
  fit_na <- .voom_fit(expr_matrix, ext_wt)

  # Step 2: Identify failed proteins (any NA coefficient)
  failed <- which(rowSums(is.na(fit_na$coefficients)) > 0)

  if (length(failed) == 0) {
    return(ModelLimma$new(
      fit = fit_na,
      design = design,
      formula = strategy$formula,
      subject_Id = setup$subject_Id,
      modelName = modelName,
      rowdata = setup$rowdata,
      trend = strategy$trend,
      robust = strategy$robust,
      dummy_model = setup$dummy_model,
      p.adjust = prolfqua::adjust_p_values
    ))
  }

  # Step 3: Compute LOD if not provided
  if (is.null(lod)) {
    mh <- MissingHelpers$new(lfqdata$data, lfqdata$config)
    lod <- mh$get_LOD()
  }

  # Step 4: Borrow variance from successful proteins
  borrowed <- compute_borrowed_variance_limma(fit_na)

  # Step 5: Impute expression matrix with LOD
  na_mask <- is.na(expr_matrix)
  expr_imputed <- expr_matrix
  expr_imputed[na_mask] <- lod
  expr_imputed <- pmax(expr_imputed, lod)

  # Impute external weight matrix: NA positions get weight = 1
  ext_wt_imputed <- ext_wt
  if (is.matrix(ext_wt_imputed)) {
    ext_wt_imputed[na_mask] <- 1
  }

  # Step 6: Fit on imputed data with vooma weights
  fit_lod <- .voom_fit(expr_imputed, ext_wt_imputed)

  # Step 7: Hybrid fit — replace failed rows
  fit_na$coefficients[failed, ] <- fit_lod$coefficients[failed, ]
  fit_na$stdev.unscaled[failed, ] <- fit_lod$stdev.unscaled[failed, ]
  fit_na$sigma[failed] <- borrowed$sigma
  fit_na$Amean[failed] <- fit_lod$Amean[failed]
  n_observed <- rowSums(!is.na(expr_matrix))
  if (df_method == "observed") {
    fit_na$df.residual[failed] <- pmax(n_observed[failed] - p, 1)
  } else {
    fit_na$df.residual[failed] <- borrowed$df
  }

  # Dummy model from imputed data (all rows complete after imputation)
  dummy_data <- setup$annotation
  dummy_data$.response <- as.numeric(expr_imputed[1, ])
  dummy_formula <- update(setup$rhs_formula, .response ~ .)
  dummy_model <- lm(dummy_formula, data = dummy_data)

  ModelLimma$new(
    fit = fit_na,
    design = design,
    formula = strategy$formula,
    subject_Id = setup$subject_Id,
    modelName = modelName,
    rowdata = setup$rowdata,
    trend = strategy$trend,
    robust = strategy$robust,
    dummy_model = dummy_model,
    p.adjust = prolfqua::adjust_p_values
  )
}


# ModelLimma -----

#' R6 class representing a limma modelling result
#'
#' Same API as \code{\link{Model}}: \code{get_anova()}, \code{get_coefficients()},
#' \code{coef_histogram()}, \code{coef_volcano()}, \code{anova_histogram()}.
#'
#' @export
#' @family modelling
#' @examples
#'
#' istar <- sim_lfq_data_protein_config(Nprot = 50)
#' lProt <- LFQData$new(istar$data, istar$config)
#' lProt$rename_response("transformedIntensity")
#' strat <- strategy_limma("transformedIntensity ~ group_")
#' mod <- build_model_limma(lProt, strat)
#'
#' coeffs <- mod$get_coefficients()
#' head(coeffs)
#' anova_tbl <- mod$get_anova()
#' head(anova_tbl)
#' mod$coef_histogram()
#' mod$coef_volcano()
#' mod$anova_histogram()
#'
ModelLimma <- R6::R6Class(
  "ModelLimma",
  inherit = ModelInterface,
  public = list(
    #' @field fit limma MArrayLM object from lmFit
    fit = NULL,
    #' @field design design matrix
    design = NULL,
    #' @field formula model formula
    formula = NULL,
    #' @field subject_Id protein ID column name(s)
    subject_Id = character(),
    #' @field modelName model name
    modelName = character(),
    #' @field rowdata protein ID mapping from to_wide()$rowdata
    rowdata = NULL,
    #' @field trend passed to eBayes
    trend = FALSE,
    #' @field robust passed to eBayes
    robust = FALSE,
    #' @field dummy_model one fitted lm for linfct extraction
    dummy_model = NULL,
    #' @field p.adjust function to adjust p-values
    p.adjust = NULL,
    #' @description
    #' initialize ModelLimma
    #' @param fit limma MArrayLM from lmFit
    #' @param design design matrix
    #' @param formula model formula
    #' @param subject_Id protein ID column name(s)
    #' @param modelName model name
    #' @param rowdata protein ID mapping
    #' @param trend passed to eBayes
    #' @param robust passed to eBayes
    #' @param dummy_model one fitted lm for linfct extraction
    #' @param p.adjust function to adjust p-values
    initialize = function(
      fit,
      design,
      formula,
      subject_Id,
      modelName,
      rowdata,
      trend = FALSE,
      robust = FALSE,
      dummy_model = NULL,
      p.adjust = prolfqua::adjust_p_values
    ) {
      self$fit <- fit
      self$design <- design
      self$formula <- formula
      self$subject_Id <- subject_Id
      self$modelName <- modelName
      self$rowdata <- rowdata
      self$trend <- trend
      self$robust <- robust
      self$dummy_model <- dummy_model
      self$p.adjust <- p.adjust
    },
    #' @description
    #' return model coefficient table in long format
    #' @return data.frame
    get_coefficients = function() {
      fit_eb <- limma::eBayes(self$fit, trend = self$trend, robust = self$robust)

      coef_mat <- fit_eb$coefficients
      se_mat <- fit_eb$stdev.unscaled * fit_eb$s2.post^0.5
      t_mat <- fit_eb$t
      p_mat <- fit_eb$p.value

      ncoef <- ncol(coef_mat)
      coef_names <- colnames(coef_mat)
      ngenes <- nrow(coef_mat)

      res_list <- vector("list", ncoef)
      for (j in seq_len(ncoef)) {
        df_j <- self$rowdata
        df_j$factor <- coef_names[j]
        df_j$Estimate <- coef_mat[, j]
        df_j$Std..Error <- se_mat[, j]
        df_j$t.value <- t_mat[, j]
        df_j$Pr...t.. <- p_mat[, j]
        res_list[[j]] <- df_j
      }
      result <- dplyr::bind_rows(res_list)
      return(result)
    },
    #' @description
    #' return anova table (F-test per protein across all non-intercept coefficients)
    #' @return data.frame
    get_anova = function() {
      fit_eb <- limma::eBayes(self$fit, trend = self$trend, robust = self$robust)

      coef_names <- colnames(fit_eb$coefficients)
      non_intercept <- which(coef_names != "(Intercept)")

      if (length(non_intercept) == 0) {
        warning("No non-intercept coefficients for ANOVA.")
        return(data.frame())
      }

      # Use limma topTable with all non-intercept coefficients for F-test
      tt <- limma::topTable(fit_eb, coef = non_intercept, number = Inf, sort.by = "none")

      result <- self$rowdata
      result$F.value <- tt$F
      result$p.value <- tt$P.Value

      # create a factor column describing the tested term
      terms_tested <- paste(coef_names[non_intercept], collapse = "+")
      result$factor <- terms_tested

      result <- self$p.adjust(result, column = "p.value", group_by_col = "factor")
      return(dplyr::ungroup(result))
    },
    #' @description
    #' histogram of model coefficient p-values
    coef_histogram = function() {
      model_coeff <- self$get_coefficients()
      model_coeff <- tidyr::unite(model_coeff, "subject_Id", self$subject_Id)
      fname <- paste0("Coef_Histogram_", self$modelName, ".pdf")
      p <- ggplot(data = model_coeff, aes(x = .data$Pr...t.., group = .data$factor)) +
        geom_histogram(breaks = seq(0, 1, by = 0.05)) +
        facet_wrap(~factor)
      return(list(plot = p, name = fname))
    },
    #' @description
    #' volcano plot of non-intercept coefficients
    coef_volcano = function() {
      model_coeff <- self$get_coefficients()
      model_coeff <- tidyr::unite(model_coeff, "subject_Id", self$subject_Id)
      fname <- paste0("Coef_VolcanoPlot_", self$modelName, ".pdf")
      p <- model_coeff |>
        dplyr::filter(.data$factor != "(Intercept)") |>
        prolfqua::multigroup_volcano(
          effect = "Estimate",
          significance = "Pr...t..",
          contrast = "factor",
          label = "subject_Id",
          xintercept = c(-1, 1),
          colour = NULL
        )
      return(list(plot = p, name = fname))
    },
    #' @description
    #' pairs plot of coefficients
    coef_pairs = function() {
      model_coeff <- self$get_coefficients()
      model_coeff <- tidyr::unite(model_coeff, "subject_Id", self$subject_Id)
      forPairs <- model_coeff |>
        dplyr::select(all_of(c("subject_Id", "factor", "Estimate"))) |>
        tidyr::pivot_wider(names_from = "factor", values_from = "Estimate")
      fname <- paste0("Coef_Pairsplot_", self$modelName, ".pdf")
      return(list(plot = forPairs, name = fname))
    },
    #' @description
    #' histogram of ANOVA F-test p-values
    #' @param what show either "p.value" or "FDR"
    anova_histogram = function(what = c("p.value", "FDR")) {
      what <- match.arg(what)
      model_anova <- self$get_anova()
      fname <- paste0("Anova_p.values_", self$modelName, ".pdf")
      p <- model_anova |>
        ggplot(aes(x = !!sym(what), group = .data$factor)) +
        geom_histogram(breaks = seq(0, 1, by = 0.05)) +
        facet_wrap(~factor)
      return(list(plot = p, name = fname))
    }
  )
)


# ContrastsLimma -----

#' Limma-based contrasts (direct limma pipeline)
#'
#' Uses limma's \code{contrasts.fit} + \code{eBayes} pipeline directly,
#' rather than fitting per-protein lm models and then moderating.
#' Inherits from \code{\link{ContrastsInterface}} with the same API as
#' \code{\link{Contrasts}}.
#'
#' @export
#' @family modelling
#' @examples
#'
#' istar <- sim_lfq_data_protein_config(Nprot = 50)
#' lProt <- LFQData$new(istar$data, istar$config)
#' lProt$rename_response("transformedIntensity")
#'
#' strat <- strategy_limma("transformedIntensity ~ group_")
#' mod_limma <- build_model_limma(lProt, strat)
#'
#' Contr <- c("dil.b_vs_a" = "group_A - group_Ctrl")
#' contr_limma <- ContrastsLimma$new(mod_limma, Contr)
#' res <- contr_limma$get_contrasts()
#' head(res)
#' stopifnot(all(c("diff", "FDR", "p.value", "statistic") %in% colnames(res)))
#'
#' # Compare with prolfqua's own pipeline
#' modelFunction <- strategy_lm("transformedIntensity ~ group_")
#' mod <- build_model(lProt, modelFunction)
#' contr_prolfqua <- Contrasts$new(mod, Contr)
#' res_prolfqua <- contr_prolfqua$get_contrasts()
#'
#' # fold changes should be very similar
#' merged <- dplyr::inner_join(
#'   dplyr::select(res, protein_Id, diff_limma = diff),
#'   dplyr::select(res_prolfqua, protein_Id, diff_prolfqua = diff),
#'   by = "protein_Id")
#' stopifnot(cor(merged$diff_limma, merged$diff_prolfqua, use = "complete.obs") > 0.99)
#'
#' # Plotter works
#' pl <- contr_limma$get_Plotter()
#'
#' # to_wide works
#' wide <- contr_limma$to_wide()
#' head(wide)
#'
#' # merge_contrasts_results works
#' Contr2 <- c("dil.b_vs_a" = "group_A - group_Ctrl")
#' csi <- ContrastsMissing$new(lProt, contrasts = Contr2)
#' merged_res <- merge_contrasts_results(contr_limma, csi)
#'
ContrastsLimma <- R6::R6Class(
  "ContrastsLimma",
  inherit = ContrastsInterface,
  public = list(
    #' @field model ModelLimma object
    model = NULL,
    #' @field contrasts named character vector of contrasts
    contrasts = character(),
    #' @field modelName model name
    modelName = character(),
    #' @field subject_Id columns with subject_Id (proteinID)
    subject_Id = character(),
    #' @field p.adjust function to adjust p-values
    p.adjust = NULL,
    #' @field contrast_result cached contrast results
    contrast_result = NULL,
    #' @field eBayes logical, apply limma eBayes moderation (default TRUE).
    #'   Set to FALSE to return raw/unmoderated statistics, e.g. for downstream
    #'   DEqMS moderation via \code{\link{ContrastsModeratedDEqMS}}.
    eBayes = TRUE,
    #' @description
    #' initialize ContrastsLimma
    #' @param model a \code{\link{ModelLimma}} object
    #' @param contrasts named character vector of contrasts
    #' @param p.adjust function to adjust p-values
    #' @param modelName name of the contrast method
    #' @param eBayes logical, apply limma eBayes moderation (default TRUE).
    #'   Set to FALSE to return raw/unmoderated statistics suitable for
    #'   wrapping with \code{\link{ContrastsModeratedDEqMS}}.
    initialize = function(model, contrasts, p.adjust = prolfqua::adjust_p_values, modelName = NULL, eBayes = TRUE) {
      self$model <- model
      self$contrasts <- contrasts
      self$subject_Id <- model$subject_Id
      self$eBayes <- eBayes
      self$modelName <- modelName %||% if (eBayes) "limma" else "limma_raw"
      self$p.adjust <- p.adjust
    },
    #' @description
    #' get both sides of contrasts
    get_contrast_sides = function() {
      tt <- self$contrasts[grep("-", self$contrasts)]
      tt <- tibble(contrast = names(tt), rhs = tt)
      tt <- tt |>
        mutate(rhs = gsub("[` ]", "", rhs)) |>
        tidyr::separate(rhs, c("group_1", "group_2"), sep = "-")
      return(tt)
    },
    #' @description
    #' get linear functions from contrasts
    #' @param global ignored (for API compatibility)
    #' @param avg logical, also compute avgAbd linfct
    get_linfct = function(global = TRUE, avg = TRUE) {
      linfct <- .linfct(self$model$dummy_model, self$contrasts, avg = avg)
      return(linfct)
    },
    #' @description
    #' get table with contrast estimates via limma pipeline
    #' @param all should all columns be returned (default FALSE)
    #' @return data.frame with contrasts
    get_contrasts = function(all = FALSE) {
      if (!is.null(self$contrast_result)) {
        return(self$contrast_result)
      }

      # Build the contrast matrix via linfct_from_model + linfct_matrix_contrasts
      linfct <- linfct_from_model(self$model$dummy_model, as_list = FALSE)
      linfct <- unique(linfct) # needed for single factor models

      # Add avg contrasts for avgAbd computation
      namtmp <- paste0("avg_", names(self$contrasts))
      cntr_avg <- paste0("(", gsub(" - ", " + ", self$contrasts), ")/2")
      names(cntr_avg) <- namtmp
      all_contrasts <- c(self$contrasts, cntr_avg)

      linfct_A <- linfct_matrix_contrasts(linfct, all_contrasts)
      # linfct_A: rows = contrast names, cols = model coefficients

      # Split into difference contrasts and avg contrasts
      diff_names <- names(self$contrasts)
      avg_names <- namtmp

      # Transpose for limma: rows = coefficients, cols = contrasts
      contrast_matrix <- t(linfct_A[diff_names, , drop = FALSE])

      # limma pipeline: contrasts.fit, optionally + eBayes
      fit2 <- limma::contrasts.fit(self$model$fit, contrast_matrix)

      if (self$eBayes) {
        fit2 <- limma::eBayes(fit2, trend = self$model$trend, robust = self$model$robust)
      }

      # Extract results per contrast
      res_list <- vector("list", length(diff_names))
      for (i in seq_along(diff_names)) {
        df_i <- self$model$rowdata
        df_i$contrast <- diff_names[i]

        if (self$eBayes) {
          tt <- limma::topTable(fit2, coef = i, number = Inf, sort.by = "none", confint = TRUE)
          df_i$diff <- tt$logFC
          df_i$std.error <- sqrt(fit2$s2.post) * fit2$stdev.unscaled[, i]
          df_i$statistic <- tt$t
          df_i$p.value <- tt$P.Value
          df_i$sigma <- sqrt(fit2$s2.post)
          df_i$df <- fit2$df.total
          df_i$conf.low <- tt$CI.L
          df_i$conf.high <- tt$CI.R
        } else {
          # Raw/unmoderated statistics from the contrasts.fit object
          sigma_raw <- fit2$sigma # per-protein residual SD
          df_raw <- fit2$df.residual # per-protein residual df
          diff_i <- fit2$coefficients[, i]
          se_i <- sigma_raw * fit2$stdev.unscaled[, i]
          t_raw <- diff_i / se_i
          p_raw <- 2 * pt(abs(t_raw), df = df_raw, lower.tail = FALSE)
          alpha <- 0.05
          prqt <- -qt(alpha / 2, df = df_raw)

          df_i$diff <- diff_i
          df_i$std.error <- se_i
          df_i$statistic <- t_raw
          df_i$p.value <- p_raw
          df_i$sigma <- sigma_raw
          df_i$df <- df_raw
          df_i$conf.low <- diff_i - prqt * se_i
          df_i$conf.high <- diff_i + prqt * se_i
        }

        res_list[[i]] <- df_i
      }
      contrast_result <- dplyr::bind_rows(res_list)

      # Compute avgAbd from the avg linfct applied to the fit coefficients
      avg_matrix <- t(linfct_A[avg_names, , drop = FALSE])
      avg_vals <- self$model$fit$coefficients %*% avg_matrix
      # avg_vals: rows = proteins, cols = avg contrasts

      avg_df_list <- vector("list", length(avg_names))
      for (i in seq_along(avg_names)) {
        df_avg <- self$model$rowdata
        df_avg$contrast <- diff_names[i] # map back to original contrast name
        df_avg$avgAbd <- avg_vals[, i]
        avg_df_list[[i]] <- df_avg
      }
      avg_df <- dplyr::bind_rows(avg_df_list)

      contrast_result <- dplyr::left_join(
        contrast_result,
        avg_df,
        by = c(self$subject_Id, "contrast")
      )

      # Adjust p-values per contrast
      contrast_result <- self$p.adjust(contrast_result, column = "p.value", group_by_col = "contrast", newname = "FDR")
      contrast_result <- contrast_result |> dplyr::relocate("FDR", .after = "diff")
      contrast_result <- dplyr::mutate(contrast_result, modelName = self$modelName, .before = 1)
      contrast_result <- dplyr::ungroup(contrast_result)
      self$contrast_result <- contrast_result

      stopifnot(all(super$column_description()$column_name %in% colnames(contrast_result)))
      return(contrast_result)
    },
    #' @description
    #' return \code{\link{ContrastsPlotter}}
    #' @param FCthreshold fold change threshold to show in plots
    #' @param FDRthreshold FDR threshold to show in plots
    #' @return \code{\link{ContrastsPlotter}}
    get_Plotter = function(FCthreshold = 1, FDRthreshold = 0.1) {
      contrast_result <- self$get_contrasts()
      res <- ContrastsPlotter$new(
        contrast_result,
        subject_Id = self$subject_Id,
        fcthresh = FCthreshold,
        volcano = list(
          list(score = "p.value", thresh = FDRthreshold),
          list(score = "FDR", thresh = FDRthreshold)
        ),
        histogram = list(
          list(score = "p.value", xlim = c(0, 1, 0.05)),
          list(score = "FDR", xlim = c(0, 1, 0.05))
        ),
        score = list(list(score = "statistic", thresh = 5)),
        modelName = "modelName",
        diff = "diff",
        contrast = "contrast"
      )
      return(res)
    },
    #' @description
    #' convert to wide format
    #' @param columns value column default p.value
    #' @return data.frame
    to_wide = function(columns = c("p.value", "FDR", "statistic")) {
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(
        contrast_minimal,
        subject_Id = self$subject_Id,
        columns = c("diff", columns),
        contrast = "contrast"
      )
      return(contrasts_wide)
    }
  )
)
