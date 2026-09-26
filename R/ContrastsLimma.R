# StrategyLimma -----

#' R6 class for limma modelling strategy
#'
#' Analogous to \code{\link{strategy_lm}} but for limma's matrix-based pipeline.
#' Consumed by \code{\link{build_model_limma}}.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' strat <- StrategyLimma$new("abundance ~ group_")
#' strat$formula
#' strat$model_name
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
    #' @param modelstr,model_name,trend,robust,weights see \code{\link{strategy_limma}}
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
#' @return The computed result.
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
#' formula, resolves the subject_id / isotopeLabel, and creates a dummy lm for
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
#'   \item{subject_id}{character vector of hierarchy keys (possibly including isotopeLabel)}
#'   \item{rowdata}{data.frame with one row per feature, columns = subject_id}
#'   \item{rhs_formula}{the RHS-only formula}
#'   \item{dummy_model}{a dummy \code{lm} fitted on one complete row}
#' }
#' @importFrom methods new
#' @keywords internal
.lfqdata_to_elist <- function(lfqdata, formula) {
  wide <- lfqdata$data_wide(as.matrix = TRUE)
  expr_matrix <- wide$data
  annotation <- wide$annotation
  subject_id <- lfqdata$hierarchy_keys()
  rowdata <- wide$rowdata |> dplyr::select(dplyr::all_of(subject_id))
  if (anyDuplicated(rowdata) && !is.null(lfqdata$isotope_label())) {
    rowdata <- wide$rowdata |>
      dplyr::select(dplyr::all_of(unique(c(subject_id, lfqdata$isotope_label()))))
    subject_id <- colnames(rowdata)
  }

  rhs_formula <- formula(delete.response(terms(formula)))
  design <- model.matrix(rhs_formula, data = annotation)

  elist <- new("EList", list(E = expr_matrix))

  # Dummy model for linfct extraction: fit lm on one complete row
  complete_rows <- which(rowSums(is.na(expr_matrix)) == 0)
  if (length(complete_rows) == 0) {
    complete_rows <- which.min(rowSums(is.na(expr_matrix)))
  }
  dummy_model <- .limma_dummy_model(annotation, rhs_formula, expr_matrix[complete_rows[1], ])

  list(
    elist = elist,
    expr_matrix = expr_matrix,
    design = design,
    annotation = annotation,
    subject_id = subject_id,
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
    if (wcol %in% colnames(lfqdata$data_long())) {
      if (wcol %in% lfqdata$get_config()$value_vars()) {
        wt_wide <- lfqdata$data_wide(as.matrix = TRUE, value = wcol)
        return(wt_wide$data)
      } else {
        fname_col <- lfqdata$file_name()
        wt_df <- unique(lfqdata$data_long()[, c(fname_col, wcol)])
        wt_df <- wt_df[match(annotation[[fname_col]], wt_df[[fname_col]]), ]
        return(wt_df[[wcol]])
      }
    }
  }

  NULL
}

.new_model_limma <- function(
  fit,
  setup,
  strategy,
  model_name,
  dummy_model = setup$dummy_model,
  imputed_proteins = character(0)
) {
  ModelLimma$new(
    fit = fit,
    subject_id = setup$subject_id,
    model_name = model_name,
    rowdata = setup$rowdata,
    trend = strategy$trend,
    robust = strategy$robust,
    dummy_model = dummy_model,
    p.adjust = prolfqua::adjust_p_values,
    imputed_proteins = imputed_proteins
  )
}

# lm of one expression row on the design, used to derive the linear functions of the contrasts.
.limma_dummy_model <- function(annotation, rhs_formula, response) {
  dummy_data <- annotation
  dummy_data$.response <- as.numeric(response)
  dummy_formula <- update(rhs_formula, .response ~ .)
  lm(dummy_formula, data = dummy_data)
}

.impute_limma_data <- function(expr_matrix, weights, lod) {
  na_mask <- is.na(expr_matrix)
  expr_imputed <- expr_matrix
  expr_imputed[na_mask] <- lod
  expr_imputed <- pmax(expr_imputed, lod)

  weights_imputed <- weights
  if (is.matrix(weights_imputed)) {
    weights_imputed[na_mask] <- 1
  }
  list(
    expression = expr_imputed,
    weights = weights_imputed
  )
}

# Shared body of the build_model_limma* builders. `fit_fun(expression, design, weights = )` fits the
# expression matrix. With `impute = TRUE`, proteins with NA coefficients are refit on LOD-imputed data
# and get the variance borrowed from the successful fits.
.build_model_limma_core <- function(
  lfqdata,
  strategy,
  model_name,
  fit_fun,
  impute = FALSE,
  lod = NULL,
  df_method = "observed"
) {
  setup <- .lfqdata_to_elist(lfqdata, strategy$formula)
  weights <- .resolve_weights(lfqdata, strategy, setup$annotation)
  fit <- fit_fun(setup$expr_matrix, setup$design, weights = weights)
  failed <- if (impute) which(rowSums(is.na(fit$coefficients)) > 0) else integer(0)
  if (length(failed) == 0) {
    return(.new_model_limma(fit, setup, strategy, model_name))
  }

  lod <- .resolve_lod(lod, lfqdata)
  borrowed <- compute_borrowed_variance_limma(fit)
  imputed <- .impute_limma_data(setup$expr_matrix, weights, lod)
  fit_lod <- fit_fun(imputed$expression, setup$design, weights = imputed$weights)
  fit$coefficients[failed, ] <- fit_lod$coefficients[failed, ]
  fit$stdev.unscaled[failed, ] <- fit_lod$stdev.unscaled[failed, ]
  fit$sigma[failed] <- borrowed$sigma
  fit$Amean[failed] <- fit_lod$Amean[failed]
  if (df_method == "observed") {
    fit$df.residual[failed] <- pmax(rowSums(!is.na(setup$expr_matrix))[failed] - ncol(setup$design), 1)
  } else {
    fit$df.residual[failed] <- borrowed$df
  }
  .new_model_limma(
    fit,
    setup,
    strategy,
    model_name,
    dummy_model = .limma_dummy_model(setup$annotation, setup$rhs_formula, imputed$expression[1, ]),
    imputed_proteins = setup$rowdata[failed, setup$subject_id[1], drop = TRUE]
  )
}

.combine_vooma_weights <- function(vooma_weights, external_weights) {
  if (is.null(external_weights)) {
    return(vooma_weights)
  }
  if (is.matrix(external_weights)) {
    return(vooma_weights * external_weights)
  }
  vooma_weights *
    rep(
      external_weights,
      each = nrow(vooma_weights)
    )
}

.plot_vooma_trend <- function(mean_expression, sigma, lowess_trend) {
  graphics::plot(
    mean_expression,
    sigma,
    xlab = "Average log2 expression",
    ylab = expression(sqrt(sigma)),
    main = "vooma: Mean-variance trend",
    pch = 16,
    cex = 0.3
  )
  graphics::lines(lowess_trend, col = "red", lwd = 2)
}

.fit_vooma <- function(
  expression,
  design,
  weights,
  span,
  plot
) {
  preliminary_fit <- limma::lmFit(
    expression,
    design,
    weights = weights
  )
  mean_expression <- preliminary_fit$Amean
  sigma <- sqrt(preliminary_fit$sigma)
  lowess_trend <- stats::lowess(
    mean_expression,
    sigma,
    f = span
  )
  trend_function <- stats::approxfun(
    lowess_trend,
    rule = 2,
    ties = list("ordered", mean)
  )
  fitted_values <- preliminary_fit$coefficients %*%
    t(preliminary_fit$design)
  vooma_weights <- 1 / trend_function(fitted_values)^4
  dim(vooma_weights) <- dim(fitted_values)

  if (plot) {
    .plot_vooma_trend(
      mean_expression,
      sigma,
      lowess_trend
    )
  }
  limma::lmFit(
    expression,
    design,
    weights = .combine_vooma_weights(
      vooma_weights,
      weights
    )
  )
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
#' @param model_name name of model (default from strategy)
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
build_model_limma <- function(lfqdata, strategy, model_name = strategy$model_name) {
  .build_model_limma_core(lfqdata, strategy, model_name, limma::lmFit)
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
#' @param model_name name of model (default: strategy name + "Imputed")
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
  model_name = paste0(strategy$model_name, "Imputed"),
  lod = NULL,
  df_method = c("observed", "borrowed")
) {
  df_method <- match.arg(df_method)
  .build_model_limma_core(lfqdata, strategy, model_name, limma::lmFit, impute = TRUE, lod = lod, df_method = df_method)
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
#' \code{../TODO/prolfqua/Archive/TODO_limma_voom_integration.md} for the mathematical basis.
#'
#' @param lfqdata an \code{\link{LFQData}} object (aggregated to protein level)
#' @param strategy output of \code{\link{strategy_limma}}
#' @param model_name name of model (default from strategy)
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
  model_name = strategy$model_name,
  span = 0.5,
  plot = FALSE
) {
  fit_vooma <- function(expression, design, weights) .fit_vooma(expression, design, weights, span, plot)
  .build_model_limma_core(lfqdata, strategy, model_name, fit_vooma)
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
#' @param model_name name of model (default: strategy name + "Imputed")
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
  model_name = paste0(strategy$model_name, "Imputed"),
  lod = NULL,
  df_method = c("observed", "borrowed"),
  span = 0.5,
  plot = FALSE
) {
  df_method <- match.arg(df_method)
  fit_vooma <- function(expression, design, weights) .fit_vooma(expression, design, weights, span, plot)
  .build_model_limma_core(lfqdata, strategy, model_name, fit_vooma, impute = TRUE, lod = lod, df_method = df_method)
}


# ModelLimma -----

#' R6 class representing a limma modelling result
#'
#' Same API as \code{\link{Model}}: \code{get_anova()}, \code{get_coefficients()},
#' \code{coef_histogram()}, \code{coef_volcano()}, \code{anova_histogram()}.
#'
#' @return An R6 class generator.
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
    #' @field subject_id protein ID column name(s)
    subject_id = character(),
    #' @field model_name model name
    model_name = character(),
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
    #' @field imputed_proteins character vector of protein ids (matching
    #'   \code{rownames(fit$coefficients)}) that were rescued via LOD
    #'   imputation. \code{ContrastsLimma$get_contrasts()} reads this to
    #'   tag refit rows with a \code{"_imputed"} modelName suffix so the
    #'   rescue is visible downstream. Empty for builders that do not
    #'   impute.
    imputed_proteins = character(0),
    #' @description
    #' initialize ModelLimma
    #' @param fit limma MArrayLM from lmFit
    #' @param subject_id protein ID column name(s)
    #' @param model_name model name
    #' @param rowdata protein ID mapping
    #' @param trend passed to eBayes
    #' @param robust passed to eBayes
    #' @param dummy_model one fitted lm for linfct extraction
    #' @param p.adjust function to adjust p-values
    #' @param imputed_proteins character vector of LOD-rescued protein ids
    initialize = function(
      fit,
      subject_id,
      model_name,
      rowdata,
      trend = FALSE,
      robust = FALSE,
      dummy_model = NULL,
      p.adjust = prolfqua::adjust_p_values,
      imputed_proteins = character(0)
    ) {
      self$fit <- fit
      self$subject_id <- subject_id
      self$model_name <- model_name
      self$rowdata <- rowdata
      self$trend <- trend
      self$robust <- robust
      self$dummy_model <- dummy_model
      self$p.adjust <- p.adjust
      self$imputed_proteins <- imputed_proteins
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
    }
  ),
  private = list(
    volcano_colour = NULL,
    volcano_prefix = "Coef_VolcanoPlot_"
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
#' @return An R6 class generator.
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
    #' @field model_name model name
    model_name = character(),
    #' @field subject_id columns with subject_id (proteinID)
    subject_id = character(),
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
    #' @param model_name name of the contrast method
    #' @param eBayes logical, apply limma eBayes moderation (default TRUE).
    #'   Set to FALSE to return raw/unmoderated statistics suitable for
    #'   wrapping with \code{\link{ContrastsModeratedDEqMS}}.
    initialize = function(model, contrasts, p.adjust = prolfqua::adjust_p_values, model_name = NULL, eBayes = TRUE) {
      self$model <- model
      self$contrasts <- contrasts
      self$subject_id <- model$subject_id
      self$eBayes <- eBayes
      self$model_name <- model_name %||% if (eBayes) "limma" else "limma_raw"
      self$p.adjust <- p.adjust
    },
    #' @description
    #' get both sides of contrasts
    get_contrast_sides = function() {
      parse_contrast_sides(self$contrasts)
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

      # linfct_a: rows = contrasts and their avg_ counterparts (for avgAbd), cols = model coefficients
      linfct_a <- self$get_linfct()
      diff_names <- names(self$contrasts)
      avg_names <- paste0("avg_", diff_names)

      # Transpose for limma: rows = coefficients, cols = contrasts
      contrast_matrix <- t(linfct_a[diff_names, , drop = FALSE])

      # limma pipeline: contrasts.fit, optionally + eBayes
      fit2 <- limma::contrasts.fit(self$model$fit, contrast_matrix)

      if (self$eBayes) {
        fit2 <- limma::eBayes(fit2, trend = self$model$trend, robust = self$model$robust)
      }

      # avgAbd: the avg_ linear functions applied to the fit coefficients (rows = proteins)
      avg_vals <- self$model$fit$coefficients %*% t(linfct_a[avg_names, , drop = FALSE])

      # Extract results per contrast
      res_list <- vector("list", length(diff_names))
      for (i in seq_along(diff_names)) {
        df_i <- self$model$rowdata
        df_i$contrast <- diff_names[i]

        if (self$eBayes) {
          tt <- limma::topTable(fit2, coef = i, number = Inf, sort.by = "none", confint = TRUE)
          df_i$diff <- tt$logFC
          df_i$std.error.unmoderated <- fit2$sigma * fit2$stdev.unscaled[, i]
          df_i$df.unmoderated <- fit2$df.residual
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
          df_i$std.error.unmoderated <- se_i
          df_i$statistic <- t_raw
          df_i$p.value <- p_raw
          df_i$sigma <- sigma_raw
          df_i$df <- df_raw
          df_i$df.unmoderated <- df_raw
          df_i$conf.low <- diff_i - prqt * se_i
          df_i$conf.high <- diff_i + prqt * se_i
        }
        df_i$avgAbd <- avg_vals[, i]
        res_list[[i]] <- df_i
      }
      contrast_result <- dplyr::bind_rows(res_list)

      # Adjust p-values per contrast
      contrast_result <- self$p.adjust(contrast_result, column = "p.value", group_by_col = "contrast", newname = "FDR")
      contrast_result <- contrast_result |> dplyr::relocate("FDR", .after = "diff")
      # Stamp modelName uniformly with the model identity. Refit proteins
      # (recorded as bare subject_id values on ModelLimma$imputed_proteins by
      # build_model_limma_impute and build_model_limma_voom_impute) are flagged
      # in the separate `estimate_type` column with "lod_imputed"; everything
      # else is "observed". Match on the primary subject_id column; no current
      # limma backend uses multi-key subject_ids.
      is_imputed <- contrast_result[[self$subject_id[1]]] %in% self$model$imputed_proteins
      contrast_result$estimate_type <- dplyr::if_else(is_imputed, "lod_imputed", "observed")
      contrast_result <- dplyr::mutate(contrast_result, modelName = self$model_name, .before = 1)
      contrast_result <- dplyr::relocate(contrast_result, "estimate_type", .after = "modelName")
      contrast_result <- dplyr::ungroup(contrast_result)
      self$contrast_result <- contrast_result

      stopifnot(all(super$column_description()$column_name %in% colnames(contrast_result)))
      return(contrast_result)
    }
  )
)
