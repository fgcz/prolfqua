# ContrastsChildToParentFacades -----
#
# Facades that take child-hierarchy (peptide/precursor) LFQData and return
# parent-level (protein) fold-change estimates. The hierarchy roll-up happens
# either inside the modelling backend (Lmer, ROPECA, Firth) or by an explicit
# pre-aggregation step driven by the facade itself (Limpa).
#
# Same-level (protein -> protein) facades live in R/ContrastsFacades.R. Shared
# helpers .assert_nested_facade_input(), .compute_missing(), and
# .stamp_facade_identity() are defined there and used here as package-internal
# functions.

# If the caller-supplied modelstr has no random-effect bars, append
# (1 | <deepest child hierarchy key>) + (1 | <sample name>) derived from the
# LFQData. Caller-supplied random effects are left untouched.
.augment_lmer_nested_modelstr <- function(modelstr, lfqdata) {
  if (grepl("\\|", modelstr)) {
    return(modelstr)
  }
  child_keys <- setdiff(lfqdata$hierarchy_keys(), lfqdata$subject_id())
  if (length(child_keys) == 0) {
    stop(
      "ContrastsLmerNestedFacade: cannot derive random effects - ",
      "no child hierarchy below subject_id().",
      call. = FALSE
    )
  }
  deepest_child <- child_keys[length(child_keys)]
  paste0(modelstr, " + (1 | ", deepest_child, ") + (1 | ", lfqdata$sample_name(), ")")
}

#' Lmer (mixed-model) contrast analysis facade for nested input
#'
#' Encapsulates the pipeline: \code{\link{strategy_lmer}} ->
#' \code{\link{build_model}} -> \code{\link{Contrasts}} ->
#' \code{\link{ContrastsModerated}}.
#'
#' Takes child-hierarchy LFQData (e.g. peptide-level measurements nested
#' within proteins) and returns protein-level fold-change estimates.
#'
#' \code{modelstr} may be either a fixed-effects-only formula (e.g.
#' \code{"~ group_"}) - in which case the facade appends
#' \code{(1 | <deepest child hierarchy key>) + (1 | <sample name>)} derived
#' from the LFQData - or a full mixed-model formula already containing random
#' effects, which is passed through unchanged.
#'
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_peptide_config()
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' # Random effects derived from lfqdata; equivalent to passing the full
#' # mixed-model formula explicitly.
#' fa <- ContrastsLmerNestedFacade$new(lfqdata, "~ group_", contrasts)
#' head(fa$get_contrasts())
#' fa$to_wide()
#' @return An R6 class generator.
ContrastsLmerNestedFacade <- R6::R6Class(
  "ContrastsLmerNestedFacade",
  inherit = ContrastsFacadeBase,
  public = list(
    #' @field model Model object
    model = NULL,
    #' @field contrast ContrastsModerated object
    contrast = NULL,
    #' @field .lfqdata stored reference to input LFQData
    .lfqdata = NULL,
    #' @field .contrast_names names of the requested contrasts
    .contrast_names = NULL,
    #' @description
    #' initialize
    #' @param lfqdata nested LFQData (subject_id strict subset of hierarchy_keys)
    #' @param modelstr model formula RHS as string. Fixed effects only
    #'   (e.g. "~ group_") - random effects are derived from \code{lfqdata}
    #'   as \code{(1 | <deepest child key>) + (1 | <sample name>)}. A
    #'   formula already containing random-effect bars is passed through
    #'   unchanged.
    #' @param contrasts named character vector of contrasts
    #' @param ... passed to \code{\link{strategy_lmer}}
    initialize = function(lfqdata, modelstr, contrasts, ...) {
      .assert_nested_facade_input(lfqdata, "ContrastsLmerNestedFacade")
      self$.lfqdata <- lfqdata
      self$.contrast_names <- names(contrasts)
      response <- lfqdata$response()
      modelstr <- .augment_lmer_nested_modelstr(modelstr, lfqdata)
      full_formula <- paste(response, modelstr)
      strat <- strategy_lmer(full_formula, ...)
      self$model <- build_model(lfqdata, strat)
      self$contrast <- ContrastsModerated$new(Contrasts$new(self$model, contrasts, model_name = "lmer_nested"))
      self$facade_name <- "lmer_nested"
      self$config <- self$contrast$get_config()
    }
  )
)


#' ROPECA contrast analysis facade for nested input
#'
#' Encapsulates the pipeline: \code{\link{strategy_lm}} ->
#' \code{\link{build_model}} -> \code{\link{Contrasts}} ->
#' \code{\link{ContrastsROPECA}}.
#'
#' ROPECA fits per-peptide linear models and aggregates peptide-level
#' p-values to the protein level. Input must be nested LFQData (peptide
#' subject_id strictly below the protein hierarchy).
#'
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_peptide_config()
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' fa <- ContrastsROPECANestedFacade$new(lfqdata, "~ group_", contrasts)
#' head(fa$get_contrasts())
#' fa$to_wide()
#' @return An R6 class generator.
ContrastsROPECANestedFacade <- R6::R6Class(
  "ContrastsROPECANestedFacade",
  inherit = ContrastsFacadeBase,
  public = list(
    #' @field model Model object (peptide-level)
    model = NULL,
    #' @field contrast ContrastsROPECA object
    contrast = NULL,
    #' @field .lfqdata stored reference to input LFQData
    .lfqdata = NULL,
    #' @field .contrast_names names of the requested contrasts
    .contrast_names = NULL,
    #' @description
    #' initialize
    #' @param lfqdata nested LFQData (peptide-level)
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param ... passed to \code{\link{strategy_lm}}
    initialize = function(lfqdata, modelstr, contrasts, ...) {
      .assert_nested_facade_input(lfqdata, "ContrastsROPECANestedFacade")
      self$.lfqdata <- lfqdata
      self$.contrast_names <- names(contrasts)
      response <- lfqdata$response()
      full_formula <- paste(response, modelstr)
      strat <- strategy_lm(full_formula, ...)
      subject_id <- lfqdata$hierarchy_keys()
      self$model <- build_model(lfqdata, strat, subject_id = subject_id)
      self$contrast <- ContrastsROPECA$new(Contrasts$new(self$model, contrasts))
      self$config <- self$contrast$get_config()
      self$config$subject_id <- self$contrast$subject_id[1]
      self$facade_name <- "ropeca_nested"
    },
    #' @description
    #' Get contrast results with standardized column names.
    #'
    #' ROPECA's \code{beta.based.significance} is mapped to \code{p.value} and
    #' \code{FDR.beta.based.significance} to \code{FDR}.
    #'
    #' Columns not directly produced by ROPECA are derived heuristically:
    #' \itemize{
    #'   \item \code{std.error = diff / statistic} (algebraic: t = estimate / SE)
    #'   \item \code{sigma = mad.estimate} (MAD of peptide fold changes)
    #'   \item \code{df = n_not_na} (number of contributing peptides)
    #'   \item \code{conf.low/conf.high} via \code{diff +/- qt(0.975, df) * |std.error|}
    #' }
    get_contrasts = function() {
      res <- self$contrast$get_contrasts(all = TRUE)
      res <- dplyr::rename(res, p.value = "beta.based.significance", FDR = "FDR.beta.based.significance")
      res$std.error <- ifelse(res$statistic != 0, res$diff / res$statistic, NA_real_)
      res$sigma <- res$mad.estimate
      res$df <- res$n_not_na
      prqt <- stats::qt(0.975, df = pmax(res$df, 1))
      res$conf.low <- res$diff - prqt * abs(res$std.error)
      res$conf.high <- res$diff + prqt * abs(res$std.error)
      protein_Id <- self$contrast$subject_id[1]
      standard_cols <- c(
        protein_Id,
        "modelName",
        "estimate_type",
        "contrast",
        "avgAbd",
        "diff",
        "FDR",
        "statistic",
        "std.error",
        "df",
        "p.value",
        "conf.low",
        "conf.high",
        "sigma"
      )
      res <- res[, standard_cols, drop = FALSE]
      .stamp_facade_identity(res, "ropeca_nested")
    },
    #' @description get ContrastsPlotter (uses standardized column names)
    #' @param fc_threshold fold change threshold
    #' @param fdr_threshold FDR threshold
    get_Plotter = function(fc_threshold = 2, fdr_threshold = 0.1) {
      contrast_result <- self$get_contrasts()
      protein_Id <- self$contrast$subject_id[1]
      ContrastsPlotter$new(
        contrast_result,
        subject_id = protein_Id,
        fcthresh = fc_threshold,
        volcano = list(list(score = "FDR", thresh = fdr_threshold)),
        histogram = list(list(score = "p.value", xlim = c(0, 1, 0.05)), list(score = "FDR", xlim = c(0, 1, 0.05))),
        modelName = "modelName",
        diff = "diff",
        contrast = "contrast"
      )
    },
    #' @description convert results to wide format
    #' @param columns value columns to pivot, default \code{c("p.value", "FDR", "statistic")}
    to_wide = function(columns = c("p.value", "FDR", "statistic")) {
      contrast_result <- self$get_contrasts()
      protein_Id <- self$contrast$subject_id[1]
      pivot_model_contrasts_to_wide(
        contrast_result,
        subject_id = protein_Id,
        columns = c("diff", columns),
        contrast = "contrast"
      )
    }
  )
)


#' Firth logistic missingness contrast analysis facade for nested input
#'
#' Encapsulates the pipeline: encode missingness ->
#' \code{\link{build_model_glm_peptide}} -> \code{\link{ContrastsFirth}}.
#'
#' Takes nested (peptide-level) LFQData and returns protein-level
#' fold-change estimates. For protein-level (aggregated) input use
#' \code{\link{ContrastsFirthFacade}} instead.
#'
#' Supports \code{options(prolfqua.vectorize = TRUE)} for faster contrast
#' computation. See \code{\link{build_contrast_analysis}} for details.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_peptide_config(Nprot = 20, weight_missing = 0.5)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' fa <- ContrastsFirthNestedFacade$new(lfqdata, "~ group_", contrasts)
#' head(fa$get_contrasts())
#' fa$to_wide()
ContrastsFirthNestedFacade <- R6::R6Class(
  "ContrastsFirthNestedFacade",
  inherit = ContrastsFacadeBase,
  public = list(
    #' @field model ModelFirth object
    model = NULL,
    #' @field contrast ContrastsFirth object
    contrast = NULL,
    #' @field .lfqdata stored reference to input LFQData
    .lfqdata = NULL,
    #' @field .contrast_names names of the requested contrasts
    .contrast_names = NULL,
    #' @description
    #' initialize
    #' @param lfqdata nested LFQData (peptide-level)
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param ... ignored; accepted for dispatch compatibility with other facades
    initialize = function(lfqdata, modelstr, contrasts, ...) {
      .assert_nested_facade_input(lfqdata, "ContrastsFirthNestedFacade")
      self$.lfqdata <- lfqdata
      self$.contrast_names <- names(contrasts)
      self$model <- build_model_glm_peptide(lfqdata, modelstr)
      self$contrast <- ContrastsFirth$new(self$model, contrasts, model_name = "firth_nested")
      self$facade_name <- "firth_nested"
      self$config <- self$contrast$get_config()
    }
  )
)


#' Limpa contrast analysis facade for nested input
#'
#' Encapsulates the full precursor -> protein pipeline:
#' \code{\link{AggregateLimpa}} -> \code{\link{strategy_limpa}} ->
#' \code{\link{build_model_limpa}} -> \code{\link{ContrastsLimma}}.
#'
#' Takes nested (precursor/peptide-level) LFQData, runs limpa's
#' DPC-based aggregation internally to produce protein-level expression
#' with posterior standard errors, then fits a vooma model with
#' imputation-aware precision weights.
#'
#' For protein-level input that was already aggregated upstream via
#' \code{\link{AggregateLimpa}}, use \code{\link{ContrastsLimpaFacade}}
#' instead.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' if (requireNamespace("limpa", quietly = TRUE)) {
#'   istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 10)
#'   lfqdata <- LFQData$new(istar$data, istar$config)
#'   lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#'   contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#'   fa <- ContrastsLimpaNestedFacade$new(lfqdata, "~ group_", contrasts)
#'   head(fa$get_contrasts())
#' }
ContrastsLimpaNestedFacade <- R6::R6Class(
  "ContrastsLimpaNestedFacade",
  inherit = ContrastsFacadeBase,
  public = list(
    #' @field model ModelLimma object (from build_model_limpa)
    model = NULL,
    #' @field contrast ContrastsLimma object
    contrast = NULL,
    #' @field .lfqdata stored reference to the aggregated protein-level LFQData
    .lfqdata = NULL,
    #' @field .lfqdata_nested stored reference to the original nested input
    .lfqdata_nested = NULL,
    #' @field .contrast_names names of the requested contrasts
    .contrast_names = NULL,
    #' @description
    #' initialize
    #' @param lfqdata nested LFQData (precursor/peptide-level, log2-transformed)
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param prefix prefix for the aggregated hierarchy level (default "protein")
    #' @param dpc_slope DPC slope parameter passed to \code{\link{AggregateLimpa}}
    #'   (default 0.8)
    #' @param plot logical; if TRUE, plot the vooma mean-variance trend
    #' @param span lowess smoother span (NULL = auto)
    #' @param ... passed to \code{\link{strategy_limpa}} (e.g. trend, robust)
    initialize = function(
      lfqdata,
      modelstr,
      contrasts,
      prefix = "protein",
      dpc_slope = 0.8,
      plot = FALSE,
      span = NULL,
      ...
    ) {
      .assert_nested_facade_input(lfqdata, "ContrastsLimpaNestedFacade")
      self$.lfqdata_nested <- lfqdata
      self$.contrast_names <- names(contrasts)
      lfq_agg <- AggregateLimpa$new(lfqdata, prefix = prefix, dpc_slope = dpc_slope)$aggregate()
      self$.lfqdata <- lfq_agg
      response <- lfq_agg$response()
      full_formula <- paste(response, modelstr)
      strat <- strategy_limpa(full_formula, plot = plot, span = span, ...)
      self$model <- build_model_limpa(lfq_agg, strat)
      self$contrast <- ContrastsLimma$new(self$model, contrasts, model_name = "limpa_nested")
      self$facade_name <- "limpa_nested"
      self$.drop_na_diff <- TRUE
      self$config <- self$contrast$get_config()
    }
  )
)
