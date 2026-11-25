# ContrastsChildToParentFacades -----
#
# Facades that take child-hierarchy (peptide/precursor) LFQData and return
# parent-level (protein) fold-change estimates. The hierarchy roll-up happens
# either inside the modelling backend (Lmer, ROPECA, Firth) or by an explicit
# pre-aggregation step driven by the facade itself (Limpa).
#
# Same-level (protein -> protein) facades live in R/ContrastsFacades.R. The
# shared helper .stamp_facade_identity() is defined there and used here as a
# package-internal function.

#' If the caller-supplied modelstr has no random-effect bars, append
#'
#' (1 | <deepest child hierarchy key>) + (1 | <sample name>) derived from the
#' LFQData. Caller-supplied random effects are left untouched.
#' @param modelstr right-hand-side model formula string
#' @param lfqdata [LFQData] object with the data and its configuration
#' @noRd
.augment_lmer_nested_modelstr <- function(modelstr, lfqdata) {
  if (grepl("|", modelstr, fixed = TRUE)) {
    return(modelstr)
  }
  child_keys <- setdiff(lfqdata$hierarchy_keys(), lfqdata$subject_id())
  if (length(child_keys) == 0) {
    abort_prolfqua(
      "Can't derive random effects: there is no hierarchy level below `subject_id()`."
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
#' @examplesIf requireNamespace("lme4", quietly = TRUE) && requireNamespace("lmerTest", quietly = TRUE)
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
    #' @description
    #' initialize
    #' @param lfqdata nested LFQData (subject_id strict subset of hierarchy_keys)
    #' @param modelstr model formula RHS as string. Fixed effects only
    #'   (e.g. "~ group_") - random effects are derived from \code{lfqdata}
    #'   as \code{(1 | <deepest child key>) + (1 | <sample name>)}. A
    #'   formula already containing random-effect bars is passed through
    #'   unchanged.
    #' @param contrasts named character vector of contrasts
    #' @param control named list of convergence controls, see \code{\link{strategy_lmer}}
    #' @param weights column names whose product gives the per-observation weights (default: the identification
    #'   confidence). Pass \code{NULL} for unweighted.
    initialize = function(lfqdata, modelstr, contrasts, control = list(), weights = .confidence_weights(lfqdata)) {
      private$.setup(lfqdata, contrasts, "lmer_nested", modelstr)
      full_formula <- paste(lfqdata$response(), .augment_lmer_nested_modelstr(modelstr, lfqdata))
      strat <- strategy_lmer(full_formula, control = control, weights = weights)
      self$model <- build_model(lfqdata, strat)
      self$contrast <- ContrastsModerated$new(Contrasts$new(self$model, contrasts, model_name = "lmer_nested"))
    }
  ),
  private = list(
    .needs = "nested"
  )
)


#' ROPECA contrast analysis facade for nested input
#'
#' Encapsulates the pipeline: \code{\link{strategy_lm}} ->
#' \code{\link{build_model}} -> \code{\link{Contrasts}} ->
#' \code{\link{ContrastsROPECA}}.
#'
#' Input must be nested LFQData (peptide subject_id strictly below the protein hierarchy).
#'
#' `ropeca_nested` is a derivation of PECA and ROPECA (Suomi and Elo 2017, Scientific Reports 7, 5869; Bioconductor
#' package PECA, \doi{10.18129/B9.bioc.PECA}). It fits one linear model per peptide, computes the contrast, takes
#' the median of the signed peptide p-values per protein and corrects it with the beta distribution of the median
#' of `n` uniform p-values, `n` being the number of peptides with a p-value, at most 10. The protein fold change is
#' the median peptide fold change. What differs from PECA and ROPECA:
#' \itemize{
#'   \item no ROTS (reproducibility-optimized test statistic) step: peptide p-values come from ordinary linear
#'     models;
#'   \item any fixed-effect design and contrast, not only two groups;
#'   \item `n` is capped at 10, so proteins with more peptides get larger p-values than in PECA;
#'   \item with an even number of peptides, the median is taken of the signed p-values, not of the t-statistics.
#' }
#' The test is conservative under the null hypothesis; `ropecaV2_nested` is the calibrated variant.
#'
#' @section Comparison with PECA:
#' On the IonStar data (`data_ionstar`, dilutions a and e, 767 peptides observed in all 8 samples, 123 proteins),
#' `ropeca_nested` and `PECA::PECA_df(test = "t", type = "median")` give the same fold changes (largest difference
#' 1e-14) and the same p-values for the 58 proteins with an odd number of peptides up to 10 (largest difference
#' 3e-14). For the 50 proteins with an even number, the median ratio of the p-values is 1.02 (Spearman correlation
#' 0.99); for the 15 proteins with more than 10 peptides, prolfqua's p-values are larger by a median factor of 1.43
#' (Spearman correlation 0.87), because of the cap. Code: `data-raw/compare_peca.R` in the source repository.
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
    #' @description
    #' initialize
    #' @param lfqdata nested LFQData (peptide-level)
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param weights column names whose product weights the per-peptide fits (default: the identification
    #'   confidence). Unless \code{NULL}, the peptides are also combined with weights
    #'   \code{identification_confidence_feature}. Pass \code{NULL} for unweighted.
    initialize = function(lfqdata, modelstr, contrasts, weights = .confidence_weights(lfqdata)) {
      full_formula <- private$.setup(lfqdata, contrasts, "ropeca_nested", modelstr)
      strat <- strategy_lm(full_formula, weights = weights)
      subject_id <- lfqdata$hierarchy_keys()
      self$model <- build_model(lfqdata, strat, subject_id = subject_id)
      self$contrast <- ContrastsROPECA$new(
        Contrasts$new(self$model, contrasts),
        child_weights = if (is.null(weights)) NULL else .child_weights(lfqdata)
      )
      self$config <- self$contrast$get_config()
      self$config$subject_id <- self$contrast$subject_id[1]
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
    #' @param all \code{FALSE} (default) returns the standard columns;
    #'   \code{TRUE} appends the remaining ROPECA columns
    get_contrasts = function(all = FALSE) {
      res <- .ropeca_standard_columns(
        self$contrast$get_contrasts(all = TRUE),
        self$contrast$subject_id[1],
        all = all
      )
      .stamp_facade_identity(res, self$facade_name)
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
        volcano = list(list(score = "FDR", thresh = fdr_threshold))
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
  ),
  private = list(
    .needs = "nested"
  )
)


#' ROPECA V2 contrast analysis facade for nested input
#'
#' Encapsulates the pipeline: \code{\link{strategy_lm}} ->
#' \code{\link{build_model}} -> \code{\link{Contrasts}} ->
#' \code{\link{ContrastsModerated}} -> \code{\link{ContrastsROPECAV2}}.
#'
#' Like \code{\link{ContrastsROPECANestedFacade}}, it fits one linear model per
#' peptide and reports the median peptide fold change per protein. The protein
#' p-value tests the median of the moderated peptide z-scores against its exact
#' null distribution, allowing for the correlation of the peptides of a protein,
#' which is estimated once from the peptide model residuals. Input must be nested
#' LFQData (peptide subject_id strictly below the protein hierarchy).
#'
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_peptide_config()
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' fa <- ContrastsROPECAV2NestedFacade$new(lfqdata, "~ group_", contrasts)
#' head(fa$get_contrasts())
#' fa$contrast$rho
#' fa$to_wide()
#' @return An R6 class generator.
ContrastsROPECAV2NestedFacade <- R6::R6Class(
  "ContrastsROPECAV2NestedFacade",
  inherit = ContrastsFacadeBase,
  public = list(
    #' @description
    #' initialize
    #' @param lfqdata nested LFQData (peptide-level)
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param weights column names whose product weights the per-peptide fits (default: the identification
    #'   confidence). Unless \code{NULL}, the peptides are also combined with weights
    #'   \code{identification_confidence_feature}. Pass \code{NULL} for unweighted.
    initialize = function(lfqdata, modelstr, contrasts, weights = .confidence_weights(lfqdata)) {
      full_formula <- private$.setup(lfqdata, contrasts, "ropecaV2_nested", modelstr)
      strat <- strategy_lm(full_formula, weights = weights)
      subject_id <- lfqdata$hierarchy_keys()
      self$model <- build_model(lfqdata, strat, subject_id = subject_id)
      rho <- .ropeca_v2_rho(self$model$model_df, subject_id, lfqdata$sample_name())
      peptide_contrasts <- ContrastsModerated$new(Contrasts$new(self$model, contrasts))
      self$contrast <- ContrastsROPECAV2$new(
        peptide_contrasts,
        rho = rho,
        child_weights = if (is.null(weights)) NULL else .child_weights(lfqdata)
      )
      self$config <- self$contrast$get_config()
      self$config$subject_id <- self$contrast$subject_id[1]
    },
    #' @description
    #' Get protein-level contrasts with the standard columns.
    #' @param all \code{FALSE} (default) returns the standard columns;
    #'   \code{TRUE} appends \code{n_peptides}, \code{median_z} and \code{rho}
    get_contrasts = function(all = FALSE) {
      .stamp_facade_identity(self$contrast$get_contrasts(all = all), self$facade_name)
    },
    #' @description get ContrastsPlotter (uses standardized column names)
    #' @param fc_threshold fold change threshold
    #' @param fdr_threshold FDR threshold
    get_Plotter = function(fc_threshold = 2, fdr_threshold = 0.1) {
      ContrastsPlotter$new(
        self$get_contrasts(),
        subject_id = self$contrast$subject_id[1],
        fcthresh = fc_threshold,
        volcano = list(list(score = "FDR", thresh = fdr_threshold))
      )
    },
    #' @description convert results to wide format
    #' @param columns value columns to pivot, default \code{c("p.value", "FDR", "statistic")}
    to_wide = function(columns = c("p.value", "FDR", "statistic")) {
      pivot_model_contrasts_to_wide(
        self$get_contrasts(),
        subject_id = self$contrast$subject_id[1],
        columns = c("diff", columns),
        contrast = "contrast"
      )
    }
  ),
  private = list(
    .needs = "nested"
  )
)


#' Quasibinomial detection-count facade for nested input
#'
#' Completes the child-feature rows, encodes detection as 0/1, then collapses
#' the binary rows into
#' detected and undetected counts per parent and sample. The resulting
#' quasibinomial model reports parent-level log odds ratios in \code{diff};
#' \code{avgAbd} is the average linear predictor on the log-odds scale.
#'
#' Child features are treated as exchangeable binomial trials. The symmetric
#' pseudo-count stabilizes complete separation but is not equivalent to Firth's
#' bias-reducing penalty. Empirical-Bayes dispersion moderation uses
#' \code{\link{ContrastsModerated}}; by default posterior dispersion is bounded
#' below by one.
#'
#' @return An R6 class generator.
#' @export
#' @examples
#' istar <- sim_lfq_data_peptide_config(n_proteins = 20, weight_missing = 0.5, seed = 3)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' facade <- ContrastsBinomialNestedFacade$new(lfqdata, "~ group_", contrasts)
#' head(facade$get_contrasts())
ContrastsBinomialNestedFacade <- R6::R6Class(
  "ContrastsBinomialNestedFacade",
  inherit = ContrastsFacadeBase,
  public = list(
    #' @description Fit the nested detection-count model and its contrasts.
    #' @param lfqdata nested \code{\link{LFQData}}
    #' @param modelstr right-hand-side model formula
    #' @param contrasts named contrast expressions
    #' @param prior_count non-negative symmetric pseudo-count
    #' @param binomial_bound bound posterior dispersion below by one
    #' @param control named list of convergence controls, see \code{\link{strategy_binomial}}
    initialize = function(
      lfqdata,
      modelstr,
      contrasts,
      prior_count = 0.1,
      binomial_bound = TRUE,
      control = list()
    ) {
      private$.setup(lfqdata, contrasts, "binomial_nested", modelstr)
      strategy <- strategy_binomial(modelstr, prior_count = prior_count, control = control)
      fitted <- .fit_detection_counts(lfqdata, strategy, contrasts, binomial_bound)
      self$model <- fitted$model
      self$contrast <- fitted$contrast
    }
  ),
  private = list(
    .needs = "nested"
  )
)


#' Firth-penalized detection-count facade for nested input
#'
#' Same pipeline as \code{\link{ContrastsBinomialNestedFacade}}: child-feature
#' detection is collapsed into detected and undetected counts per parent and
#' sample, so the dispersion and the residual degrees of freedom come from the
#' samples, not from the number of child features. The fixed pseudo-count is
#' replaced by Firth's bias reduction (\code{\link{StrategyBinomialFirth}}),
#' which gives finite log odds ratios under complete separation. Child
#' identity is not modelled.
#'
#' @return An R6 class generator.
#' @export
#' @examples
#' istar <- sim_lfq_data_peptide_config(n_proteins = 20, weight_missing = 0.5, seed = 3)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' facade <- ContrastsBinomialFirthNestedFacade$new(lfqdata, "~ group_", contrasts)
#' head(facade$get_contrasts())
ContrastsBinomialFirthNestedFacade <- R6::R6Class(
  "ContrastsBinomialFirthNestedFacade",
  inherit = ContrastsFacadeBase,
  public = list(
    #' @description Fit the Firth-penalized nested detection-count model and its contrasts.
    #' @param lfqdata nested \code{\link{LFQData}}
    #' @param modelstr right-hand-side model formula
    #' @param contrasts named contrast expressions
    #' @param binomial_bound bound posterior dispersion below by one
    #' @param control named list of convergence controls, see \code{\link{strategy_binomial_firth}}
    initialize = function(lfqdata, modelstr, contrasts, binomial_bound = TRUE, control = list()) {
      private$.setup(lfqdata, contrasts, "binomial_firth_nested", modelstr)
      strategy <- strategy_binomial_firth(modelstr, control = control)
      fitted <- .fit_detection_counts(lfqdata, strategy, contrasts, binomial_bound)
      self$model <- fitted$model
      self$contrast <- fitted$contrast
    }
  ),
  private = list(
    .needs = "nested"
  )
)


#' Shared pipeline of the detection-count facades: collapse child detection into
#'
#' counts per parent and sample, fit `strategy`, and moderate the Wald contrasts.
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param strategy modelling strategy
#' @param contrasts named character vector of contrast expressions
#' @param binomial_bound bound applied to binomial detection estimates
#' @noRd
.fit_detection_counts <- function(lfqdata, strategy, contrasts, binomial_bound) {
  counts <- .summarize_detection_counts(.prepare_detection_lfqdata(lfqdata))
  model <- build_model(
    counts,
    strategy,
    subject_id = lfqdata$subject_id(),
    model_name = strategy$model_name
  )
  base_contrast <- Contrasts$new(model, contrasts, model_name = strategy$model_name)
  base_contrast$config <- ContrastConfiguration$new(
    subject_id = lfqdata$subject_id(),
    supports_dea_qc = FALSE
  )
  variance_floor <- if (isTRUE(binomial_bound)) 1 else NULL
  contrast <- ContrastsModerated$new(
    base_contrast,
    model_name = strategy$model_name,
    variance_floor = variance_floor
  )
  list(model = model, contrast = contrast)
}


#' One weight per child feature: its feature-level identification confidence
#' @param lfqdata nested LFQData
#' @return data frame with the hierarchy key columns and `weight`
#' @noRd
.child_weights <- function(lfqdata) {
  feature_col <- lfqdata$get_config()$identification_confidence_feature
  weights <- dplyr::distinct(lfqdata$data_long()[, c(lfqdata$hierarchy_keys(), feature_col)])
  dplyr::rename(weights, weight = dplyr::all_of(feature_col))
}

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
#'   istar <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 10)
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
    #' @param trend,robust passed to \code{\link[limma]{eBayes}}, see \code{\link{strategy_limpa}}
    #' @param weights column names whose product gives the prior weights of the protein-level fit (default: the
    #'   identification confidence, the maximum over the observed children). Pass \code{NULL} for unweighted.
    initialize = function(
      lfqdata,
      modelstr,
      contrasts,
      prefix = "protein",
      dpc_slope = 0.8,
      plot = FALSE,
      span = NULL,
      trend = FALSE,
      robust = FALSE,
      weights = .confidence_weights(lfqdata)
    ) {
      private$.setup(lfqdata, contrasts, "limpa_nested", modelstr)
      lfq_agg <- AggregateLimpa$new(lfqdata, prefix = prefix, dpc_slope = dpc_slope)$aggregate()
      # get_missing() compares against the aggregated protein-level data.
      self$.lfqdata <- lfq_agg
      strat <- strategy_limpa(
        paste(lfq_agg$response(), modelstr),
        plot = plot,
        span = span,
        trend = trend,
        robust = robust,
        weights = weights
      )
      self$model <- build_model_limpa(lfq_agg, strat)
      self$contrast <- ContrastsLimma$new(self$model, contrasts, model_name = "limpa_nested")
      self$.drop_na_diff <- TRUE
    }
  ),
  private = list(
    .needs = "nested"
  )
)
