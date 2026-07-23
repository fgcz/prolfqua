# ContrastsFacades -----

.assert_aggregated_facade_input <- function(lfqdata, facade_name) {
  subject_id <- lfqdata$subject_id()
  hierarchy_keys <- lfqdata$hierarchy_keys()
  if (!identical(subject_id, hierarchy_keys)) {
    stop(
      facade_name,
      " requires aggregated LFQData. ",
      "`lfqdata$subject_id()` must equal `lfqdata$hierarchy_keys()`. ",
      "Aggregate first.",
      call. = FALSE
    )
  }
}

.assert_nested_facade_input <- function(lfqdata, facade_name) {
  subject_id <- lfqdata$subject_id()
  hierarchy_keys <- lfqdata$hierarchy_keys()
  if (!(all(subject_id %in% hierarchy_keys) && length(subject_id) < length(hierarchy_keys))) {
    stop(
      facade_name,
      " requires LFQData with additional hierarchy below `subject_id()`. ",
      "`lfqdata$subject_id()` must be a strict subset of ",
      "`lfqdata$hierarchy_keys()`. Do not aggregate first.",
      call. = FALSE
    )
  }
}

# Compute protein × contrast pairs present in input but absent from output.
# Returns a data.frame with hierarchy columns + contrast, or 0 rows if nothing missing.
.compute_missing <- function(lfqdata, contrast_names, contrast_result) {
  subject_id <- lfqdata$subject_id()
  all_subjects <- lfqdata$data_long() |>
    dplyr::select(dplyr::all_of(subject_id)) |>
    dplyr::distinct()
  expected <- tidyr::crossing(all_subjects, contrast = contrast_names)
  estimated <- contrast_result |>
    dplyr::select(dplyr::all_of(c(subject_id, "contrast"))) |>
    dplyr::distinct()
  dplyr::anti_join(expected, estimated, by = c(subject_id, "contrast"))
}

# Stamp the selected facade key as the row's model identity. `modelName` is the
# single identity column: the inner contrast classes set it to their own model
# name (e.g. "WaldTest", "limma") and this overwrites it with the facade key
# (e.g. "rfit", "lm_impute") that the user actually selected. Rescue/imputation
# state is carried separately in `estimate_type`, set by the inner classes.
.stamp_facade_identity <- function(res, facade_name) {
  res$modelName <- facade_name
  if ("facade" %in% colnames(res)) {
    res$facade <- NULL
  }
  dplyr::relocate(res, dplyr::any_of(c("modelName", "estimate_type")))
}

#' Limma contrast analysis facade
#'
#' Encapsulates the pipeline: \code{\link{strategy_limma}} ->
#' \code{\link{build_model_limma}} -> \code{\link{ContrastsLimma}}.
#'
#' Supports \code{options(prolfqua.vectorize = TRUE)} for faster
#' \code{\link{linfct_matrix_contrasts}} evaluation.
#' See \code{\link{build_contrast_analysis}} for details.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_protein_config()
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$rename_response("transformedIntensity")
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' fa <- ContrastsLimmaFacade$new(lfqdata, "~ group_", contrasts)
#' head(fa$get_contrasts())
#' fa$to_wide()
ContrastsLimmaFacade <- R6::R6Class(
  "ContrastsLimmaFacade",
  inherit = ContrastsFacadeBase,
  public = list(
    #' @field model ModelLimma object
    model = NULL,
    #' @field contrast ContrastsLimma object
    contrast = NULL,
    #' @field .lfqdata stored reference to input LFQData
    .lfqdata = NULL,
    #' @field .contrast_names names of the requested contrasts
    .contrast_names = NULL,
    #' @description
    #' initialize
    #' @param lfqdata LFQData object
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param weights column name for per-observation weights (default:
    #'   \code{lfqdata$nr_children_col()}). Pass \code{NULL} for unweighted.
    #' @param ... passed to \code{\link{strategy_limma}} (e.g. trend, robust)
    initialize = function(lfqdata, modelstr, contrasts, weights = lfqdata$nr_children_col(), ...) {
      .assert_aggregated_facade_input(lfqdata, "ContrastsLimmaFacade")
      self$.lfqdata <- lfqdata
      self$.contrast_names <- names(contrasts)
      response <- lfqdata$response()
      full_formula <- paste(response, modelstr)
      strat <- strategy_limma(full_formula, weights = weights, ...)
      self$model <- build_model_limma(lfqdata, strat)
      self$contrast <- ContrastsLimma$new(self$model, contrasts, model_name = "limma")
      self$facade_name <- "limma"
      self$.drop_na_diff <- TRUE
      self$config <- self$contrast$get_config()
    }
  )
)


#' Limma contrast analysis with LOD imputation facade
#'
#' Encapsulates the pipeline: \code{\link{strategy_limma}} ->
#' \code{\link{build_model_limma_impute}} -> \code{\link{ContrastsLimma}}.
#'
#' Proteins whose limma fit produces NA coefficients (typically from entire
#' missing groups) are recovered by imputing missing values with the limit of
#' detection (LOD) and refitting. The variance is borrowed from successful fits
#' and degrees of freedom are corrected so that inference is not artificially
#' precise from the constant imputation.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_protein_config(Nprot = 30, weight_missing = 0.5)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$rename_response("transformedIntensity")
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' fa <- ContrastsLimmaImputeFacade$new(lfqdata, "~ group_", contrasts)
#' head(fa$get_contrasts())
#' fa$to_wide()
ContrastsLimmaImputeFacade <- R6::R6Class(
  "ContrastsLimmaImputeFacade",
  inherit = ContrastsFacadeBase,
  public = list(
    #' @field model ModelLimma object (with imputed proteins)
    model = NULL,
    #' @field contrast ContrastsLimma object
    contrast = NULL,
    #' @field .lfqdata stored reference to input LFQData
    .lfqdata = NULL,
    #' @field .contrast_names names of the requested contrasts
    .contrast_names = NULL,
    #' @description
    #' initialize
    #' @param lfqdata LFQData object (aggregated to protein level)
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param lod numeric limit of detection; if NULL, auto-computed from data
    #' @param df_method "observed" uses max(n_observed - p, 1);
    #'   "borrowed" uses median df from successful fits
    #' @param weights column name for per-observation weights (default:
    #'   \code{lfqdata$nr_children_col()}). Pass \code{NULL} for unweighted.
    #' @param ... passed to \code{\link{strategy_limma}} (e.g. trend, robust)
    initialize = function(
      lfqdata,
      modelstr,
      contrasts,
      lod = NULL,
      df_method = c("observed", "borrowed"),
      weights = lfqdata$nr_children_col(),
      ...
    ) {
      .assert_aggregated_facade_input(lfqdata, "ContrastsLimmaImputeFacade")
      self$.lfqdata <- lfqdata
      self$.contrast_names <- names(contrasts)
      response <- lfqdata$response()
      full_formula <- paste(response, modelstr)
      strat <- strategy_limma(full_formula, weights = weights, ...)
      self$model <- build_model_limma_impute(lfqdata, strat, lod = lod, df_method = df_method)
      self$contrast <- ContrastsLimma$new(self$model, contrasts, model_name = "limma_impute")
      self$facade_name <- "limma_impute"
      self$.drop_na_diff <- TRUE
      self$config <- self$contrast$get_config()
    }
  )
)


#' Limma-voom contrast analysis facade
#'
#' Encapsulates the pipeline: \code{\link{strategy_limma}} ->
#' \code{\link{build_model_limma_voom}} -> \code{\link{ContrastsLimma}}.
#'
#' Uses vooma-style precision weights derived from a mean-variance trend,
#' optionally combined with external weights (e.g. peptide/precursor counts).
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_protein_config()
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$rename_response("transformedIntensity")
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' fa <- ContrastsLimmaVoomFacade$new(lfqdata, "~ group_", contrasts)
#' head(fa$get_contrasts())
#' fa$to_wide()
ContrastsLimmaVoomFacade <- R6::R6Class(
  "ContrastsLimmaVoomFacade",
  inherit = ContrastsFacadeBase,
  public = list(
    #' @field model ModelLimma object
    model = NULL,
    #' @field contrast ContrastsLimma object
    contrast = NULL,
    #' @field .lfqdata stored reference to input LFQData
    .lfqdata = NULL,
    #' @field .contrast_names names of the requested contrasts
    .contrast_names = NULL,
    #' @description
    #' initialize
    #' @param lfqdata LFQData object
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param weights column name for per-observation weights (default:
    #'   \code{lfqdata$nr_children_col()}). Pass \code{NULL} for unweighted.
    #' @param span lowess smoother span for vooma trend (default 0.5)
    #' @param plot logical; if TRUE, plot the mean-variance trend
    #' @param ... passed to \code{\link{strategy_limma}} (e.g. trend, robust)
    initialize = function(
      lfqdata,
      modelstr,
      contrasts,
      weights = lfqdata$nr_children_col(),
      span = 0.5,
      plot = FALSE,
      ...
    ) {
      .assert_aggregated_facade_input(lfqdata, "ContrastsLimmaVoomFacade")
      self$.lfqdata <- lfqdata
      self$.contrast_names <- names(contrasts)
      response <- lfqdata$response()
      full_formula <- paste(response, modelstr)
      strat <- strategy_limma(full_formula, weights = weights, ...)
      self$model <- build_model_limma_voom(lfqdata, strat, span = span, plot = plot)
      self$contrast <- ContrastsLimma$new(self$model, contrasts, model_name = "limma_voom")
      self$facade_name <- "limma_voom"
      self$.drop_na_diff <- TRUE
      self$config <- self$contrast$get_config()
    }
  )
)


#' Limma-voom contrast analysis with LOD imputation facade
#'
#' Encapsulates the pipeline: \code{\link{strategy_limma}} ->
#' \code{\link{build_model_limma_voom_impute}} -> \code{\link{ContrastsLimma}}.
#'
#' Combines vooma precision weights with LOD imputation for proteins whose
#' fit produces NA coefficients (typically from entire missing groups).
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_protein_config(Nprot = 30, weight_missing = 0.5)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$rename_response("transformedIntensity")
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' fa <- ContrastsLimmaVoomImputeFacade$new(lfqdata, "~ group_", contrasts)
#' head(fa$get_contrasts())
#' fa$to_wide()
ContrastsLimmaVoomImputeFacade <- R6::R6Class(
  "ContrastsLimmaVoomImputeFacade",
  inherit = ContrastsFacadeBase,
  public = list(
    #' @field model ModelLimma object (with imputed proteins)
    model = NULL,
    #' @field contrast ContrastsLimma object
    contrast = NULL,
    #' @field .lfqdata stored reference to input LFQData
    .lfqdata = NULL,
    #' @field .contrast_names names of the requested contrasts
    .contrast_names = NULL,
    #' @description
    #' initialize
    #' @param lfqdata LFQData object (aggregated to protein level)
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param lod numeric limit of detection; if NULL, auto-computed from data
    #' @param df_method "observed" uses max(n_observed - p, 1);
    #'   "borrowed" uses median df from successful fits
    #' @param weights column name for per-observation weights (default:
    #'   \code{lfqdata$nr_children_col()}). Pass \code{NULL} for unweighted.
    #' @param span lowess smoother span for vooma trend (default 0.5)
    #' @param plot logical; if TRUE, plot the mean-variance trend
    #' @param ... passed to \code{\link{strategy_limma}} (e.g. trend, robust)
    initialize = function(
      lfqdata,
      modelstr,
      contrasts,
      lod = NULL,
      df_method = c("observed", "borrowed"),
      weights = lfqdata$nr_children_col(),
      span = 0.5,
      plot = FALSE,
      ...
    ) {
      .assert_aggregated_facade_input(lfqdata, "ContrastsLimmaVoomImputeFacade")
      self$.lfqdata <- lfqdata
      self$.contrast_names <- names(contrasts)
      response <- lfqdata$response()
      full_formula <- paste(response, modelstr)
      strat <- strategy_limma(full_formula, weights = weights, ...)
      self$model <- build_model_limma_voom_impute(
        lfqdata,
        strat,
        lod = lod,
        df_method = df_method,
        span = span,
        plot = plot
      )
      self$contrast <- ContrastsLimma$new(self$model, contrasts, model_name = "limma_voom_impute")
      self$facade_name <- "limma_voom_impute"
      self$.drop_na_diff <- TRUE
      self$config <- self$contrast$get_config()
    }
  )
)


#' Limpa contrast analysis facade
#'
#' Encapsulates the pipeline: \code{\link{strategy_limpa}} ->
#' \code{\link{build_model_limpa}} -> \code{\link{ContrastsLimma}}.
#'
#' Operates as a \code{"same"}-level facade: it consumes whatever aggregated
#' LFQData prolfqua's normal aggregation pipeline produced and fits limpa at
#' that hierarchy level. If \code{config$opt_se} is set (e.g. when the input
#' came from \code{\link{AggregateLimpa}}), the per-observation standard error
#' is used as a vooma precision-weight predictor; otherwise plain vooma is
#' fit. The \code{config$nr_children} column is required and is used to flag
#' imputed observations (\code{nr_children == 0}) for vooma's
#' imputation-aware DF correction.
#'
#' For nested (peptide/precursor) input that should be rolled up to proteins
#' via limpa's DPC quantification, use \code{\link{ContrastsLimpaNestedFacade}}
#' instead — that facade owns the \code{\link{AggregateLimpa}} pre-step.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' istar <- prolfqua::sim_lfq_data_peptide_config()
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#' agg <- AggregateLimpa$new(lfqdata, "protein")
#' lfq_agg <- agg$aggregate()
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' fa <- ContrastsLimpaFacade$new(lfq_agg, "~ group_", contrasts)
#' head(fa$get_contrasts())
ContrastsLimpaFacade <- R6::R6Class(
  "ContrastsLimpaFacade",
  inherit = ContrastsFacadeBase,
  public = list(
    #' @field model ModelLimma object (from build_model_limpa)
    model = NULL,
    #' @field contrast ContrastsLimma object
    contrast = NULL,
    #' @field .lfqdata stored reference to input LFQData
    .lfqdata = NULL,
    #' @field .contrast_names names of the requested contrasts
    .contrast_names = NULL,
    #' @description
    #' initialize
    #' @param lfqdata aggregated LFQData. If \code{config$opt_se} is set
    #'   (e.g. from \code{\link{AggregateLimpa}}), the SE column is used as a
    #'   vooma precision-weight predictor; otherwise vooma is fit without an
    #'   external predictor.
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param plot logical; if TRUE, plot the vooma mean-variance trend
    #' @param span lowess smoother span (NULL = auto)
    #' @param ... passed to \code{\link{strategy_limpa}} (e.g. trend, robust)
    initialize = function(lfqdata, modelstr, contrasts, plot = FALSE, span = NULL, ...) {
      .assert_aggregated_facade_input(lfqdata, "ContrastsLimpaFacade")
      self$.lfqdata <- lfqdata
      self$.contrast_names <- names(contrasts)
      response <- lfqdata$response()
      full_formula <- paste(response, modelstr)
      strat <- strategy_limpa(full_formula, plot = plot, span = span, ...)
      self$model <- build_model_limpa(lfqdata, strat)
      self$contrast <- ContrastsLimma$new(self$model, contrasts, model_name = "limpa")
      self$facade_name <- "limpa"
      self$.drop_na_diff <- TRUE
      self$config <- self$contrast$get_config()
    }
  )
)


#' LM contrast analysis facade
#'
#' Encapsulates the pipeline: \code{\link{strategy_lm}} ->
#' \code{\link{build_model}} -> \code{\link{Contrasts}} ->
#' \code{\link{ContrastsModerated}}.
#'
#' Supports \code{options(prolfqua.vectorize = TRUE)} for faster contrast
#' computation. See \code{\link{build_contrast_analysis}} for details.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_protein_config()
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$rename_response("transformedIntensity")
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' fa <- ContrastsLMFacade$new(lfqdata, "~ group_", contrasts)
#' head(fa$get_contrasts())
#' fa$to_wide()
ContrastsLMFacade <- R6::R6Class(
  "ContrastsLMFacade",
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
    #' @param lfqdata LFQData object
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param weights column name for per-observation weights (default:
    #'   \code{lfqdata$nr_children_col()}). Pass \code{NULL} for unweighted.
    #' @param ... passed to \code{\link{strategy_lm}}
    initialize = function(lfqdata, modelstr, contrasts, weights = lfqdata$nr_children_col(), ...) {
      .assert_aggregated_facade_input(lfqdata, "ContrastsLMFacade")
      self$.lfqdata <- lfqdata
      self$.contrast_names <- names(contrasts)
      response <- lfqdata$response()
      full_formula <- paste(response, modelstr)
      strat <- strategy_lm(full_formula, weights = weights, ...)
      self$model <- build_model(lfqdata, strat)
      self$contrast <- ContrastsModerated$new(Contrasts$new(self$model, contrasts, model_name = "lm"))
      self$facade_name <- "lm"
      self$config <- self$contrast$get_config()
    }
  )
)


#' RLM contrast analysis facade
#'
#' Encapsulates the pipeline: \code{\link{strategy_rlm}} ->
#' \code{\link{build_model}} -> \code{\link{Contrasts}} ->
#' \code{\link{ContrastsModerated}}.
#'
#' Supports \code{options(prolfqua.vectorize = TRUE)} for faster contrast
#' computation. See \code{\link{build_contrast_analysis}} for details.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_protein_config()
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$rename_response("transformedIntensity")
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' fa <- ContrastsRLMFacade$new(lfqdata, "~ group_", contrasts)
#' head(fa$get_contrasts())
#' fa$to_wide()
ContrastsRLMFacade <- R6::R6Class(
  "ContrastsRLMFacade",
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
    #' @param lfqdata LFQData object
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param ... passed to \code{\link{strategy_rlm}}
    initialize = function(lfqdata, modelstr, contrasts, ...) {
      .assert_aggregated_facade_input(lfqdata, "ContrastsRLMFacade")
      self$.lfqdata <- lfqdata
      self$.contrast_names <- names(contrasts)
      response <- lfqdata$response()
      full_formula <- paste(response, modelstr)
      strat <- strategy_rlm(full_formula, ...)
      self$model <- build_model(lfqdata, strat)
      self$contrast <- ContrastsModerated$new(Contrasts$new(self$model, contrasts, model_name = "rlm"))
      self$facade_name <- "rlm"
      self$config <- self$contrast$get_config()
    }
  )
)


#' Rfit rank-based regression contrast analysis facade
#'
#' Encapsulates the pipeline: \code{\link{strategy_rfit}} ->
#' \code{\link{build_model}} -> \code{\link{Contrasts}} ->
#' \code{\link{ContrastsModerated}}, using \code{\link[Rfit]{rfit}} rank-based
#' linear regression. Output is the standard Wald contrast schema; the
#' empirical-Bayes moderation shrinks the rank-based scale across proteins.
#'
#' Unlike \code{\link{ContrastsLMFacade}} this backend takes no observation
#' weights (\code{rfit} has no \code{weights} argument), so \code{nr_children}
#' weighting is not applied.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_protein_config()
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$rename_response("transformedIntensity")
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' fa <- ContrastsRfitFacade$new(lfqdata, "~ group_", contrasts)
#' head(fa$get_contrasts())
#' fa$to_wide()
ContrastsRfitFacade <- R6::R6Class(
  "ContrastsRfitFacade",
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
    #' @param lfqdata LFQData object
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param ... passed to \code{\link{strategy_rfit}}
    initialize = function(lfqdata, modelstr, contrasts, ...) {
      .assert_aggregated_facade_input(lfqdata, "ContrastsRfitFacade")
      self$.lfqdata <- lfqdata
      self$.contrast_names <- names(contrasts)
      response <- lfqdata$response()
      full_formula <- paste(response, modelstr)
      strat <- strategy_rfit(full_formula, ...)
      self$model <- build_model(lfqdata, strat)
      self$contrast <- ContrastsModerated$new(Contrasts$new(self$model, contrasts, model_name = "rfit"))
      self$facade_name <- "rfit"
      self$config <- self$contrast$get_config()
    }
  )
)


#'
#' LM + missing-value imputation contrast analysis facade
#'
#' Encapsulates the pipeline: \code{\link{strategy_lm}} ->
#' \code{\link{build_model}} -> \code{\link{Contrasts}} ->
#' merge with \code{\link{ContrastsMissing}} ->
#' \code{\link{ContrastsModerated}}.
#'
#' Proteins without a fitted model get their contrasts filled in from the
#' group-mean imputation method (\code{\link{ContrastsMissing}}).
#'
#' @section Deprecated:
#' This facade is deprecated because its missing-data leg uses
#' \code{\link{ContrastsMissing}} — a pre-fitting group-mean
#' substitution, not a refit. Prefer the \code{lm_impute} facade
#' (\code{\link{ContrastsLMImputeFacade}}), which fits an LM for every
#' protein and refits failed/singular fits with LOD imputation and
#' borrowed per-protein variance, flagging rescued rows as
#' \code{estimate_type = "lod_imputed"} in the output.
#' Construction emits a \code{.Deprecated} warning; the entry is kept
#' in \code{\link{FACADE_REGISTRY}} so historical YAMLs continue to
#' work.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @seealso \code{\link{ContrastsLMImputeFacade}},
#'   \code{\link{build_model_impute}}, \code{\link{ContrastsMissing}}
#' @examples
#' # ContrastsMissing requires protein-level data (hierarchyDepth == len(hierarchy_keys()))
#' istar <- sim_lfq_data_protein_config(Nprot = 30)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$rename_response("transformedIntensity")
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' suppressWarnings(
#'   fa <- ContrastsLMMissingFacade$new(lfqdata, "~ group_", contrasts)
#' )
#' head(fa$get_contrasts())
#' fa$to_wide()
ContrastsLMMissingFacade <- R6::R6Class(
  "ContrastsLMMissingFacade",
  inherit = ContrastsFacadeBase,
  public = list(
    #' @field model Model object
    model = NULL,
    #' @field contrast ContrastsModerated object (merged with ContrastsMissing)
    contrast = NULL,
    #' @field missing_contrast ContrastsMissing object
    missing_contrast = NULL,
    #' @field merged merged contrast result list from merge_contrasts_results
    merged = NULL,
    #' @field .lfqdata stored reference to input LFQData
    .lfqdata = NULL,
    #' @field .contrast_names names of the requested contrasts
    .contrast_names = NULL,
    #' @description
    #' initialize
    #' @param lfqdata LFQData object
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param weights column name for per-observation weights (default:
    #'   \code{lfqdata$nr_children_col()}). Pass \code{NULL} for unweighted.
    #' @param ... passed to \code{\link{strategy_lm}}
    initialize = function(lfqdata, modelstr, contrasts, weights = lfqdata$nr_children_col(), ...) {
      .Deprecated(
        new = "ContrastsLMImputeFacade",
        msg = paste(
          "ContrastsLMMissingFacade (method = 'lm_missing') is",
          "deprecated: its second leg uses ContrastsMissing (group-mean",
          "substitution, no model fit). Prefer 'lm_impute' which refits",
          "failed/singular proteins with LOD imputation and borrowed",
          "variance, flagging rescued rows as estimate_type 'lod_imputed'.",
          "See ?ContrastsLMMissingFacade for migration."
        )
      )
      .assert_aggregated_facade_input(lfqdata, "ContrastsLMMissingFacade")
      self$.lfqdata <- lfqdata
      self$.contrast_names <- names(contrasts)
      response <- lfqdata$response()
      full_formula <- paste(response, modelstr)
      strat <- strategy_lm(full_formula, weights = weights, ...)
      self$model <- build_model(lfqdata, strat)
      base_contrast <- ContrastsModerated$new(Contrasts$new(self$model, contrasts))
      self$missing_contrast <- suppressWarnings(
        ContrastsMissing$new(lfqdata, contrasts = contrasts)
      )
      self$merged <- merge_contrasts_results(base_contrast, self$missing_contrast)
      self$contrast <- self$merged$merged
      self$facade_name <- "lm_missing"
      self$config <- self$contrast$get_config()
    },
    #' @description get \code{\link{ContrastsPlotter}} built from the stamped
    #'   facade output (so modelName is the facade key, not the merged leg names)
    #' @param fc_threshold fold change threshold
    #' @param fdr_threshold FDR threshold
    get_Plotter = function(fc_threshold = 1, fdr_threshold = 0.1) {
      ContrastsPlotter$new(
        self$get_contrasts(),
        subject_id = self$contrast$subject_id,
        fcthresh = fc_threshold,
        volcano = list(
          list(score = "p.value", thresh = fdr_threshold),
          list(score = "FDR", thresh = fdr_threshold)
        ),
        histogram = list(
          list(score = "p.value", xlim = c(0, 1, 0.05)),
          list(score = "FDR", xlim = c(0, 1, 0.05))
        ),
        score = list(list(score = "statistic", thresh = 5)),
        diff = "diff",
        contrast = "contrast"
      )
    },
    #' @description convert results to wide format from the stamped facade output
    #' @param columns value columns to pivot
    to_wide = function(columns = c("p.value", "FDR", "statistic")) {
      pivot_model_contrasts_to_wide(
        self$get_contrasts(),
        subject_id = self$contrast$subject_id,
        columns = c("diff", columns),
        contrast = "contrast"
      )
    }
  )
)


#' LM contrast analysis with LOD imputation facade
#'
#' Encapsulates the pipeline: \code{\link{strategy_lm}} ->
#' \code{\link{build_model}} (with \code{impute = TRUE}) ->
#' \code{\link{Contrasts}} -> \code{\link{ContrastsModerated}}.
#'
#' Proteins whose initial lm fit fails or produces NA coefficients are
#' re-fitted after imputing missing values with the limit of detection (LOD).
#' The covariance matrix is borrowed from successful fits so that the
#' variance is not underestimated by the constant imputation.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_protein_config(Nprot = 30, weight_missing = 0.5)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$rename_response("transformedIntensity")
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' fa <- ContrastsLMImputeFacade$new(lfqdata, "~ group_", contrasts)
#' head(fa$get_contrasts())
#' fa$to_wide()
ContrastsLMImputeFacade <- R6::R6Class(
  "ContrastsLMImputeFacade",
  inherit = ContrastsFacadeBase,
  public = list(
    #' @field model Model object (with imputed proteins)
    model = NULL,
    #' @field contrast ContrastsModerated object
    contrast = NULL,
    #' @field .lfqdata stored reference to input LFQData
    .lfqdata = NULL,
    #' @field .contrast_names names of the requested contrasts
    .contrast_names = NULL,
    #' @description
    #' initialize
    #' @param lfqdata LFQData object (aggregated to protein level)
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param lod numeric limit of detection; if NULL, auto-computed from data
    #' @param borrow_method "sigma" borrows scalar sigma and uses per-protein
    #'   (X'X)^-1; "vcov" borrows element-wise median of full vcov matrices
    #' @param df_method "observed" uses max(n_observed - p, 1);
    #'   "borrowed" uses median df from successful fits
    #' @param weights column name for per-observation weights (default:
    #'   \code{lfqdata$nr_children_col()}). Pass \code{NULL} for unweighted.
    #' @param ... passed to \code{\link{strategy_lm}}
    initialize = function(
      lfqdata,
      modelstr,
      contrasts,
      lod = NULL,
      borrow_method = c("sigma", "vcov"),
      df_method = c("observed", "borrowed"),
      weights = lfqdata$nr_children_col(),
      ...
    ) {
      .assert_aggregated_facade_input(lfqdata, "ContrastsLMImputeFacade")
      self$.lfqdata <- lfqdata
      self$.contrast_names <- names(contrasts)
      response <- lfqdata$response()
      full_formula <- paste(response, modelstr)
      strat <- strategy_lm(full_formula, weights = weights, ...)
      self$model <- build_model_impute(
        lfqdata,
        strat,
        lod = lod,
        borrow_method = borrow_method,
        df_method = df_method
      )
      self$contrast <- ContrastsModerated$new(Contrasts$new(self$model, contrasts, model_name = "lm_impute"))
      self$facade_name <- "lm_impute"
      self$config <- self$contrast$get_config()
    }
  )
)


#' Rfit rank-regression contrast analysis with LOD imputation facade
#'
#' The rank-based analogue of \code{\link{ContrastsLMImputeFacade}}.
#' Encapsulates the pipeline: \code{\link{strategy_rfit}} ->
#' \code{\link{build_model_impute}} -> \code{\link{Contrasts}} ->
#' \code{\link{ContrastsModerated}}, using \code{\link[Rfit]{rfit}}.
#'
#' Proteins whose initial \code{rfit} fit fails, is singular, or has an
#' incomplete design are re-fitted after imputing missing values with the
#' limit of detection (LOD). Uncertainty is borrowed from successful donor
#' fits so it is not underestimated by the constant imputation.
#'
#' Differences from \code{\link{ContrastsLMImputeFacade}}:
#' \itemize{
#'   \item Covariance borrowing is \strong{vcov-only} (element-wise median of
#'     the donors' full named covariance matrices). \code{rfit} exposes no
#'     \code{lm}-style \code{cov.unscaled}, so the scalar-sigma borrowing mode
#'     is not available and is not exposed.
#'   \item Borrowing is \strong{fail-hard} (\code{on_misalign = "fail"}): if the
#'     donor covariance matrices cannot be aligned by coefficient name, the
#'     rescue is skipped and those proteins remain in \code{get_missing()}
#'     rather than carrying an \code{lm}-specific approximation.
#'   \item No observation weights (\code{rfit} has no \code{weights} argument),
#'     matching plain \code{rfit}.
#' }
#'
#' \code{modelName} is the facade key \code{"rfit_impute"}; rescued rows are
#' flagged in the \code{estimate_type} column as \code{"lod_imputed"}, the same
#' convention used by \code{\link{ContrastsLMImputeFacade}}.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @seealso \code{\link{ContrastsLMImputeFacade}}, \code{\link{ContrastsRfitFacade}},
#'   \code{\link{build_model_impute}}
#' @examples
#' istar <- sim_lfq_data_protein_config(Nprot = 30, weight_missing = 0.5)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$rename_response("transformedIntensity")
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' fa <- ContrastsRfitImputeFacade$new(lfqdata, "~ group_", contrasts)
#' head(fa$get_contrasts())
#' fa$to_wide()
ContrastsRfitImputeFacade <- R6::R6Class(
  "ContrastsRfitImputeFacade",
  inherit = ContrastsFacadeBase,
  public = list(
    #' @field model Model object (with imputed proteins)
    model = NULL,
    #' @field contrast ContrastsModerated object
    contrast = NULL,
    #' @field .lfqdata stored reference to input LFQData
    .lfqdata = NULL,
    #' @field .contrast_names names of the requested contrasts
    .contrast_names = NULL,
    #' @description
    #' initialize
    #' @param lfqdata LFQData object (aggregated to protein level)
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param lod numeric limit of detection; if NULL, auto-computed from data
    #' @param df_method "observed" uses max(n_observed - p, 1);
    #'   "borrowed" uses median df from successful fits
    #' @param ... passed to \code{\link{strategy_rfit}}
    initialize = function(
      lfqdata,
      modelstr,
      contrasts,
      lod = NULL,
      df_method = c("observed", "borrowed"),
      ...
    ) {
      .assert_aggregated_facade_input(lfqdata, "ContrastsRfitImputeFacade")
      df_method <- match.arg(df_method)
      self$.lfqdata <- lfqdata
      self$.contrast_names <- names(contrasts)
      response <- lfqdata$response()
      full_formula <- paste(response, modelstr)
      strat <- strategy_rfit(full_formula, ...)
      self$model <- build_model_impute(
        lfqdata,
        strat,
        lod = lod,
        borrow_method = "vcov",
        df_method = df_method,
        on_misalign = "fail"
      )
      self$contrast <- ContrastsModerated$new(Contrasts$new(self$model, contrasts, model_name = "rfit_impute"))
      self$facade_name <- "rfit_impute"
      self$config <- self$contrast$get_config()
    }
  )
)


#' Firth logistic missingness contrast analysis facade
#'
#' Encapsulates the pipeline: encode missingness ->
#' \code{\link{build_model_glm_protein}} -> \code{\link{ContrastsFirth}}.
#'
#' Requires aggregated (protein-level) LFQData. For nested peptide-level
#' input use \code{\link{ContrastsFirthNestedFacade}}.
#'
#' Supports \code{options(prolfqua.vectorize = TRUE)} for faster contrast
#' computation. See \code{\link{build_contrast_analysis}} for details.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_protein_config(Nprot = 20, weight_missing = 0.5)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' fa <- ContrastsFirthFacade$new(lfqdata, "~ group_", contrasts)
#' head(fa$get_contrasts())
#' fa$to_wide()
ContrastsFirthFacade <- R6::R6Class(
  "ContrastsFirthFacade",
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
    #' @param lfqdata LFQData object
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param ... ignored; accepted for dispatch compatibility with other facades
    initialize = function(lfqdata, modelstr, contrasts, ...) {
      .assert_aggregated_facade_input(lfqdata, "ContrastsFirthFacade")
      self$.lfqdata <- lfqdata
      self$.contrast_names <- names(contrasts)
      self$model <- build_model_glm_protein(lfqdata, modelstr)
      self$contrast <- ContrastsFirth$new(self$model, contrasts, model_name = "firth")
      self$facade_name <- "firth"
      self$config <- self$contrast$get_config()
    }
  )
)


#' DEqMS contrast analysis facade
#'
#' Encapsulates the pipeline: \code{\link{strategy_lm}} ->
#' \code{\link{build_model}} -> \code{\link{Contrasts}} ->
#' merge with \code{\link{ContrastsMissing}} ->
#' \code{\link{ContrastsModeratedDEqMS}}.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_protein_config(Nprot = 50)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$rename_response("transformedIntensity")
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' fa <- ContrastsDEqMSFacade$new(lfqdata, "~ group_", contrasts)
#' head(fa$get_contrasts())
#' fa$to_wide()
ContrastsDEqMSFacade <- R6::R6Class(
  "ContrastsDEqMSFacade",
  inherit = ContrastsFacadeBase,
  public = list(
    #' @field model Model object
    model = NULL,
    #' @field contrast ContrastsModeratedDEqMS object
    contrast = NULL,
    #' @field .lfqdata stored reference to input LFQData
    .lfqdata = NULL,
    #' @field .contrast_names names of the requested contrasts
    .contrast_names = NULL,
    #' @description
    #' initialize
    #' @param lfqdata LFQData object
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param weights column name for per-observation weights (default:
    #'   \code{lfqdata$nr_children_col()}). Pass \code{NULL} for unweighted.
    #' @param ... passed to \code{\link{strategy_lm}}
    initialize = function(lfqdata, modelstr, contrasts, weights = lfqdata$nr_children_col(), ...) {
      .assert_aggregated_facade_input(lfqdata, "ContrastsDEqMSFacade")
      self$.lfqdata <- lfqdata
      self$.contrast_names <- names(contrasts)
      response <- lfqdata$response()
      full_formula <- paste(response, modelstr)
      strat <- strategy_lm(full_formula, weights = weights, ...)
      self$model <- build_model(lfqdata, strat)
      base_contrast <- Contrasts$new(self$model, contrasts, model_name = "deqms")
      count_column <- lfqdata$nr_children_col()
      count_df <- lfqdata$data_long() |>
        dplyr::select(dplyr::all_of(c(base_contrast$subject_id, count_column)))
      self$contrast <- ContrastsModeratedDEqMS$new(base_contrast, count_df = count_df, count_column = count_column)
      self$facade_name <- "deqms"
      self$config <- self$contrast$get_config()
    }
  )
)


#' DEqMS contrast analysis with vooma weights facade
#'
#' Combines vooma precision weights (mean-variance trend) with DEqMS
#' count-dependent variance moderation. Vooma runs \strong{without} external
#' weights so it captures only the mean-variance relationship; the peptide
#' count enters solely through DEqMS moderation, avoiding double-counting.
#'
#' Pipeline: \code{\link{strategy_limma}} (weights = NULL) ->
#' \code{\link{build_model_limma_voom}} -> \code{\link{ContrastsLimma}}
#' (eBayes = FALSE) -> \code{\link{ContrastsModeratedDEqMS}}.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_protein_config(Nprot = 50)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$rename_response("transformedIntensity")
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' fa <- ContrastsDEqMSVoomFacade$new(lfqdata, "~ group_", contrasts)
#' head(fa$get_contrasts())
ContrastsDEqMSVoomFacade <- R6::R6Class(
  "ContrastsDEqMSVoomFacade",
  inherit = ContrastsFacadeBase,
  public = list(
    #' @field model ModelLimma object
    model = NULL,
    #' @field contrast ContrastsModeratedDEqMS object
    contrast = NULL,
    #' @field .lfqdata stored reference to input LFQData
    .lfqdata = NULL,
    #' @field .contrast_names names of the requested contrasts
    .contrast_names = NULL,
    #' @description
    #' initialize
    #' @param lfqdata LFQData object
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param span lowess smoother span for vooma trend (default 0.5)
    #' @param plot logical; if TRUE, plot the mean-variance trend
    #' @param ... passed to \code{\link{strategy_limma}} (e.g. trend, robust)
    initialize = function(lfqdata, modelstr, contrasts, span = 0.5, plot = FALSE, ...) {
      .assert_aggregated_facade_input(lfqdata, "ContrastsDEqMSVoomFacade")
      self$.lfqdata <- lfqdata
      self$.contrast_names <- names(contrasts)
      response <- lfqdata$response()
      full_formula <- paste(response, modelstr)
      # No external weights — vooma handles mean-variance only
      strat <- strategy_limma(full_formula, weights = NULL, ...)
      self$model <- build_model_limma_voom(lfqdata, strat, span = span, plot = plot)
      # eBayes = FALSE: DEqMS replaces eBayes moderation
      base_contrast <- ContrastsLimma$new(self$model, contrasts, eBayes = FALSE, model_name = "deqms_voom")
      count_column <- lfqdata$nr_children_col()
      count_df <- lfqdata$data_long() |>
        dplyr::select(dplyr::all_of(c(base_contrast$subject_id, count_column)))
      self$contrast <- ContrastsModeratedDEqMS$new(base_contrast, count_df = count_df, count_column = count_column)
      self$facade_name <- "deqms_voom"
      self$config <- self$contrast$get_config()
    }
  )
)


.facade_registry_env <- new.env(parent = emptyenv())

.seed_facade_registry <- function() {
  entries <- list(
    lm = list(class = "ContrastsLMFacade", needs = "same"),
    lm_impute = list(class = "ContrastsLMImputeFacade", needs = "same"),
    lm_missing = list(class = "ContrastsLMMissingFacade", needs = "same"),
    limma = list(class = "ContrastsLimmaFacade", needs = "same"),
    limma_impute = list(class = "ContrastsLimmaImputeFacade", needs = "same"),
    limma_voom = list(class = "ContrastsLimmaVoomFacade", needs = "same"),
    limma_voom_impute = list(class = "ContrastsLimmaVoomImputeFacade", needs = "same"),
    limpa = list(class = "ContrastsLimpaFacade", needs = "same"),
    limpa_nested = list(class = "ContrastsLimpaNestedFacade", needs = "nested"),
    rlm = list(class = "ContrastsRLMFacade", needs = "same"),
    rfit = list(class = "ContrastsRfitFacade", needs = "same"),
    rfit_impute = list(class = "ContrastsRfitImputeFacade", needs = "same"),
    deqms = list(class = "ContrastsDEqMSFacade", needs = "same"),
    deqms_voom = list(class = "ContrastsDEqMSVoomFacade", needs = "same"),
    binomial_nested = list(class = "ContrastsBinomialNestedFacade", needs = "nested"),
    firth = list(class = "ContrastsFirthFacade", needs = "same"),
    firth_nested = list(class = "ContrastsFirthNestedFacade", needs = "nested"),
    lmer_nested = list(class = "ContrastsLmerNestedFacade", needs = "nested"),
    ropeca_nested = list(class = "ContrastsROPECANestedFacade", needs = "nested")
  )
  for (nm in names(entries)) {
    entry <- entries[[nm]]
    entry$package <- "prolfqua"
    entry$needs_saint_annotation <- FALSE
    assign(nm, entry, envir = .facade_registry_env)
  }
}
.seed_facade_registry()

#' Register a contrast facade class
#'
#' Adds an entry to the prolfqua facade registry so that
#' \code{\link{lookup_facade}} (and downstream consumers such as
#' \code{prolfquapp::R6_DEAnalyse$build_facade()}) can resolve it by
#' short name. Intended for downstream packages that ship additional
#' modelling backends (e.g. \code{prolfquasaint}) and want to plug into
#' the same dispatch path as the built-in facades.
#'
#' @param name short method name, e.g. \code{"saint"}.
#' @param class character; R6 facade class name to instantiate.
#' @param needs one of \code{"same"} (facade emits contrasts at the same
#'   hierarchy level as its input \emph{e.g.} protein -> protein FC,
#'   peptide -> peptide FC) or \code{"nested"} (facade takes child-level
#'   input and emits parent-level contrasts, \emph{e.g.}
#'   peptide -> protein FC).
#' @param package package the facade class lives in. Defaults to
#'   \code{"prolfqua"}.
#' @param needs_saint_annotation \code{TRUE} if the backend requires
#'   annotation reading in SAINT mode (bait / control columns).
#'   Default \code{FALSE}.
#' @param ... additional fields stored on the registry entry, available
#'   to consumers via \code{\link{lookup_facade}}.
#' @return The registered entry (invisibly).
#' @export
#' @family modelling
#' @examples
#' lookup_facade("lm")$class
#' # downstream packages call this from .onLoad():
#' # prolfqua::register_facade("saint", class = "ContrastsSAINTFacade",
#' #                           needs = "same", package = "prolfquasaint",
#' #                           needs_saint_annotation = TRUE)
register_facade <- function(name, class, needs, package = "prolfqua", needs_saint_annotation = FALSE, ...) {
  stopifnot(
    is.character(name),
    length(name) == 1,
    nzchar(name),
    is.character(class),
    length(class) == 1,
    nzchar(class),
    is.character(needs),
    length(needs) == 1,
    is.character(package),
    length(package) == 1
  )
  entry <- list(
    class = class,
    needs = needs,
    package = package,
    needs_saint_annotation = isTRUE(needs_saint_annotation),
    ...
  )
  assign(name, entry, envir = .facade_registry_env)
  invisible(entry)
}

#' Unregister a contrast facade class
#'
#' @param name short method name registered via
#'   \code{\link{register_facade}}.
#' @return Invisibly \code{TRUE} if the entry was removed, \code{FALSE}
#'   if no entry existed.
#' @export
#' @family modelling
#' @examples
#' register_facade("example", class = "ExampleFacade", needs = "same")
#' unregister_facade("example")
unregister_facade <- function(name) {
  if (exists(name, envir = .facade_registry_env, inherits = FALSE)) {
    rm(list = name, envir = .facade_registry_env)
    invisible(TRUE)
  } else {
    invisible(FALSE)
  }
}

#' Look up a contrast facade by short name
#'
#' @param name short method name (e.g. \code{"lm"}, \code{"saint"}).
#' @return Registry entry list with at least \code{class},
#'   \code{needs}, \code{package}, \code{needs_saint_annotation}; or
#'   \code{NULL} if no entry exists.
#' @export
#' @family modelling
#' @examples
#' lookup_facade("lm")$class
#' is.null(lookup_facade("does-not-exist"))
lookup_facade <- function(name) {
  if (length(name) != 1 || is.na(name) || !nzchar(name)) {
    return(NULL)
  }
  if (!exists(name, envir = .facade_registry_env, inherits = FALSE)) {
    return(NULL)
  }
  get(name, envir = .facade_registry_env, inherits = FALSE)
}

#' List currently registered contrast facades
#'
#' @return A list mirroring \code{\link{FACADE_REGISTRY}} but reflecting
#'   any entries added via \code{\link{register_facade}} from downstream
#'   packages.
#' @export
#' @family modelling
#' @examples
#' head(names(list_facades()))
list_facades <- function() {
  names_ <- ls(envir = .facade_registry_env, sorted = TRUE)
  out <- lapply(names_, function(nm) get(nm, envir = .facade_registry_env, inherits = FALSE))
  names(out) <- names_
  out
}

#' Registry of available contrast facade classes
#'
#' Read-only snapshot of the built-in prolfqua facade registry, taken at
#' package load time. It is derived from the single source of truth
#' (\code{.seed_facade_registry()}); use \code{\link{register_facade}} to add
#' entries from downstream packages, \code{\link{lookup_facade}} to resolve a
#' single entry, and \code{\link{list_facades}()} to enumerate the live
#' registry (which reflects downstream registrations).
#'
#' Each entry has fields \code{class}, \code{needs}, \code{package},
#' and \code{needs_saint_annotation}.
#'
#' @return A named list of built-in facade registry entries.
#' @export
#' @examples
#' names(FACADE_REGISTRY)
#' FACADE_REGISTRY$limma$class
FACADE_REGISTRY <- structure(
  list_facades(),
  class = c("facade_registry", "list")
)
