# ContrastsFacades -----

.assert_aggregated_facade_input <- function(lfqdata, facade_name) {
  subject_id <- lfqdata$subject_Id()
  hierarchy_keys <- lfqdata$config$hierarchy_keys()
  if (!identical(subject_id, hierarchy_keys)) {
    stop(
      facade_name,
      " requires aggregated LFQData. ",
      "`lfqdata$subject_Id()` must equal `lfqdata$config$hierarchy_keys()`. ",
      "Aggregate first.",
      call. = FALSE
    )
  }
}

.assert_nested_facade_input <- function(lfqdata, facade_name) {
  subject_id <- lfqdata$subject_Id()
  hierarchy_keys <- lfqdata$config$hierarchy_keys()
  if (!(all(subject_id %in% hierarchy_keys) && length(subject_id) < length(hierarchy_keys))) {
    stop(
      facade_name,
      " requires LFQData with additional hierarchy below `subject_Id()`. ",
      "`lfqdata$subject_Id()` must be a strict subset of ",
      "`lfqdata$config$hierarchy_keys()`. Do not aggregate first.",
      call. = FALSE
    )
  }
}

.add_facade_column <- function(res, facade_name) {
  if (!("facade" %in% colnames(res))) {
    res <- dplyr::mutate(res, facade = facade_name, .before = 1)
  }
  res
}

#' Limma contrast analysis facade
#'
#' Encapsulates the pipeline: \code{\link{strategy_limma}} ->
#' \code{\link{build_model_limma}} -> \code{\link{ContrastsLimma}}.
#'
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
  public = list(
    #' @field model ModelLimma object
    model = NULL,
    #' @field contrast ContrastsLimma object
    contrast = NULL,
    #' @description
    #' initialize
    #' @param lfqdata LFQData object
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param ... passed to \code{\link{strategy_limma}} (e.g. trend, robust)
    initialize = function(lfqdata, modelstr, contrasts, ...) {
      .assert_aggregated_facade_input(lfqdata, "ContrastsLimmaFacade")
      response <- lfqdata$config$get_response()
      full_formula <- paste(response, modelstr)
      strat          <- strategy_limma(full_formula, ...)
      self$model     <- build_model_limma(lfqdata, strat)
      self$contrast  <- ContrastsLimma$new(self$model, contrasts)
    },
    #' @description get contrast results
    #' @param ... passed to ContrastsLimma$get_contrasts
    get_contrasts = function(...) {
      .add_facade_column(self$contrast$get_contrasts(...), "limma")
    },
    #' @description get ContrastsPlotter
    #' @param ... passed to ContrastsLimma$get_Plotter
    get_Plotter   = function(...) self$contrast$get_Plotter(...),
    #' @description convert results to wide format
    #' @param ... passed to ContrastsLimma$to_wide
    to_wide       = function(...) self$contrast$to_wide(...)
  )
)


#' LM contrast analysis facade
#'
#' Encapsulates the pipeline: \code{\link{strategy_lm}} ->
#' \code{\link{build_model}} -> \code{\link{Contrasts}} ->
#' \code{\link{ContrastsModerated}}.
#'
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
  public = list(
    #' @field model Model object
    model = NULL,
    #' @field contrast ContrastsModerated object
    contrast = NULL,
    #' @description
    #' initialize
    #' @param lfqdata LFQData object
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param ... passed to \code{\link{strategy_lm}}
    initialize = function(lfqdata, modelstr, contrasts, ...) {
      .assert_aggregated_facade_input(lfqdata, "ContrastsLMFacade")
      response <- lfqdata$config$get_response()
      full_formula <- paste(response, modelstr)
      strat          <- strategy_lm(full_formula, ...)
      self$model     <- build_model(lfqdata, strat)
      self$contrast  <- ContrastsModerated$new(Contrasts$new(self$model, contrasts))
    },
    #' @description get contrast results
    #' @param ... passed to ContrastsModerated$get_contrasts
    get_contrasts = function(...) {
      .add_facade_column(self$contrast$get_contrasts(...), "lm")
    },
    #' @description get ContrastsPlotter
    #' @param ... passed to ContrastsModerated$get_Plotter
    get_Plotter   = function(...) self$contrast$get_Plotter(...),
    #' @description convert results to wide format
    #' @param ... passed to ContrastsModerated$to_wide
    to_wide       = function(...) self$contrast$to_wide(...)
  )
)


#' Lmer contrast analysis facade
#'
#' Encapsulates the pipeline: \code{\link{strategy_lmer}} ->
#' \code{\link{build_model}} -> \code{\link{Contrasts}} ->
#' \code{\link{ContrastsModerated}}.
#'
#' This facade requires data with hierarchy below the analysis subject, for
#' example peptide-level measurements nested within proteins.
#'
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_peptide_config()
#' istar$config <- old2new(istar$config)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' fa <- ContrastsLmerFacade$new(
#'   lfqdata,
#'   "~ group_ + (1 | peptide_Id) + (1 | sampleName)",
#'   contrasts
#' )
#' head(fa$get_contrasts())
#' fa$to_wide()
ContrastsLmerFacade <- R6::R6Class(
  "ContrastsLmerFacade",
  public = list(
    #' @field model Model object
    model = NULL,
    #' @field contrast ContrastsModerated object
    contrast = NULL,
    #' @description
    #' initialize
    #' @param lfqdata LFQData object
    #' @param modelstr model formula string (e.g. "~ group_ + (1 | peptide_Id)")
    #' @param contrasts named character vector of contrasts
    #' @param ... passed to \code{\link{strategy_lmer}}
    initialize = function(lfqdata, modelstr, contrasts, ...) {
      .assert_nested_facade_input(lfqdata, "ContrastsLmerFacade")
      response <- lfqdata$config$get_response()
      full_formula <- paste(response, modelstr)
      strat          <- strategy_lmer(full_formula, ...)
      self$model     <- build_model(lfqdata, strat)
      self$contrast  <- ContrastsModerated$new(Contrasts$new(self$model, contrasts))
    },
    #' @description get contrast results
    #' @param ... passed to ContrastsModerated$get_contrasts
    get_contrasts = function(...) {
      .add_facade_column(self$contrast$get_contrasts(...), "lmer")
    },
    #' @description get ContrastsPlotter
    #' @param ... passed to ContrastsModerated$get_Plotter
    get_Plotter   = function(...) self$contrast$get_Plotter(...),
    #' @description convert results to wide format
    #' @param ... passed to ContrastsModerated$to_wide
    to_wide       = function(...) self$contrast$to_wide(...)
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
#' @export
#' @family modelling
#' @examples
#' # ContrastsMissing requires protein-level data (hierarchyDepth == len(hierarchy_keys()))
#' istar <- sim_lfq_data_protein_config(Nprot = 30)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$rename_response("transformedIntensity")
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' fa <- ContrastsLMMissingFacade$new(lfqdata, "~ group_", contrasts)
#' head(fa$get_contrasts())
#' fa$to_wide()
ContrastsLMMissingFacade <- R6::R6Class(
  "ContrastsLMMissingFacade",
  public = list(
    #' @field model Model object
    model = NULL,
    #' @field contrast ContrastsModerated object (merged with ContrastsMissing)
    contrast = NULL,
    #' @field missing_contrast ContrastsMissing object
    missing_contrast = NULL,
    #' @field merged merged contrast result list from merge_contrasts_results
    merged = NULL,
    #' @description
    #' initialize
    #' @param lfqdata LFQData object
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param ... passed to \code{\link{strategy_lm}}
    initialize = function(lfqdata, modelstr, contrasts, ...) {
      .assert_aggregated_facade_input(lfqdata, "ContrastsLMMissingFacade")
      response <- lfqdata$config$get_response()
      full_formula <- paste(response, modelstr)
      strat                  <- strategy_lm(full_formula, ...)
      self$model             <- build_model(lfqdata, strat)
      base_contrast          <- ContrastsModerated$new(Contrasts$new(self$model, contrasts))
      self$missing_contrast  <- ContrastsMissing$new(lfqdata, contrasts = contrasts)
      self$merged            <- merge_contrasts_results(base_contrast, self$missing_contrast)
      self$contrast          <- self$merged$merged
    },
    #' @description get contrast results
    #' @param ... passed to ContrastsTable$get_contrasts
    get_contrasts = function(...) {
      .add_facade_column(self$contrast$get_contrasts(...), "lm_missing")
    },
    #' @description get ContrastsPlotter
    #' @param ... passed to ContrastsTable$get_Plotter
    get_Plotter   = function(...) self$contrast$get_Plotter(...),
    #' @description convert results to wide format
    #' @param ... passed to ContrastsTable$to_wide
    to_wide       = function(...) self$contrast$to_wide(...)
  )
)


#' DEqMS contrast analysis facade
#'
#' Encapsulates the pipeline: \code{\link{strategy_lm}} ->
#' \code{\link{build_model}} -> \code{\link{Contrasts}} ->
#' merge with \code{\link{ContrastsMissing}} ->
#' \code{\link{ContrastsModeratedDEqMS}}.
#'
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
  public = list(
    #' @field model Model object
    model = NULL,
    #' @field contrast ContrastsModeratedDEqMS object
    contrast = NULL,
    #' @description
    #' initialize
    #' @param lfqdata LFQData object
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param ... passed to \code{\link{strategy_lm}}
    initialize = function(lfqdata, modelstr, contrasts, ...) {
      .assert_aggregated_facade_input(lfqdata, "ContrastsDEqMSFacade")
      response <- lfqdata$config$get_response()
      full_formula <- paste(response, modelstr)
      strat         <- strategy_lm(full_formula, ...)
      self$model    <- build_model(lfqdata, strat)
      base_contrast <- Contrasts$new(self$model, contrasts)
      count_column <- lfqdata$config$nr_children
      count_df <- lfqdata$data |>
        dplyr::select(dplyr::all_of(c(base_contrast$subject_Id, count_column)))
      self$contrast <- ContrastsModeratedDEqMS$new(base_contrast,
                                                    count_df = count_df,
                                                    count_column = count_column)
    },
    #' @description get contrast results
    #' @param ... passed to ContrastsModeratedDEqMS$get_contrasts
    get_contrasts = function(...) {
      .add_facade_column(self$contrast$get_contrasts(...), "deqms")
    },
    #' @description get ContrastsPlotter
    #' @param ... passed to ContrastsModeratedDEqMS$get_Plotter
    get_Plotter   = function(...) self$contrast$get_Plotter(...),
    #' @description convert results to wide format
    #' @param ... passed to ContrastsModeratedDEqMS$to_wide
    to_wide       = function(...) self$contrast$to_wide(...)
  )
)


#' ROPECA contrast analysis facade
#'
#' Encapsulates the pipeline: \code{\link{strategy_lm}} ->
#' \code{\link{build_model}} -> \code{\link{Contrasts}} ->
#' \code{\link{ContrastsROPECA}}.
#'
#' ROPECA operates on peptide-level data and aggregates peptide-level
#' p-values to the protein level. The \code{lfqdata} object must contain
#' peptide-level data (i.e. \code{hierarchyDepth >= 2}).
#'
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_peptide_config()
#' istar$config <- old2new(istar$config)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' fa <- ContrastsROPECAFacade$new(lfqdata, "~ group_", contrasts)
#' head(fa$get_contrasts())
#' fa$to_wide()
ContrastsROPECAFacade <- R6::R6Class(
  "ContrastsROPECAFacade",
  public = list(
    #' @field model Model object (peptide-level)
    model = NULL,
    #' @field contrast ContrastsROPECA object
    contrast = NULL,
    #' @description
    #' initialize
    #' @param lfqdata LFQData object with peptide-level data (hierarchyDepth >= 2)
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param ... passed to \code{\link{strategy_lm}}
    initialize = function(lfqdata, modelstr, contrasts, ...) {
      .assert_nested_facade_input(lfqdata, "ContrastsROPECAFacade")
      response <- lfqdata$config$get_response()
      full_formula <- paste(response, modelstr)
      strat         <- strategy_lm(full_formula, ...)
      # ROPECA requires subject_Id with >1 element (protein + peptide level).
      # Use all hierarchy keys, not just those at hierarchyDepth.
      subject_Id    <- lfqdata$config$hierarchy_keys()
      self$model    <- build_model(lfqdata, strat, subject_Id = subject_Id)
      self$contrast <- ContrastsROPECA$new(Contrasts$new(self$model, contrasts))
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
    #'   \item \code{conf.low/conf.high} via \code{diff ± qt(0.975, df) * |std.error|}
    #' }
    get_contrasts = function() {
      res <- self$contrast$get_contrasts(all = TRUE)
      res <- dplyr::rename(res,
        p.value = "beta.based.significance",
        FDR     = "FDR.beta.based.significance"
      )
      # Derive std.error from median fold change / median t-statistic
      res$std.error <- ifelse(res$statistic != 0,
                              res$diff / res$statistic, NA_real_)
      # sigma: MAD of peptide-level fold changes as robust dispersion proxy
      res$sigma <- res$mad.estimate
      # df: number of contributing peptides
      res$df <- res$n_not_na
      # Confidence intervals from derived std.error and df
      prqt <- stats::qt(0.975, df = pmax(res$df, 1))
      res$conf.low  <- res$diff - prqt * abs(res$std.error)
      res$conf.high <- res$diff + prqt * abs(res$std.error)
      # Select standard columns only
      protein_Id <- self$contrast$subject_Id[1]
      standard_cols <- c(protein_Id, "modelName", "contrast", "avgAbd",
                         "diff", "FDR", "statistic", "std.error", "df",
                         "p.value", "conf.low", "conf.high", "sigma")
      res <- res[, standard_cols, drop = FALSE]
      .add_facade_column(res, "ropeca")
    },
    #' @description get ContrastsPlotter (uses standardized column names)
    #' @param FCthreshold fold change threshold
    #' @param FDRthreshold FDR threshold
    get_Plotter = function(FCthreshold = 2, FDRthreshold = 0.1) {
      contrast_result <- self$get_contrasts()
      protein_Id <- self$contrast$subject_Id[1]
      ContrastsPlotter$new(
        contrast_result,
        subject_Id = protein_Id,
        fcthresh   = FCthreshold,
        volcano    = list(list(score = "FDR",     thresh = FDRthreshold)),
        histogram  = list(list(score = "p.value", xlim   = c(0, 1, 0.05)),
                          list(score = "FDR",     xlim   = c(0, 1, 0.05))),
        modelName  = "modelName",
        diff       = "diff",
        contrast   = "contrast"
      )
    },
    #' @description convert results to wide format
    #' @param columns value columns to pivot, default \code{c("p.value", "FDR", "statistic")}
    to_wide = function(columns = c("p.value", "FDR", "statistic")) {
      contrast_result <- self$get_contrasts()
      protein_Id <- self$contrast$subject_Id[1]
      pivot_model_contrasts_2_Wide(
        contrast_result,
        subject_Id = protein_Id,
        columns    = c("diff", columns),
        contrast   = "contrast"
      )
    }
  )
)
