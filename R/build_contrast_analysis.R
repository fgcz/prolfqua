# build_contrast_analysis -----

#' Build a contrast analysis using one of several statistical methods
#'
#' A builder function that dispatches to the appropriate facade class based on
#' the chosen method. Each facade encapsulates the full pipeline from strategy
#' construction through modelling to contrast computation.
#'
#' @section Vectorized mode:
#' Set \code{options(prolfqua.vectorize = TRUE)} before calling this function
#' to activate vectorized implementations of \code{\link{compute_contrast}} and
#' \code{\link{linfct_matrix_contrasts}}. This affects all methods that use the
#' Wald test path (lm, rlm, firth, lmer) and can give a significant speed-up
#' for large datasets. Results are numerically identical. Example:
#' \preformatted{
#' options(prolfqua.vectorize = TRUE)
#' fa <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "lm")
#' options(prolfqua.vectorize = FALSE)  # restore default
#' }
#'
#' @param lfqdata an \code{\link{LFQData}} object
#' @param modelstr model formula string without the response variable
#'   (e.g. \code{"~ group_"}). The response is taken automatically from
#'   \code{lfqdata$get_config()$get_response()}.
#' @param contrasts named character vector of contrasts
#'   (e.g. \code{c("A_vs_B" = "group_A - group_B")})
#' @param method a registered facade key. The built-in keys are \code{"lm"},
#'   \code{"lm_impute"}, \code{"lm_missing"}, \code{"limma"},
#'   \code{"limma_impute"}, \code{"limma_voom"}, \code{"limma_voom_impute"},
#'   \code{"limpa"}, \code{"limpa_nested"}, \code{"rlm"}, \code{"rfit"},
#'   \code{"rfit_impute"}, \code{"deqms"}, \code{"deqms_voom"}, \code{"firth"},
#'   \code{"firth_nested"}, \code{"binomial_nested"}, \code{"lmer_nested"},
#'   \code{"ropeca_nested"};
#'   downstream packages may add more via \code{\link{register_facade}}. The
#'   authoritative list is \code{names(\link{list_facades}())}. Defaults to
#'   \code{"lm"}.
#' @param ... additional arguments forwarded to the underlying strategy function
#'   (e.g. \code{trend}, \code{robust} for \code{strategy_limma})
#' @return one of \code{\link{ContrastsLimmaFacade}},
#'   \code{\link{ContrastsLMFacade}}, \code{\link{ContrastsRLMFacade}},
#'   \code{\link{ContrastsRfitFacade}}, \code{\link{ContrastsRfitImputeFacade}},
#'   \code{\link{ContrastsLmerNestedFacade}},
#'   \code{\link{ContrastsLMMissingFacade}}, \code{\link{ContrastsLMImputeFacade}},
#'   \code{\link{ContrastsDEqMSFacade}},
#'   \code{\link{ContrastsROPECANestedFacade}}, \code{\link{ContrastsFirthFacade}},
#'   \code{\link{ContrastsFirthNestedFacade}}, \code{\link{ContrastsLimpaFacade}},
#'   \code{\link{ContrastsBinomialNestedFacade}}, or
#'   \code{\link{ContrastsLimpaNestedFacade}}
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_protein_config(Nprot = 20)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$rename_response("transformedIntensity")
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#'
#' fa_lm    <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "lm")
#' head(fa_lm$get_contrasts())
#'
#' fa_limma <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "limma")
#' head(fa_limma$get_contrasts())
#'
#' fa_miss <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "lm_missing")
#' head(fa_miss$get_contrasts())
#'
#' fa_deqms <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "deqms")
#' head(fa_deqms$get_contrasts())
#'
#' istar_pep <- sim_lfq_data_peptide_config()
#' lfqdata_pep <- LFQData$new(istar_pep$data, istar_pep$config)
#' lfqdata_pep <- lfqdata_pep$get_Transformer()$log2()$lfq
#'
#' fa_lmer <- build_contrast_analysis(
#'   lfqdata_pep,
#'   "~ group_ + (1 | peptide_Id) + (1 | sampleName)",
#'   contrasts,
#'   method = "lmer_nested"
#' )
#' head(fa_lmer$get_contrasts())
#'
#' fa_ropeca <- build_contrast_analysis(lfqdata_pep, "~ group_", contrasts, method = "ropeca_nested")
#' head(fa_ropeca$get_contrasts())
#'
#' fa_firth <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "firth")
#' head(fa_firth$get_contrasts())
build_contrast_analysis <- function(
  lfqdata,
  modelstr,
  contrasts,
  method = "lm",
  ...
) {
  # The facade registry (seeded by .seed_facade_registry(), extensible via
  # register_facade()) is the single source of truth for method -> class.
  choices <- names(list_facades())
  method <- tryCatch(
    match.arg(method, choices),
    error = function(e) {
      abort_bad_argument(
        "method",
        must = paste0("be one of: ", paste(choices, collapse = ", ")),
        not = paste(method, collapse = ", "),
        parent = e
      )
    }
  )
  entry <- lookup_facade(method)
  # Decoys must not enter the fit or the shared variance pool (limma prior /
  # DEqMS variance-count trend), where they would perturb *target* q-values.
  # Opt-in: only when the config declares `pattern_decoys` (NULL leaves existing
  # behaviour untouched). Normalization is expected to have already run on the
  # full data upstream; here we drop decoys so the model sees targets only.
  cfg <- lfqdata$get_config()
  if (!is.null(cfg$pattern_decoys)) {
    top <- lfqdata$hierarchy_keys()[1]
    n_before <- length(unique(lfqdata$data_long()[[top]]))
    lfqdata <- lfqdata$remove_decoys()
    n_dropped <- n_before - length(unique(lfqdata$data_long()[[top]]))
    if (n_dropped > 0) {
      message(
        "build_contrast_analysis: dropped ",
        n_dropped,
        " decoy ",
        top,
        " before modelling (targets-only fit)."
      )
    }
  }
  if (!is.null(entry$builder)) {
    return(entry$builder(lfqdata, modelstr, contrasts, ...))
  }
  facade_class <- utils::getFromNamespace(entry$class, entry$package %||% "prolfqua")
  facade_class$new(lfqdata, modelstr, contrasts, ...)
}
