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
#'   \code{lfqdata$config$get_response()}.
#' @param contrasts named character vector of contrasts
#'   (e.g. \code{c("A_vs_B" = "group_A - group_B")})
#' @param method one of \code{"lm"}, \code{"rlm"}, \code{"lmer"}, \code{"lm_missing"},
#'   \code{"lm_impute"}, \code{"limma"}, \code{"limma_impute"}, \code{"limma_voom"},
#'   \code{"limma_voom_impute"}, \code{"deqms"}, \code{"deqms_voom"},
#'   \code{"ropeca"}, \code{"firth"}, \code{"limpa"}
#' @param ... additional arguments forwarded to the underlying strategy function
#'   (e.g. \code{trend}, \code{robust} for \code{strategy_limma})
#' @return one of \code{\link{ContrastsLimmaFacade}},
#'   \code{\link{ContrastsLMFacade}}, \code{\link{ContrastsRLMFacade}},
#'   \code{\link{ContrastsLmerFacade}},
#'   \code{\link{ContrastsLMMissingFacade}}, \code{\link{ContrastsLMImputeFacade}},
#'   \code{\link{ContrastsDEqMSFacade}},
#'   \code{\link{ContrastsROPECAFacade}}, or \code{\link{ContrastsFirthFacade}}
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
#' istar_pep$config <- old2new(istar_pep$config)
#' lfqdata_pep <- LFQData$new(istar_pep$data, istar_pep$config)
#' lfqdata_pep <- lfqdata_pep$get_Transformer()$log2()$lfq
#'
#' fa_lmer <- build_contrast_analysis(
#'   lfqdata_pep,
#'   "~ group_ + (1 | peptide_Id) + (1 | sampleName)",
#'   contrasts,
#'   method = "lmer"
#' )
#' head(fa_lmer$get_contrasts())
#'
#' fa_ropeca <- build_contrast_analysis(lfqdata_pep, "~ group_", contrasts, method = "ropeca")
#' head(fa_ropeca$get_contrasts())
#'
#' fa_firth <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "firth")
#' head(fa_firth$get_contrasts())
build_contrast_analysis <- function(
  lfqdata,
  modelstr,
  contrasts,
  method = c(
    "lm",
    "rlm",
    "lmer",
    "lm_missing",
    "lm_impute",
    "limma",
    "limma_impute",
    "limma_voom",
    "limma_voom_impute",
    "deqms",
    "deqms_voom",
    "ropeca",
    "firth",
    "limpa"
  ),
  ...
) {
  method <- match.arg(method)
  switch(
    method,
    lm = ContrastsLMFacade$new(lfqdata, modelstr, contrasts, ...),
    rlm = ContrastsRLMFacade$new(lfqdata, modelstr, contrasts, ...),
    lmer = ContrastsLmerFacade$new(lfqdata, modelstr, contrasts, ...),
    lm_missing = ContrastsLMMissingFacade$new(lfqdata, modelstr, contrasts, ...),
    lm_impute = ContrastsLMImputeFacade$new(lfqdata, modelstr, contrasts, ...),
    limma = ContrastsLimmaFacade$new(lfqdata, modelstr, contrasts, ...),
    limma_impute = ContrastsLimmaImputeFacade$new(lfqdata, modelstr, contrasts, ...),
    limma_voom = ContrastsLimmaVoomFacade$new(lfqdata, modelstr, contrasts, ...),
    limma_voom_impute = ContrastsLimmaVoomImputeFacade$new(lfqdata, modelstr, contrasts, ...),
    deqms = ContrastsDEqMSFacade$new(lfqdata, modelstr, contrasts, ...),
    deqms_voom = ContrastsDEqMSVoomFacade$new(lfqdata, modelstr, contrasts, ...),
    ropeca = ContrastsROPECAFacade$new(lfqdata, modelstr, contrasts, ...),
    firth = ContrastsFirthFacade$new(lfqdata, modelstr, contrasts),
    limpa = ContrastsLimpaFacade$new(lfqdata, modelstr, contrasts, ...)
  )
}
