# build_contrast_analysis -----

#' Build a contrast analysis using one of several statistical methods
#'
#' A builder function that dispatches to the appropriate facade class based on
#' the chosen method. Each facade encapsulates the full pipeline from strategy
#' construction through modelling to contrast computation.
#'
#' @param lfqdata an \code{\link{LFQData}} object
#' @param modelstr model formula string without the response variable
#'   (e.g. \code{"~ group_"}). The response is taken automatically from
#'   \code{lfqdata$config$get_response()}.
#' @param contrasts named character vector of contrasts
#'   (e.g. \code{c("A_vs_B" = "group_A - group_B")})
#' @param method one of \code{"lm"}, \code{"lm_missing"}, \code{"limma"},
#'   \code{"deqms"}, \code{"ropeca"}
#' @param count_df data.frame with subject_Id columns plus a count column;
#'   required when \code{method = "deqms"}
#' @param count_column name of the count column in \code{count_df};
#'   required when \code{method = "deqms"}
#' @param ... additional arguments forwarded to the underlying strategy function
#'   (e.g. \code{trend}, \code{robust} for \code{strategy_limma})
#' @return one of \code{\link{ContrastsLimmaFacade}},
#'   \code{\link{ContrastsLMFacade}}, \code{\link{ContrastsLMMissingFacade}},
#'   \code{\link{ContrastsDEqMSFacade}}, or \code{\link{ContrastsROPECAFacade}}
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_peptide_config()
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#'
#' fa_lm    <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "lm")
#' head(fa_lm$get_contrasts())
#'
#' fa_limma <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "limma")
#' head(fa_limma$get_contrasts())
#'
#' fa_ropeca <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "ropeca")
#' head(fa_ropeca$get_contrasts())
#'
#' # lm_missing requires protein-level data
#' istar_prot <- sim_lfq_data_protein_config(Nprot = 20)
#' lfqdata_prot <- LFQData$new(istar_prot$data, istar_prot$config)
#' lfqdata_prot$rename_response("transformedIntensity")
#' fa_miss <- build_contrast_analysis(lfqdata_prot, "~ group_", contrasts, method = "lm_missing")
#' head(fa_miss$get_contrasts())
build_contrast_analysis <- function(lfqdata,
                                    modelstr,
                                    contrasts,
                                    method = c("lm", "lm_missing", "limma", "deqms", "ropeca"),
                                    count_df     = NULL,
                                    count_column = NULL,
                                    ...) {
  method <- match.arg(method)
  switch(method,
    lm         = ContrastsLMFacade$new(lfqdata, modelstr, contrasts, ...),
    lm_missing = ContrastsLMMissingFacade$new(lfqdata, modelstr, contrasts, ...),
    limma      = ContrastsLimmaFacade$new(lfqdata, modelstr, contrasts, ...),
    deqms      = {
      if (is.null(count_df) || is.null(count_column)) {
        stop("count_df and count_column are required for method = 'deqms'")
      }
      ContrastsDEqMSFacade$new(lfqdata, modelstr, contrasts,
                                count_df = count_df,
                                count_column = count_column, ...)
    },
    ropeca     = ContrastsROPECAFacade$new(lfqdata, modelstr, contrasts, ...)
  )
}
