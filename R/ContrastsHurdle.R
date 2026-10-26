# ContrastsHurdle ----

#' Hurdle contrasts: intensity and detection combined with Fisher's method
#'
#' Joins the contrasts of a fitted intensity facade (e.g. \code{lm_impute}) and a
#' fitted detection facade (\code{firth_nested}, \code{binomial_nested} or
#' \code{binomial_firth_nested}) by
#' subject and contrast, and replaces \code{p.value} and \code{FDR} with Fisher's
#' combined test, as in msqrob2's \code{hypothesisTestHurdle}. No model is fitted.
#'
#' For each subject and contrast, \eqn{X^2 = -2 \sum \log p} over the available
#' component p-values, with two degrees of freedom per component, so a row with a
#' single component keeps that component's p-value. The detection p-value is left
#' out for subjects whose child features are detected in every sample, because
#' detection carries no information there. \code{FDR} is the Benjamini-Hochberg
#' adjustment within each contrast.
#'
#' \code{diff}, \code{avgAbd}, \code{std.error}, \code{statistic}, \code{df},
#' the confidence interval, \code{sigma} and \code{estimate_type} come from the
#' intensity component. They are \code{NA} for subjects only the detection
#' component estimates; with \code{lm} these are the proteins detected in one
#' group only, so \code{filter_significant()} and \code{get_ora()} drop them.
#' Pair with \code{lm_impute}, which estimates their \code{diff}.
#'
#' Added columns: \code{p.value.intensity}, \code{FDR.intensity},
#' \code{diff.detection} (log odds ratio), \code{p.value.detection},
#' \code{FDR.detection} (all as reported by the components),
#' \code{chisq.hurdle}, \code{df.hurdle}, and \code{hurdle_components}
#' (\code{both}, \code{intensity} or \code{detection}: the p-values combined).
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_peptide_config(Nprot = 20, weight_missing = 0.5, seed = 3)
#' lfq_peptide <- LFQData$new(istar$data, istar$config)
#' lfq_peptide <- lfq_peptide$get_Transformer()$log2()$lfq
#' lfq_protein <- lfq_peptide$get_Aggregator()$aggregate()
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' intensity <- build_contrast_analysis(lfq_protein, "~ group_", contrasts, method = "lm_impute")
#' detection <- build_contrast_analysis(lfq_peptide, "~ group_", contrasts, method = "binomial_nested")
#' hurdle <- ContrastsHurdleFacade$new(intensity, detection)
#' res <- hurdle$get_contrasts()
#' table(res$hurdle_components)
#' head(res)
#' hurdle$to_wide()
ContrastsHurdleFacade <- R6::R6Class(
  "ContrastsHurdleFacade",
  inherit = ContrastsInterface,
  public = list(
    #' @field intensity fitted intensity facade
    intensity = NULL,
    #' @field detection fitted detection facade (\code{firth_nested}, \code{binomial_nested} or
    #'   \code{binomial_firth_nested})
    detection = NULL,
    #' @field model_name name stamped into the \code{modelName} column
    model_name = NULL,
    #' @field subject_id columns identifying the subject (e.g. protein_Id)
    subject_id = NULL,
    #' @description
    #' Combine a fitted intensity facade and a fitted detection facade.
    #' @param intensity fitted intensity facade, e.g. from
    #'   \code{build_contrast_analysis(method = "lm_impute")}
    #' @param detection fitted \code{firth_nested}, \code{binomial_nested} or
    #'   \code{binomial_firth_nested} facade with the same subject id and contrasts
    #' @param model_name name stamped into the \code{modelName} column
    initialize = function(intensity, detection, model_name = "hurdle") {
      stopifnot(inherits(intensity, "ContrastsInterface"), inherits(detection, "ContrastsInterface"))
      subject_id <- intensity$get_config()$subject_id
      if (!identical(subject_id, detection$get_config()$subject_id)) {
        stop(
          "intensity and detection must have the same subject_id: ",
          paste(subject_id, collapse = ", "),
          " vs ",
          paste(detection$get_config()$subject_id, collapse = ", ")
        )
      }
      intensity_contrasts <- dplyr::ungroup(intensity$get_contrasts())
      detection_contrasts <- dplyr::ungroup(detection$get_contrasts())
      if (!setequal(intensity_contrasts$contrast, detection_contrasts$contrast)) {
        stop(
          "intensity and detection must have the same contrasts: ",
          paste(unique(intensity_contrasts$contrast), collapse = ", "),
          " vs ",
          paste(unique(detection_contrasts$contrast), collapse = ", ")
        )
      }
      detection_model <- unique(as.character(detection_contrasts$modelName))
      if (!all(detection_model %in% c("firth_nested", "binomial_nested", "binomial_firth_nested"))) {
        stop(
          "detection must be a firth_nested, binomial_nested or binomial_firth_nested facade, not: ",
          paste(detection_model, collapse = ", ")
        )
      }
      self$intensity <- intensity
      self$detection <- detection
      self$model_name <- model_name
      self$subject_id <- subject_id
      self$config <- intensity$get_config()
      private$intensity_contrasts <- intensity_contrasts
      private$detection_contrasts <- detection_contrasts
    },
    #' @description
    #' get contrast results with Fisher-combined \code{p.value} and \code{FDR}
    get_contrasts = function() {
      by <- c(self$subject_id, "contrast")
      intensity <- private$intensity_contrasts |>
        dplyr::select(-"modelName") |>
        dplyr::rename(p.value.intensity = "p.value", FDR.intensity = "FDR")
      detection <- private$detection_contrasts |>
        dplyr::select(
          dplyr::all_of(by),
          diff.detection = "diff",
          p.value.detection = "p.value",
          FDR.detection = "FDR"
        )
      res <- dplyr::full_join(intensity, detection, by = by) |>
        dplyr::left_join(.informative_detection(self$detection$.lfqdata), by = self$subject_id)

      p_detection <- dplyr::if_else(res$.detection_informative, res$p.value.detection, NA_real_)
      log_p <- cbind(log(res$p.value.intensity), log(p_detection))
      n_components <- rowSums(!is.na(log_p))
      res$chisq.hurdle <- dplyr::if_else(n_components > 0, -2 * rowSums(log_p, na.rm = TRUE), NA_real_)
      res$df.hurdle <- dplyr::if_else(n_components > 0, 2 * n_components, NA_real_)
      res$p.value <- stats::pchisq(res$chisq.hurdle, res$df.hurdle, lower.tail = FALSE)
      res$hurdle_components <- dplyr::case_when(
        n_components == 2 ~ "both",
        !is.na(res$p.value.intensity) ~ "intensity",
        !is.na(p_detection) ~ "detection"
      )
      res <- adjust_p_values(res) |>
        dplyr::ungroup() |>
        dplyr::select(-".detection_informative") |>
        dplyr::mutate(modelName = self$model_name, .before = 1) |>
        dplyr::relocate("estimate_type", .after = "modelName") |>
        dplyr::relocate(
          dplyr::all_of(c(
            "p.value.intensity",
            "FDR.intensity",
            "diff.detection",
            "p.value.detection",
            "FDR.detection",
            "chisq.hurdle",
            "df.hurdle",
            "hurdle_components"
          )),
          .after = dplyr::last_col()
        )
      stopifnot(all(super$column_description()$column_name %in% colnames(res)))
      res
    }
  ),
  private = list(
    intensity_contrasts = NULL,
    detection_contrasts = NULL
  )
)

#' Hurdle contrast facade for nested input
#'
#' Takes nested (peptide-level) LFQData and builds both hurdle components
#' internally: median polish to protein level and \code{lm_impute} on the protein
#' intensities, \code{binomial_firth_nested} on the peptide detection. The two are
#' combined with \code{\link{ContrastsHurdleFacade}}. Use
#' \code{ContrastsHurdleFacade} directly to pair other intensity or detection
#' facades.
#'
#' \code{get_missing()} reports against the median-polish protein data. The two
#' fitted models are \code{facade$contrast$intensity$model} and
#' \code{facade$contrast$detection$model}.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_peptide_config(Nprot = 20, weight_missing = 0.5, seed = 3)
#' lfq_peptide <- LFQData$new(istar$data, istar$config)
#' lfq_peptide <- lfq_peptide$get_Transformer()$log2()$lfq
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' facade <- build_contrast_analysis(lfq_peptide, "~ group_", contrasts, method = "hurdle_nested")
#' head(facade$get_contrasts())
ContrastsHurdleNestedFacade <- R6::R6Class(
  "ContrastsHurdleNestedFacade",
  inherit = ContrastsFacadeBase,
  public = list(
    #' @description
    #' initialize
    #' @param lfqdata nested LFQData (peptide-level, log2-transformed)
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param ... ignored; accepted for dispatch compatibility with other facades
    initialize = function(lfqdata, modelstr, contrasts, ...) {
      private$.setup(lfqdata, contrasts, "hurdle_nested", modelstr)
      lfq_protein <- lfqdata$get_Aggregator("medpolish")$aggregate()
      # get_missing() compares against the aggregated protein-level data.
      self$.lfqdata <- lfq_protein
      intensity <- ContrastsLMImputeFacade$new(lfq_protein, modelstr, contrasts)
      detection <- ContrastsBinomialFirthNestedFacade$new(lfqdata, modelstr, contrasts)
      self$contrast <- ContrastsHurdleFacade$new(intensity, detection, model_name = "hurdle_nested")
    }
  )
)

# Subjects with at least one undetected child feature in some sample; detection
# carries no information for the others.
.informative_detection <- function(lfqdata) {
  .summarize_detection_counts(.prepare_detection_lfqdata(lfqdata)) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(lfqdata$subject_id()))) |>
    dplyr::summarise(.detection_informative = any(.data$undetected > 0), .groups = "drop")
}
