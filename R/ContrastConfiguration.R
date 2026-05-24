# ContrastConfiguration ----

#' Contrast Configuration
#'
#' @description
#' Holds the column-role mapping and behaviour flags that describe how a
#' specific contrast backend (lm, limma, SAINT, ...) lays out its result
#' table. Mirrors \code{\link{AnalysisConfiguration}} for the modelling
#' side: it carries names, not data.
#'
#' Consumers that need to operate on a contrast table generically
#' (filter, build ORA inputs, build rank lists, summarise per protein)
#' should read column names through the accessors on this object rather
#' than hard-coding them. This lets backends with different column
#' conventions (e.g. SAINT's \code{BFDR}/\code{log2_EFCs}) plug into the
#' same downstream code path as LM-style backends without renaming
#' columns.
#'
#' Defaults match the standard prolfqua contrast schema produced by
#' \code{\link{Contrasts}}, \code{\link{ContrastsLimma}},
#' \code{\link{ContrastsModerated}} etc.
#'
#' @family configuration
#' @return An R6 class generator.
#' @export
#' @examples
#' cfg <- ContrastConfiguration$new(subject_id = "protein_Id")
#' cfg$effect_col
#' cfg$fdr_col
#' cfg$has_pvalue()
#' cfg$supports_dea_qc
#'
#' # SAINT-flavoured config
#' saint <- ContrastConfiguration$new(
#'   subject_id = "protein_Id",
#'   contrast_col = "Bait",
#'   effect_col = "log2_EFCs",
#'   score_col = "SaintScore",
#'   pvalue_col = NA_character_,
#'   fdr_col = "BFDR",
#'   supports_dea_qc = FALSE,
#'   needs_saint_annotation = TRUE
#' )
#' saint$has_pvalue()
ContrastConfiguration <- R6::R6Class(
  "ContrastConfiguration",
  public = list(
    #' @field subject_id columns identifying the contrast subject (e.g. protein id)
    subject_id = character(),
    #' @field model_name_col column carrying the model name label
    model_name_col = "modelName",
    #' @field contrast_col column carrying the contrast label
    contrast_col = "contrast",
    #' @field effect_col column with the (signed) effect size; for LM this
    #'   is \code{diff}, for SAINT \code{log2_EFCs}
    effect_col = "diff",
    #' @field score_col column with the per-contrast test statistic
    score_col = "statistic",
    #' @field pvalue_col column with the raw p-value, or \code{NA_character_}
    #'   for backends that do not produce one (e.g. SAINTexpress)
    pvalue_col = "p.value",
    #' @field fdr_col column with the FDR / BFDR / adjusted p-value
    fdr_col = "FDR",
    #' @field avg_abundance_col column with the per-feature mean abundance
    avg_abundance_col = "avgAbd",
    #' @field supports_dea_qc whether the backend supports the legacy
    #'   differential-expression QC HTML report
    supports_dea_qc = TRUE,
    #' @field needs_saint_annotation whether the backend needs annotation
    #'   reading in SAINT mode (bait / control columns)
    needs_saint_annotation = FALSE,
    #' @field significance_directional if \code{TRUE}, the backend only
    #'   considers positive-effect features biologically meaningful
    #'   (e.g. SAINTexpress, where negative \code{log2_EFCs} are not
    #'   interpreted as "down" hits). Affects
    #'   \code{filter_significant()}: when \code{TRUE} the filter is
    #'   one-sided on the effect column; when \code{FALSE} it is
    #'   symmetric (\code{|effect| > threshold}).
    significance_directional = FALSE,
    #' @description
    #' Construct a ContrastConfiguration. All arguments map directly to
    #' fields of the same name.
    #' @param subject_id columns identifying the contrast subject
    #' @param model_name_col model-name column
    #' @param contrast_col contrast-label column
    #' @param effect_col signed effect column
    #' @param score_col test-statistic column
    #' @param pvalue_col raw p-value column, or \code{NA_character_}
    #' @param fdr_col FDR / adjusted-p column
    #' @param avg_abundance_col mean-abundance column
    #' @param supports_dea_qc supports DEA QC HTML
    #' @param needs_saint_annotation needs SAINT-mode annotation
    #' @param significance_directional one-sided filter on effect column
    initialize = function(
      subject_id = character(),
      model_name_col = "modelName",
      contrast_col = "contrast",
      effect_col = "diff",
      score_col = "statistic",
      pvalue_col = "p.value",
      fdr_col = "FDR",
      avg_abundance_col = "avgAbd",
      supports_dea_qc = TRUE,
      needs_saint_annotation = FALSE,
      significance_directional = FALSE
    ) {
      self$subject_id <- subject_id
      self$model_name_col <- model_name_col
      self$contrast_col <- contrast_col
      self$effect_col <- effect_col
      self$score_col <- score_col
      self$pvalue_col <- pvalue_col
      self$fdr_col <- fdr_col
      self$avg_abundance_col <- avg_abundance_col
      self$supports_dea_qc <- supports_dea_qc
      self$needs_saint_annotation <- needs_saint_annotation
      self$significance_directional <- significance_directional
    },
    #' @description Does the backend produce a usable p-value column?
    has_pvalue = function() {
      pv <- self$pvalue_col
      length(pv) == 1 && !is.na(pv) && nzchar(pv)
    }
  )
)
