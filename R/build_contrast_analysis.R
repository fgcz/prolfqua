# build_contrast_analysis -----

#' Build a contrast analysis using one of several statistical methods
#'
#' The factory of prolfqua's modelling: \code{method} selects a facade class
#' from the facade registry, and the facade runs the full pipeline from
#' strategy construction through modelling to contrast computation. Every
#' facade returns its results through the same interface
#' (\code{get_contrasts()}, \code{get_Plotter()}, \code{get_missing()},
#' \code{to_wide()}), so the results of different methods have the same
#' columns and can be compared or combined.
#'
#' Methods ending in \code{_nested} take peptide- or precursor-level input and
#' report protein-level contrasts; the others take input at the level they
#' report, usually aggregated proteins. \code{"limpa"} needs input aggregated
#' with \code{\link{AggregateLimpa}}.
#'
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param modelstr one right-hand-side formula string without the response variable
#'   (e.g. \code{"~ group_"}). The response is taken automatically from
#'   \code{lfqdata$get_config()$get_response()}.
#' @param contrasts nonempty named character vector of contrast expressions
#'   (e.g. \code{c("A_vs_B" = "group_A - group_B")})
#' @param method a registered facade key. The built-in keys are:
#'   linear models \code{"lm"}, \code{"rlm"}, \code{"rfit"};
#'   limma \code{"limma"}, \code{"limma_voom"}, \code{"deqms"},
#'   \code{"deqms_voom"}, \code{"limpa"}, \code{"limpa_nested"};
#'   limit-of-detection imputation \code{"lm_impute"}, \code{"limma_impute"},
#'   \code{"limma_voom_impute"}, \code{"rfit_impute"};
#'   peptide-level models \code{"lmer_nested"}, \code{"ropeca_nested"},
#'   \code{"ropecaV2_nested"};
#'   detection models \code{"firth"}, \code{"binomial_nested"},
#'   \code{"binomial_firth_nested"};
#'   the hurdle model \code{"hurdle_nested"}, combining intensity and detection.
#'   Downstream packages may add more via \code{\link{register_facade}}. The
#'   Keys and column names are case sensitive and must match exactly. The
#'   authoritative list is \code{names(\link{list_facades}())}. Defaults to
#'   \code{"lm"}.
#' @param ... additional arguments forwarded to the underlying strategy function
#'   (e.g. \code{trend}, \code{robust} for \code{strategy_limma})
#' @details The formula is translated into a design matrix with R's current
#'   contrast options (normally treatment coding). Configured predictors are
#'   categorical labels; their reference level is normally the first sorted
#'   level. Each contrast expression combines group-level coefficient weights,
#'   rather than referring directly to design-matrix coefficient names. For
#'   example, \code{"group_A - group_Ctrl"} subtracts the control group mean
#'   from the A group mean, including the reference coefficient as needed.
#'   Interaction cell means are named with a colon and quoted with backticks
#'   inside expressions. See the \emph{Under the hood} article for worked
#'   one-factor and \code{~ A * B} design matrices and linear functions.
#'
#'   Internally, \code{limpa_nested} and \code{hurdle_nested} aggregate child
#'   features to proteins. Detection methods (\code{binomial_nested},
#'   \code{binomial_firth_nested}, \code{firth}, and the detection component
#'   of \code{hurdle_nested}) encode intensities as detected/not detected.
#'   These aggregation and detection steps have no inverse. The analysis keeps
#'   its input in \code{analysis$.lfqdata}; for aggregated input this is already
#'   at the modelled hierarchy. \code{fitted()} and \code{residuals()} expose
#'   the modelled observed responses (inspect each hurdle component separately).
#'   Intensity effects are differences on the modelled scale, hence log2 fold
#'   changes for log2 responses; detection effects are log odds ratios, with
#'   fitted responses reported as probabilities. Transformations performed
#'   before this call retain their input columns in LFQData; this call does not
#'   record a transformation history or infer inverse transformations.
#'
#'   With a fixed design and fixed observations per protein, the per-protein
#'   \code{lm} and \code{lmer_nested} fitting loops grow approximately linearly
#'   with the number of proteins and run sequentially. \code{limma} fits the
#'   expression matrix in one call and has lower per-protein overhead.
#'   Single-run elapsed seconds for 500, 1000 and 2000 proteins were respectively
#'   0.173, 0.288 and 0.592 for \code{lm}; 0.077, 0.043 and 0.050 for
#'   \code{limma}; and 7.902, 15.314 and 29.667 for \code{lmer_nested}
#'   (R 4.6.2, arm64 macOS, 12 samples, three peptides per protein for the
#'   nested model, no missing responses). These timings cover facade
#'   construction, including fitting, but exclude lazy contrast extraction.
#'   They are illustrative timings, not guarantees: sample count, design width,
#'   missingness and mixed-model convergence affect cost. Some mixed fits in
#'   this run emitted convergence warnings. Reproduction code and environment
#'   details are in \code{srr/benchmark-runtime.R} in the source repository.
#' @section Prior art:
#' \itemize{
#'   \item Improvements on existing R implementations, applied per protein behind one interface:
#'     \code{lm} (stats), \code{rlm} (MASS), \code{rfit} (Rfit), \code{limma} and \code{limma_voom} (limma),
#'     \code{deqms} (DEqMS), \code{limpa} and \code{limpa_nested} (limpa), \code{lmer_nested} (lme4, lmerTest),
#'     \code{firth} (logistf).
#'   \item Reimplementation in R of a published method: \code{ropeca_nested}, ROPECA (Suomi and Elo 2017,
#'     Bioconductor package PECA), on peptide-level models of any design.
#'   \item First implementations of new methods: the \code{_impute} methods (\code{lm_impute},
#'     \code{rfit_impute}, \code{limma_impute}, \code{limma_voom_impute}), which refit groups without
#'     observations at the limit of detection with a borrowed variance (Wolski et al. 2023);
#'     \code{ropecaV2_nested}, a median peptide z-score with a correlation-aware null; \code{binomial_nested} and
#'     \code{binomial_firth_nested}, per-protein models of peptide detection.
#'   \item \code{hurdle_nested} combines the intensity and detection tests as msqrob2's
#'     \code{hypothesisTestHurdle} does, on prolfqua's own components.
#' }
#' The Comparing methods article on the package website compares these
#' implementations on one dataset. Comparisons with MSstats, msqrob2 and proDA are in the companion package
#' prolfquabenchmark.
#' @return a facade object of the class registered for \code{method}, see
#'   \code{\link{lookup_facade}}.
#' @export
#' @srrstats {RE4.12} Modelling aggregation and detection encoding, their retained responses and lack of inverses
#'   are documented; preceding user transformations retain their input columns.
#' @family modelling
#' @examplesIf requireNamespace("lme4", quietly = TRUE) && requireNamespace("lmerTest", quietly = TRUE)
#' istar <- sim_lfq_data_protein_config(n_proteins = 20)
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
#' @srrstats {G1.1} Section "Prior art" classifies every method as improvement on an R implementation,
#'   reimplementation in R, or first implementation of a new method.
#' @srrstats {G1.6} The ComparingMethods article compares implementations from other R packages (limma, DEqMS, limpa,
#'   MASS, Rfit, lme4, logistf) wrapped by prolfqua on one dataset; MSstats, msqrob2 and proDA comparisons are in
#'   prolfquabenchmark.
#' @srrstatsTODO {G2.1a} The builder documents LFQData, formula strings and named contrasts; limit: types are not
#'   exhaustive API-wide.
#' @srrstats {G2.3b} Builder help explicitly requires case-sensitive method keys and column names.
#' @srrstats {RE1.0} The builder takes a validated RHS formula and strategies retain the full response formula.
#' @srrstats {RE2.1} `method` is the missing-value control (drop per fit, LOD imputation, detection models); NA is a
#'   missing observation, Inf and -Inf raise. Limit: NaN raises with Inf instead of being treated as NA.
#' @srrstats {RE5.0} Builder details report measured scaling for three methods; source benchmark states design and
#'   timing limits.
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
  if (!inherits(lfqdata, "LFQData")) {
    abort_bad_argument("lfqdata", "be an LFQData object")
  }
  method <- .arg_match_prolfqua(method, choices, "method")
  entry <- lookup_facade(method)
  # Decoys must not enter the fit or the shared variance pool (limma prior /
  # DEqMS variance-count trend), where they would perturb *target* q-values.
  # Opt-in: only when the config declares `pattern_decoys` (NULL leaves existing
  # behaviour untouched). Normalization is expected to have already run on the
  # full data upstream; here we drop decoys so the model sees targets only.
  cfg <- lfqdata$get_config()
  if (!is.null(cfg$pattern_decoys)) {
    # remove_decoys() records and reports the dropped decoys
    lfqdata <- lfqdata$remove_decoys()
  }
  facade_class <- utils::getFromNamespace(entry$class, entry$package)
  facade_class$new(lfqdata, modelstr, contrasts, ...)
}
