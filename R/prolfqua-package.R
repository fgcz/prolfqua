#' @section Options:
#' \describe{
#'   \item{`prolfqua.verbose`}{`TRUE` (default) or `FALSE`. With `FALSE`,
#'     prolfqua prints no informative messages, such as columns added, rows
#'     removed or files written. Warnings and errors are not affected. The
#'     messages carry class `prolfqua_message`, so an application can route them
#'     to its own log with `withCallingHandlers(prolfqua_message = ...)`.}
#'   \item{`prolfqua.progress`}{Progress reporting while models and contrasts
#'     are computed per protein: `NULL` (default) shows a terminal progress
#'     bar, `"message"` a throttled message with elapsed time, and a
#'     `function(i, total, label)` is called with the progress.}
#' }
#' @section Lifecycle:
#' prolfqua is stable. Its exported API follows the Bioconductor deprecation
#' cycle: a renamed or removed interface calls \code{.Deprecated()} for one
#' release before calling \code{.Defunct()} in a subsequent release. Changes
#' affecting analyses are recorded in NEWS.
#'
#' @section Glossary:
#' \describe{
#'   \item{LFQData}{A long-format measurement table together with its
#'     \code{AnalysisConfiguration}; each row describes a feature in a sample.}
#'   \item{AnalysisConfiguration}{The mapping from table columns to hierarchy
#'     identifiers, sample identifiers, categorical factors and the response.}
#'   \item{Hierarchy and hierarchy depth}{Ordered measurement levels, such as
#'     protein, peptide and precursor. Depth selects the levels identifying a
#'     modelled subject; deeper levels identify its child features.}
#'   \item{subject_id}{The identifier columns of a modelled subject, usually
#'     one protein or a protein-peptide pair.}
#'   \item{Factor}{An explanatory variable defining experimental groups.
#'     \code{setup_analysis()} converts configured predictors to character
#'     labels, including numeric and ordered-factor inputs; numeric covariates
#'     and ordered-factor contrasts are not preserved through this setup.}
#'   \item{Response}{The numeric intensity column selected by the configuration,
#'     usually log2 abundance for modelling. On this scale a difference of one
#'     corresponds to a twofold abundance ratio.}
#'   \item{Contrast}{A named linear combination of group means or coefficients,
#'     for example treatment minus control.}
#'   \item{linfct}{A matrix with one row per requested linear function and one
#'     column per model coefficient. Multiplying it by the coefficient vector
#'     gives the estimates being tested.}
#'   \item{Strategy}{The formula, fitting algorithm and methods for extracting
#'     coefficients, residual variance, degrees of freedom and contrasts.}
#'   \item{Model}{Fitted regressions and their diagnostics, either per subject
#'     or in a matrix representation used by limma.}
#'   \item{Facade}{An analysis object composing a strategy, model and contrast
#'     calculation behind the same result and plotting interface.}
#'   \item{Aggregated and nested input}{Aggregated input has the same hierarchy
#'     keys as the output subject. Nested input also has child-feature keys;
#'     the analysis reports parent-level contrasts.}
#'   \item{nr_children}{The sample-wise number of measured child features,
#'     such as peptides used to estimate a protein intensity. It may be used
#'     as a fitting weight; DEqMS separately pools counts per protein.}
#'   \item{LOD and imputation}{The limit of detection on the working response
#'     scale, supplied explicitly or estimated from low-observation group
#'     means. Imputation replaces missing responses. The LOD methods refit
#'     otherwise unestimable subjects and borrow variance from successful fits.}
#'   \item{estimate_type}{A result label identifying how an estimate was
#'     obtained, for example from observed data or an imputed refit.}
#'   \item{Moderation}{Empirical Bayes shrinkage of residual variances towards
#'     a shared prior, optionally dependent on abundance or peptide count.
#'     Moderated degrees of freedom combine residual and prior information.}
#'   \item{Null hypothesis and p-value}{For a usual two-sided contrast test,
#'     the null states that the population contrast equals zero. Its p-value
#'     is the probability, assuming the null and model assumptions hold, of
#'     a statistic at least as extreme as the observed one. It is not the
#'     probability that the null hypothesis is true.}
#'   \item{FDR}{The false discovery rate: the expected proportion of false
#'     discoveries among rejected hypotheses. The default Benjamini-Hochberg
#'     adjustment is applied within contrast, assuming independent or suitably
#'     positively dependent tests. An adjusted p-value is not the probability
#'     that an individual protein is a false discovery.}
#'   \item{Confidence interval}{A range from a procedure with the stated
#'     repeated-sampling coverage under the model assumptions. Wald intervals
#'     use the estimated standard error and a t quantile; they do not assign
#'     a probability to a fixed parameter after observing the data.}
#'   \item{Linear-model assumptions}{On the working scale, classical linear
#'     models assume a linear mean, independent errors and constant residual
#'     variance; exact t inference also assumes normal errors. Weights allow
#'     specified relative variances and mixed models allow specified dependence
#'     within subjects. Robust, detection and imputation methods have
#'     different assumptions; diagnostics and the experimental design must
#'     guide the method choice.}
#' }
#' @keywords internal
"_PACKAGE"
