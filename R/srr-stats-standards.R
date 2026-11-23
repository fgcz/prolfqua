#' Remaining statistical software standards
#'
#' Standards 0.3.0, General and Regression. Current evidence and feasibility
#' are assessed separately in srr/standards-compliance.md. Both partial and unmet
#' requirements retain TODO tags; only fully met standards use srrstats.
#'
#' @srrstatsVerbose TRUE
#' @srrstats {G2.11} Columns with non-standard classes are processed by their values: the response loses extra
#'   classes and attributes (`I()`, units, labelled) and becomes double; factor and hierarchy columns of class Date or
#'   `I()` are used by their labels; list columns in configured roles are refused (test-input-types.R,
#'   test-eda-output.R). Limit: the units package itself is not used in the tests; a response carrying a `units`
#'   attribute is.
#' @srrstats {RE2.2} Missing predictors (hierarchy, factor and sample keys) raise; missing responses are dropped per fit
#'   or imputed, and `impute_from_model()` generates values for missing responses from models fitted without them.
#' @srrstats {RE3.2} Convergence defaults are documented with their reasons in strategy_rlm(), strategy_lmer(),
#'   strategy_binomial(), strategy_binomial_firth(), strategy_logistf(), strategy_rfit() and estimate_intensity():
#'   the criteria are relative (rlm, glm, rfit, medpolish on the log2 scale) or in log odds (Firth). limpa's iterative
#'   steps (dpc(), dpcQuant()) keep limpa's defaults, documented on AggregateLimpa; limpa reports no convergence
#'   status, so the fit status of limpa analyses gives `not_converged` NA.
#' @noRd
NULL

#' NA_standards
#'
#' @srrstatsNA {G2.14a} NA intensities are missing observations in LFQ data; an option to error on them is not
#'   meaningful. Non-finite values and NA keys raise errors (`.assert_model_input()`).
#' @srrstatsNA {G5.4c} Relevant prior implementations are available as software and can be called directly.
#'   Published-paper stored numerical values are not needed as a substitute for unavailable code.
#' @srrstatsNA {G5.11} The current suite uses packaged datasets and generated simulations; no extended downloaded
#'   assets are required. Conditional download requirement is inapplicable to the existing test workflow.
#' @srrstatsNA {G5.11a} No test workflow downloads extra assets. No download-failure skip policy is necessary unless
#'   such a workflow is added.
#' @srrstatsNA {RE4.14} No public extrapolation/forecasting interface is implemented. Statistical purpose is inference
#'   on contrasts among experimental design groups; internal predict() is used only to fill missing responses for an
#'   existing fitted design. -
#' @srrstatsNA {RE4.15} No extrapolation/forecasting horizon is part of the public workflow; there is no forecast
#'   interval whose growth could be characterized. -
#' @srrstatsNA {RE6.3} No public forecasting/extrapolation workflow exists; model-versus-forecast display is outside
#'   package scope. -
#' @srrstatsNA {RE7.4} Forecast horizons and forecast error growth are outside the package public workflow, as for
#'   RE4.15. -
#' @noRd
NULL
