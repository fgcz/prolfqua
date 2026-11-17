#' Remaining statistical software standards
#'
#' Standards 0.3.0, General and Regression. Current evidence and feasibility
#' are assessed separately in srr/standards-compliance.md. Both partial and unmet
#' requirements retain TODO tags; only fully met standards use srrstats.
#'
#' @srrstatsVerbose TRUE
#' @srrstatsTODO {G2.11} Validate/preprocess such columns intentionally and test the promised behavior; deferred
#'   follow-up: TODO/prolfqua/TODO_ropensci_compliance.md section 5.
#' @srrstatsTODO {G2.14a} Deferred NA-policy design remains necessary; deferred follow-up:
#'   TODO/prolfqua/TODO_ropensci_compliance.md section 5.
#' @srrstatsTODO {G3.1a} Add and document a covariance-policy API if accepted in a later design; deferred follow-up:
#'   TODO/prolfqua/TODO_ropensci_compliance.md section 5.
#' @srrstatsTODO {G5.6b} No repeated-seed recovery experiment for the same known generative parameters; deferred
#'   follow-up: TODO/prolfqua/TODO_ropensci_compliance.md section 5.
#' @srrstatsTODO {G5.9} Deterministic fixture construction is not a noise-susceptibility test of statistical results;
#'   deferred follow-up: TODO/prolfqua/TODO_ropensci_compliance.md section 5.
#' @srrstatsTODO {G5.9a} Add fixed-problem perturbation tests in a later statistical-testing plan; deferred follow-up:
#'   TODO/prolfqua/TODO_ropensci_compliance.md section 5.
#' @srrstatsTODO {G5.9b} No committed fixed-problem test establishes result stability under different RNG seeds or
#'   initial conditions; deferred follow-up: TODO/prolfqua/TODO_ropensci_compliance.md section 5.
#' @srrstatsTODO {G5.10} Add shared-flag extended statistical recovery/performance tests later; deferred follow-up:
#'   TODO/prolfqua/TODO_ropensci_compliance.md section 5.
#' @srrstatsTODO {G5.12} Document conditions when the deferred extended suite is implemented; deferred follow-up:
#'   TODO/prolfqua/TODO_ropensci_compliance.md section 5.
#' @srrstatsTODO {RE2.2} missing predictors not handled separately from missing responses; deferred follow-up:
#'   TODO/prolfqua/TODO_ropensci_compliance.md section 5.
#' @srrstatsTODO {RE2.4b} no check for perfect predictor-response fit (sigma = 0); deferred follow-up:
#'   TODO/prolfqua/TODO_ropensci_compliance.md section 5.
#' @srrstatsTODO {RE3.2} defaults undocumented; deferred follow-up: TODO/prolfqua/TODO_ropensci_compliance.md section
#'   5.
#' @srrstatsTODO {RE3.3} convergence thresholds not user-settable; deferred follow-up:
#'   TODO/prolfqua/TODO_ropensci_compliance.md section 5.
#' @srrstatsTODO {RE4.5} nobs per protein not exposed; deferred follow-up: TODO/prolfqua/TODO_ropensci_compliance.md
#'   section 5.
#' @srrstatsTODO {RE4.9} no fitted-values accessor; deferred follow-up: TODO/prolfqua/TODO_ropensci_compliance.md
#'   section 5.
#' @srrstatsTODO {RE4.10} no residuals accessor or documentation; deferred follow-up:
#'   TODO/prolfqua/TODO_ropensci_compliance.md section 5.
#' @srrstatsTODO {RE4.16} New-group prediction is required for categorical-group models independently of
#'   forecasting. No prediction method for new categorical groups is implemented; its semantics need a later design.
#' @srrstatsTODO {RE4.17} print method summarising model inputs and coefficients; deferred follow-up:
#'   TODO/prolfqua/TODO_ropensci_compliance.md section 5.
#' @srrstatsTODO {RE6.1} plot() dispatch directing users to get_Plotter(); deferred follow-up:
#'   TODO/prolfqua/TODO_ropensci_compliance.md section 5.
#' @srrstatsTODO {RE6.2} no plot of fitted values with confidence intervals; deferred follow-up:
#'   TODO/prolfqua/TODO_ropensci_compliance.md section 5.
#' @srrstatsTODO {RE7.1a} timing comparison noiseless vs noisy; deferred follow-up:
#'   TODO/prolfqua/TODO_ropensci_compliance.md section 5.
#' @noRd
NULL

#' NA_standards
#'
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
