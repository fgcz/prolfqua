# ContrastsComparison -----

#' Compare the contrasts of several methods on their subject and contrast pairs
#'
#' Aligns the contrast results of several methods, e.g. facades of [build_contrast_analysis()] or [ContrastsTable]s
#' holding another tool's results, on the subject (e.g. protein) and the contrast, and reports which pairs each method
#' estimates and where their inferences contradict each other. The columns of each result are read through its
#' [ContrastConfiguration], so tables with other column names can be compared too.
#'
#' The status of a pair in a method is `not estimated` (no row or no effect), `not tested` (an effect but no FDR),
#' `not significant`, or `up` / `down` (FDR below `fdr_threshold`, by the sign of the effect).
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' sim <- sim_lfq_data_protein_config(n_proteins = 20, weight_missing = 0.5, seed = 3)
#' lfqdata <- LFQData$new(sim$data, sim$config)
#' contrasts <- c(A_vs_Ctrl = "group_A - group_Ctrl", B_vs_Ctrl = "group_B - group_Ctrl")
#' methods <- c("lm", "limma", "lm_impute")
#' fits <- lapply(stats::setNames(methods, methods), function(method) {
#'   build_contrast_analysis(lfqdata, "~ group_", contrasts, method = method)
#' })
#' comparison <- ContrastsComparison$new(fits)
#' head(comparison$aligned())
#' comparison$coverage()
#' conflicts <- comparison$conflicts()
#' conflicts[conflicts$coverage_conflict | conflicts$significance_conflict, ]
#' head(comparison$pairs("lm", columns = c("effect", "p.value")))
#' if (requireNamespace("UpSetR", quietly = TRUE)) comparison$upset()
ContrastsComparison <- R6::R6Class(
  "ContrastsComparison",
  public = list(
    #' @field subject_id columns identifying the subjects, shared by all methods
    subject_id = character(),
    #' @field fdr_threshold FDR below which a contrast is significant
    fdr_threshold = 0.05,
    #' @description
    #' align the contrast results of the methods
    #' @param contrasts named list of contrast objects (facades, [ContrastsInterface] objects, [ContrastsTable]s); the
    #'   names label the methods
    #' @param fdr_threshold FDR below which a contrast is significant
    initialize = function(contrasts, fdr_threshold = 0.05) {
      if (!is.numeric(fdr_threshold) || length(fdr_threshold) != 1 || !(fdr_threshold > 0 && fdr_threshold < 1)) {
        abort_bad_argument("fdr_threshold", "be one number between 0 and 1")
      }
      self$subject_id <- .comparison_subject_id(contrasts)
      self$fdr_threshold <- fdr_threshold
      private$methods <- names(contrasts)
      tables <- Map(.comparison_table, contrasts, names(contrasts))
      private$long <- .align_comparison_tables(tables, self$subject_id, fdr_threshold)
    },
    #' @description
    #' one row per subject, contrast and method over all pairs any method reports: `effect`, `statistic`, `p.value`,
    #' `FDR`, `estimate_type` (NA for tables without it) and `status`
    #' @return data.frame
    aligned = function() {
      private$long
    },
    #' @description
    #' the pairs of one method with all the given columns defined, e.g. `c("effect", "p.value")` for the pairs with a
    #' fold change and a p-value
    #' @param method method label, one of the names of `contrasts`
    #' @param columns one or more of `"effect"`, `"statistic"`, `"p.value"`
    #' @return data.frame with the subject columns and `contrast`
    pairs = function(method, columns = "effect") {
      method <- match.arg(method, private$methods)
      columns <- match.arg(columns, c("effect", "statistic", "p.value"), several.ok = TRUE)
      rows <- private$long[private$long$method == method, ]
      rows[stats::complete.cases(rows[columns]), c(self$subject_id, "contrast")]
    },
    #' @description
    #' number of pairs per set of methods that estimate them (`none`: reported, but without an effect)
    #' @return data.frame with columns `estimated_by` and `pairs`
    coverage = function() {
      private$long |>
        dplyr::summarize(
          estimated_by = paste(.data$method[.data$status != "not estimated"], collapse = " & "),
          .by = dplyr::all_of(c(self$subject_id, "contrast"))
        ) |>
        dplyr::mutate(estimated_by = dplyr::if_else(.data$estimated_by == "", "none", .data$estimated_by)) |>
        dplyr::count(.data$estimated_by, name = "pairs", sort = TRUE)
    },
    #' @description
    #' per pair, the number of methods by status, and where the methods contradict each other:
    #' `sign_conflict` (up in one method, down in another), `significance_conflict` (significant in one, tested but
    #' not significant in another), `coverage_conflict` (significant in one, not estimated by another)
    #' @return data.frame, one row per pair
    conflicts = function() {
      private$long |>
        dplyr::summarize(
          up = sum(.data$status == "up"),
          down = sum(.data$status == "down"),
          not_significant = sum(.data$status == "not significant"),
          not_tested = sum(.data$status == "not tested"),
          not_estimated = sum(.data$status == "not estimated"),
          .by = dplyr::all_of(c(self$subject_id, "contrast"))
        ) |>
        dplyr::mutate(
          sign_conflict = .data$up > 0 & .data$down > 0,
          significance_conflict = .data$up + .data$down > 0 & .data$not_significant > 0,
          coverage_conflict = .data$up + .data$down > 0 & .data$not_estimated > 0
        )
    },
    #' @description
    #' UpSet plot of the pairs each method estimates
    #' @return the plot of [UpSetR::upset()]
    upset = function() {
      rlang::check_installed("UpSetR", reason = "for UpSet plots.")
      estimated <- private$long[private$long$status != "not estimated", ]
      keys <- do.call(paste, c(estimated[c(self$subject_id, "contrast")], sep = "~"))
      sets <- split(keys, factor(estimated$method, levels = private$methods))
      UpSetR::upset(UpSetR::fromList(sets), nsets = length(sets), order.by = "freq")
    }
  ),
  private = list(
    long = NULL,
    methods = character()
  )
)

#' Check the arguments of ContrastsComparison and return the subject columns all methods share
#' @param contrasts named list of contrast objects
#' @noRd
.comparison_subject_id <- function(contrasts) {
  # unique, non-empty names, one per element (missing names give none)
  labels <- unique(names(contrasts))
  if (!is.list(contrasts) || length(contrasts) < 2 || length(labels[nzchar(labels)]) != length(contrasts)) {
    abort_bad_argument("contrasts", "be a list of at least two contrast objects with unique names")
  }
  subject_ids <- unique(lapply(contrasts, function(contrast) contrast$get_config()$subject_id))
  if (length(subject_ids) != 1) {
    abort_bad_argument(
      "contrasts",
      "share their subject columns",
      not = toString(vapply(subject_ids, function(x) paste0("(", toString(x), ")"), ""))
    )
  }
  subject_ids[[1]]
}

#' One row per subject, contrast and method over all pairs any method reports, with the status of the pair
#' @param tables one table per method from .comparison_table()
#' @param subject_id subject columns
#' @param fdr_threshold FDR below which a contrast is significant
#' @noRd
.align_comparison_tables <- function(tables, subject_id, fdr_threshold) {
  contrast_names <- unique(lapply(tables, function(table) sort(unique(table$contrast))))
  if (length(contrast_names) != 1) {
    abort_bad_argument(
      "contrasts",
      "have the same contrast names",
      not = toString(vapply(contrast_names, function(x) paste0("(", toString(x), ")"), ""))
    )
  }
  rows <- dplyr::bind_rows(tables)
  pairs <- dplyr::distinct(rows[c(subject_id, "contrast")])
  aligned <- dplyr::left_join(
    tidyr::expand_grid(pairs, method = names(tables)),
    rows,
    by = c(subject_id, "contrast", "method"),
    relationship = "one-to-one"
  )
  aligned$status <- .contrast_status(aligned$effect, aligned$FDR, fdr_threshold)
  aligned
}

#' One method's contrasts with the columns named by role, read through its ContrastConfiguration
#' @param contrast contrast object with `get_contrasts()` and `get_config()`
#' @param method method label
#' @noRd
.comparison_table <- function(contrast, method) {
  config <- contrast$get_config()
  result <- contrast$get_contrasts()
  # a backend without p-values has the p-value column ""
  roles <- c(
    effect = config$effect_col,
    statistic = config$score_col,
    p.value = config$pvalue_col,
    FDR = config$fdr_col
  )
  .stop_if_missing_columns(
    result,
    c(config$subject_id, config$contrast_col, roles[nzchar(roles)]),
    paste0("compare the contrasts of `", method, "`")
  )
  table <- result[config$subject_id]
  table$contrast <- result[[config$contrast_col]]
  table$method <- method
  for (role in names(roles)) {
    table[[role]] <- if (nzchar(roles[[role]])) result[[roles[[role]]]] else NA_real_
  }
  table$estimate_type <- if ("estimate_type" %in% colnames(result)) result$estimate_type else NA_character_
  table
}

#' Status of a pair in a method: not estimated, not tested, not significant, up or down
#' @param effect signed effect sizes
#' @param fdr false discovery rates
#' @param fdr_threshold FDR below which a contrast is significant
#' @noRd
.contrast_status <- function(effect, fdr, fdr_threshold) {
  dplyr::case_when(
    is.na(effect) ~ "not estimated",
    is.na(fdr) ~ "not tested",
    fdr < fdr_threshold & effect > 0 ~ "up",
    fdr < fdr_threshold & effect < 0 ~ "down",
    .default = "not significant"
  )
}
