# ContrastsInterface ----
#' Base class for all Contrasts classes
#' @return An R6 class generator.
#' @export
#' @examples
#' int <- ContrastsInterface$new()
#' testthat::expect_error(int$get_contrast_sides())
#' testthat::expect_error(int$get_contrasts())
#' testthat::expect_error(int$get_missing())
#' testthat::expect_error(int$get_Plotter())
#' testthat::expect_error(int$to_wide())
#' # get_rank / get_ora / filter_significant call get_contrasts() internally,
#' # which is not implemented on the bare interface, so they surface that error.
#' testthat::expect_error(int$get_rank())
#' testthat::expect_error(int$get_ora())
#' testthat::expect_error(int$filter_significant())
#' identical(int$extra_artifacts(), list())
#' int$column_description()
ContrastsInterface <- R6::R6Class(
  "ContrastsInterface",
  public = list(
    #' @field config \code{\link{ContrastConfiguration}} describing the
    #'   backend's column-role mapping and behaviour flags. Subclasses
    #'   should populate this in \code{initialize()} so generic consumers
    #'   can resolve column names without hard-coding them.
    config = NULL,
    #' @description
    #' Return the \code{\link{ContrastConfiguration}}, constructing an
    #' LM-flavoured default on demand if a subclass has not set one.
    get_config = function() {
      if (is.null(self$config)) {
        subject_id <- if (length(self$subject_id) > 0) self$subject_id else character()
        self$config <- ContrastConfiguration$new(subject_id = subject_id)
      }
      self$config
    },
    #' @description
    #' get table with sides of the contrast
    get_contrast_sides = function() {
      stop("get_contrast_sides not implemented.")
    },
    #' @description
    #' get table with contrast results (similar to limma topTable function)
    get_contrasts = function() {
      stop("get_contrasts not implemented.")
    },
    #' @description
    #' initialize plotter
    get_Plotter = function() {
      stop("get_Plotter not implemented.")
    },
    #' @description
    #' create wide representation of data.
    to_wide = function() {
      stop("to_wide not implemented.")
    },
    #' @description
    #' get protein × contrast pairs that could not be estimated.
    #' Returns a data.frame with hierarchy columns and contrast, or 0 rows.
    get_missing = function() {
      stop("get_missing not implemented.")
    },
    #' @description
    #' get rank-list input table for enrichment tools.
    #'
    #' If the backend's \code{\link{ContrastConfiguration}} has a p-value
    #' column the default rank score is signed
    #' \code{sign(effect) * -log10(p.value)}. Otherwise the default is the
    #' effect column itself.
    #' @param score column name to use as rank score, or \code{NULL} to
    #'   pick from the config
    #' @return data.frame with subject id, contrast, and score columns.
    get_rank = function(score = NULL) {
      contrasts <- self$get_contrasts()
      cfg <- self$get_config()
      subject_id <- cfg$subject_id
      if (length(subject_id) == 0) {
        subject_id <- self$subject_id
      }
      contrast_col <- cfg$contrast_col
      effect_col <- cfg$effect_col
      pvalue_col <- cfg$pvalue_col
      if (is.null(score)) {
        if (cfg$has_pvalue() && all(c(effect_col, pvalue_col) %in% colnames(contrasts))) {
          contrasts$.rank_score <- sign(contrasts[[effect_col]]) *
            (-log10(contrasts[[pvalue_col]]))
          score_col_internal <- ".rank_score"
        } else {
          score_col_internal <- effect_col
        }
      } else {
        score_col_internal <- score
      }
      required <- c(subject_id, contrast_col, score_col_internal)
      missing_columns <- setdiff(required, colnames(contrasts))
      if (length(missing_columns) > 0) {
        stop(
          "Cannot create rank table. Missing columns: ",
          paste(missing_columns, collapse = ", ")
        )
      }
      rank_table <- contrasts |>
        dplyr::select(
          dplyr::all_of(subject_id),
          contrast = dplyr::all_of(contrast_col),
          score = dplyr::all_of(score_col_internal)
        ) |>
        dplyr::filter(!is.na(.data$score))
      return(rank_table)
    },
    #' @description
    #' get ORA input table — features passing the FDR and absolute
    #' effect-size threshold in the requested direction.
    #' @param up if TRUE return positive effects, otherwise negative effects
    #' @param FDR_threshold false discovery rate threshold
    #' @param diff_threshold absolute effect-size threshold
    #' @return filtered contrast data.frame for ORA input.
    get_ora = function(up = TRUE, FDR_threshold = 0.05, diff_threshold = 1) {
      contrasts <- self$get_contrasts()
      cfg <- self$get_config()
      subject_id <- cfg$subject_id
      if (length(subject_id) == 0) {
        subject_id <- self$subject_id
      }
      contrast_col <- cfg$contrast_col
      effect_col <- cfg$effect_col
      fdr_col <- cfg$fdr_col
      required <- c(subject_id, contrast_col, fdr_col, effect_col)
      missing_columns <- setdiff(required, colnames(contrasts))
      if (length(missing_columns) > 0) {
        stop(
          "Cannot create ORA table. Missing columns: ",
          paste(missing_columns, collapse = ", ")
        )
      }
      if (up) {
        ora <- contrasts[
          contrasts[[fdr_col]] < FDR_threshold &
            !is.na(contrasts[[fdr_col]]) &
            contrasts[[effect_col]] > diff_threshold &
            !is.na(contrasts[[effect_col]]),
        ]
      } else {
        ora <- contrasts[
          contrasts[[fdr_col]] < FDR_threshold &
            !is.na(contrasts[[fdr_col]]) &
            contrasts[[effect_col]] < -diff_threshold &
            !is.na(contrasts[[effect_col]]),
        ]
      }
      return(ora)
    },
    #' @description
    #' Filter contrast rows that pass the FDR and effect-size threshold.
    #'
    #' By default this is symmetric on the effect column
    #' (\code{|effect| > diff_threshold}). For backends whose
    #' \code{\link{ContrastConfiguration}} has
    #' \code{significance_directional = TRUE} (e.g. SAINTexpress, where
    #' only positive enrichment is biologically meaningful), the filter
    #' is one-sided (\code{effect > diff_threshold}).
    #' @param FDR_threshold false discovery rate threshold
    #' @param diff_threshold effect-size threshold
    #' @return filtered contrast data.frame.
    filter_significant = function(FDR_threshold = 0.05, diff_threshold = 1) {
      contrasts <- self$get_contrasts()
      cfg <- self$get_config()
      fdr_col <- cfg$fdr_col
      effect_col <- cfg$effect_col
      required <- c(fdr_col, effect_col)
      missing_columns <- setdiff(required, colnames(contrasts))
      if (length(missing_columns) > 0) {
        stop(
          "Cannot filter significant contrasts. Missing columns: ",
          paste(missing_columns, collapse = ", ")
        )
      }
      effect_predicate <- if (isTRUE(cfg$significance_directional)) {
        contrasts[[effect_col]] > diff_threshold
      } else {
        abs(contrasts[[effect_col]]) > diff_threshold
      }
      keep <- contrasts[[fdr_col]] < FDR_threshold &
        !is.na(contrasts[[fdr_col]]) &
        effect_predicate &
        !is.na(contrasts[[effect_col]])
      contrasts[keep, , drop = FALSE]
    },
    #' @description
    #' Per-subject summary table for downstream report grobs. Returns a
    #' small data.frame with canonical column names
    #' \code{contrast}, \code{effect}, \code{score}, \code{fdr}
    #' regardless of backend.
    #' @param rounded if TRUE round numeric columns to 3 significant
    #'   digits (matches the prolfquapp report convention)
    contrast_summary_table = function(rounded = TRUE) {
      contrasts <- self$get_contrasts()
      cfg <- self$get_config()
      contrast_col <- cfg$contrast_col
      effect_col <- cfg$effect_col
      score_col <- cfg$score_col
      fdr_col <- cfg$fdr_col
      required <- c(contrast_col, effect_col, score_col, fdr_col)
      missing_columns <- setdiff(required, colnames(contrasts))
      if (length(missing_columns) > 0) {
        stop(
          "Cannot build contrast summary. Missing columns: ",
          paste(missing_columns, collapse = ", ")
        )
      }
      res <- data.frame(
        contrast = contrasts[[contrast_col]],
        effect = contrasts[[effect_col]],
        score = contrasts[[score_col]],
        fdr = contrasts[[fdr_col]],
        stringsAsFactors = FALSE
      )
      if (isTRUE(rounded)) {
        round3 <- function(x) if (is.numeric(x)) signif(x, 3) else x
        res$effect <- round3(res$effect)
        res$score <- round3(res$score)
        res$fdr <- round3(res$fdr)
      }
      res
    },
    #' @description
    #' Backend-specific extra artifacts that a report should carry
    #' along with the contrast results (e.g. SAINT input tables).
    #' Default returns an empty named list.
    extra_artifacts = function() {
      list()
    },
    #' @description
    #' column description
    column_description = function() {
      description <- c(
        "modelName" = "type of model",
        "contrast" = "name of difference e.g. group1_vs_group2",
        "avgAbd" = "mean abundance value of protein in all samples",
        "diff" = "difference among conditions",
        "FDR" = "false discovery rate",
        "statistic" = "t-statistics",
        "std.error" = "standard error",
        "df" = "degrees of freedom",
        "p.value" = "p-value",
        "conf.low" = "lower value of 95 confidence interval",
        "conf.high" = "high value of 95 confidence interval",
        "sigma" = "residual standard deviation of linear model (needed for empirical Bayes variance shrinkage)."
      )
      description <- data.frame(column_name = names(description), description = description)
      return(description)
    }
  )
)


# Merge contrasts ----
#' Merge contrast results coming from two different model.
#'
#' Typically used with results of \code{\link{Contrasts}} and \code{\link{ContrastsMissing}}
#'
#' @param prefer contrasts to use preferentially
#' @param add contrasts to add from if missing in prefer
#' @param model_name name of the merged model default "mergedModel"
#' @return Contrast definitions or contrast results.
#' @export
#' @family modelling
#'
#' @examples
#' prefer_df <- data.frame(
#'   protein_Id = c("P1", "P2"),
#'   contrast = "A_vs_B",
#'   modelName = "prefer",
#'   diff = c(1, 2),
#'   p.value = c(0.01, 0.2),
#'   FDR = c(0.02, 0.2),
#'   statistic = c(3, NA)
#' )
#' add_df <- data.frame(
#'   protein_Id = c("P1", "P2"),
#'   contrast = "A_vs_B",
#'   modelName = "add",
#'   diff = c(1.1, 2.1),
#'   p.value = c(0.02, 0.03),
#'   FDR = c(0.03, 0.04),
#'   statistic = c(2.8, 2.5)
#' )
#' prefer <- ContrastsTable$new(prefer_df, subject_id = "protein_Id", model_name = "prefer")
#' add <- ContrastsTable$new(add_df, subject_id = "protein_Id", model_name = "add")
#' merged <- merge_contrasts_results(prefer, add)
#' merged$merged$get_contrasts()
merge_contrasts_results <- function(prefer, add, model_name = "mergedModel") {
  c_a <- prefer$get_contrasts()
  c_b <- add$get_contrasts()

  if (length(colnames(c_a)) < length(colnames(c_b))) {
    c_b <- dplyr::select(c_b, colnames(c_a))
  }

  c_a <- dplyr::filter(c_a, !is.na(.data$statistic))
  more_id <- setdiff(
    distinct(select(c_b, c(prefer$subject_id, "contrast"))),
    distinct(select(c_a, c(add$subject_id, "contrast")))
  )
  more <- inner_join(more_id, c_b)

  same_id <- select(c_a, c(add$subject_id, "contrast"))
  same <- inner_join(same_id, c_b)

  merged <- bind_rows(c_a, more)

  if (prefer$model_name == add$model_name) {
    prefer_model_name <- paste0(prefer$model_name, "_prefer")
    add_model_name <- paste0(add$model_name, "_add")
    c_a$modelName <- prefer_model_name
    more$modelName <- add_model_name
  } else {
    prefer_model_name <- prefer$model_name
    add_model_name <- add$model_name
  }

  merged$modelName <- factor(merged$modelName, levels = c(levels(factor(c_a$modelName)), add_model_name))

  merged <- ContrastsTable$new(
    merged,
    subject_id = prefer$subject_id,
    model_name = paste0(prefer_model_name, "_", add_model_name)
  )
  more <- ContrastsTable$new(more, subject_id = prefer$subject_id, model_name = add_model_name)
  same <- ContrastsTable$new(same, subject_id = prefer$subject_id, model_name = add_model_name)
  return(list(merged = merged, more = more, same = same))
}
