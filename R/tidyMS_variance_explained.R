#' Variance explained by the design, per feature
#'
#' Fits a least-squares model to every feature of a log-transformed [LFQData]
#' and splits each feature's variance with its ANOVA table. Use it to check a
#' normalization: a good normalization lowers the within-group variance more than
#' the overall variance, so the share explained by the design, the model R²,
#' rises.
#'
#' @param lfqdata a log-transformed [LFQData] object.
#' @param modelstr model formula string, as in [strategy_lm()]. `NULL` (default)
#'   uses the full factorial of `lfqdata$relevant_factor_keys()`, for example
#'   `"log2_abundance ~ Treatment * Background"`.
#' @return A tibble with one row per feature and ANOVA term, residuals included:
#'   the feature keys, `factor`, `Df`, `Sum.Sq`, `SS_total` (the feature's total
#'   sum of squares), `eta2` (`Sum.Sq / SS_total`) and `R2` (the feature's model
#'   R², `1 - SS_residual / SS_total`, repeated on each of its rows).
#' @details
#' One model is fitted per lowest-level feature (all hierarchy keys and the
#' isotope label), the grouping of [summarize_stats()]. Missing responses are
#' omitted within each fit. Only full-rank fits with more than one residual
#' degree of freedom are kept, the set returned by [get_complete_model_fit()],
#' so all features share the same terms; a message reports how many features
#' were left out.
#'
#' Read the model R² first. The ANOVA uses sequential (type I) sums of squares:
#' when a feature is unbalanced, for example through missing values, the shares
#' of the individual factors depend on the order of the terms in `modelstr`.
#' R² does not.
#'
#' Without any effect R² is positive on average: E\[R²\] = p / (n - 1) for p
#' model degrees of freedom and n observed values. Compare normalizations of the
#' same samples with the same missing values; imputation changes n.
#'
#' Use R² as a check of a normalization, not as a criterion to pick the
#' normalization with the highest value: that choice favours the group
#' differences tested afterwards.
#' @srrstats {RE4.11} summarize_variance_explained returns per-feature model R² and per-term eta² as goodness-of-fit
#'   statistics.
#' @seealso [pooled_variance_explained()]
#' @family stats
#' @export
#' @examples
#' bb <- prolfqua::sim_lfq_data_protein_config(Nprot = 20)
#' lfqdata <- LFQData$new(bb$data, bb$config)
#' # the simulated intensities are already on the log2 scale
#' lfqdata$is_transformed(TRUE)
#' ve <- summarize_variance_explained(lfqdata)
#' head(ve)
#' pooled_variance_explained(ve)
summarize_variance_explained <- function(lfqdata, modelstr = NULL) {
  if (!lfqdata$is_transformed()) {
    abort_bad_argument(
      "lfqdata",
      "be log-transformed; on raw intensities the multiplicative noise dominates the variance shares"
    )
  }
  if (is.null(modelstr)) {
    modelstr <- paste0(lfqdata$response(), " ~ ", paste(lfqdata$relevant_factor_keys(), collapse = " * "))
  }
  subject_id <- c(lfqdata$hierarchy_keys(), lfqdata$isotope_label())
  model <- build_model(lfqdata, strategy_lm(modelstr), subject_id = subject_id)
  fits <- get_complete_model_fit(model$model_df)
  n_features <- nrow(model$model_df)
  if (nrow(fits) == 0L) {
    abort_prolfqua(
      "No feature has a full-rank fit with more than one residual degree of freedom for ",
      modelstr,
      "."
    )
  }
  if (nrow(fits) < n_features) {
    inform_prolfqua(
      n_features - nrow(fits),
      " of ",
      n_features,
      " features left out: no full-rank fit with more than one residual degree of freedom."
    )
  }
  fits$variance <- lapply(fits$linear_model, .variance_table)
  fits |>
    dplyr::ungroup() |>
    dplyr::select(dplyr::all_of(c(subject_id, "variance"))) |>
    tidyr::unnest("variance")
}

#' ANOVA table of one fit as variance shares
#' @param fit a fitted `lm` object
#' @return data.frame with factor, Df, Sum.Sq, SS_total, eta2 and R2
#' @keywords internal
#' @noRd
.variance_table <- function(fit) {
  tab <- stats::anova(fit)
  term <- trimws(rownames(tab))
  ss <- tab[["Sum Sq"]]
  ss_total <- sum(ss)
  data.frame(
    factor = term,
    Df = tab[["Df"]],
    Sum.Sq = ss,
    SS_total = ss_total,
    eta2 = ss / ss_total,
    R2 = 1 - ss[term == "Residuals"] / ss_total,
    stringsAsFactors = FALSE
  )
}

#' Variance explained by each term, pooled over features
#'
#' Sums the sums of squares of [summarize_variance_explained()] over features.
#' For each term the result is `sum(Sum.Sq) / sum(SS_total)`; one minus the
#' share of the residuals is the pooled R². On features without missing values
#' the pooled R² equals the PERMANOVA R² with Euclidean distances between
#' samples. Features with a large total variance weigh more than in the median
#' of the per-feature `R2`.
#'
#' @param x the tibble returned by [summarize_variance_explained()].
#' @return A data.frame with columns `factor` and `eta2`, one row per term in
#'   the order of the ANOVA table, residuals included.
#' @seealso [summarize_variance_explained()]
#' @family stats
#' @export
#' @examples
#' bb <- prolfqua::sim_lfq_data_2factor_config(Nprot = 20)
#' bb$config$factor_depth <- 2
#' lfqdata <- LFQData$new(bb$data, bb$config)
#' lfqdata$is_transformed(TRUE)
#' pooled <- pooled_variance_explained(summarize_variance_explained(lfqdata))
#' pooled
#' # pooled R²
#' 1 - pooled$eta2[pooled$factor == "Residuals"]
pooled_variance_explained <- function(x) {
  missing_cols <- setdiff(c("factor", "Sum.Sq", "SS_total"), names(x))
  if (length(missing_cols) > 0L) {
    abort_missing_columns(missing_cols, "x")
  }
  terms <- unique(x$factor)
  eta2 <- vapply(
    terms,
    function(term) {
      rows <- x$factor == term
      sum(x$Sum.Sq[rows]) / sum(x$SS_total[rows])
    },
    numeric(1)
  )
  data.frame(factor = terms, eta2 = unname(eta2), stringsAsFactors = FALSE)
}
