# Moderation and p-value adjustment ----

#' Moderate p-values - limma approach
#' @param contrast_df contrast table with the columns `sigma`, the residual degrees of freedom and the estimate
#' @param df name of the residual degrees-of-freedom column
#' @param estimate name of the effect-estimate column
#' @param robust passed to [limma::squeezeVar()]: estimate the prior robustly
#' @param confint confidence level of the moderated intervals
#' @param variance_floor optional positive lower bound for posterior variances
#' @return `contrast_df`, ungrouped, with the added columns `moderated.var.post`, `moderated.var.prior`,
#'   `moderated.df.prior`, `moderated.std.error`, `moderated.statistic`, `moderated.df.total`,
#'   `moderated.p.value`, `moderated.conf.low` and `moderated.conf.high`
#' @srrstats {G5.4} Tests compare moderation to squeezeVar, detection to glm/logistf and Hurdle to Fisher combination.
#'   test-power_fdr.R compares FDR planning with power.t.test, analytic normal mixing and ssize.fdr.
#' @srrstats {G5.4b} test-reference-implementations.R compares the re-implemented methods with the original software:
#'   lm, limma, lmer, rlm, rfit and deqms with stats::lm, limma, lmerTest, MASS, Rfit and stored DEqMS output. ROPECA
#'   is a derivation (G5.4a), not a re-implementation.
#' @noRd
.moderated_p_limma <- function(
  contrast_df,
  df = "df",
  estimate = "diff",
  robust = FALSE,
  confint = 0.95,
  variance_floor = NULL
) {
  valid_floor <- is.numeric(variance_floor) && length(variance_floor) == 1L && is.finite(variance_floor)
  if (!is.null(variance_floor) && !(valid_floor && variance_floor > 0)) {
    abort_bad_argument("variance_floor", "be NULL or one positive number")
  }

  # Empirical-Bayes variance moderation via limma. Since limma 3.x, squeezeVar()
  # estimates the prior robustly on fractional / low residual df (as produced by
  # rlm, lmer_nested and the imputation facades) without the min_df workaround
  # that prolfqua's former squeezeVarRob() fork carried over from MSqRob.
  if (all(is.na(contrast_df$sigma))) {
    # no subject has a variance (all failed or with zero residual variance): no prior, nothing to moderate
    squeezed_var <- list(var.post = rep(NA_real_, nrow(contrast_df)), var.prior = NA_real_, df.prior = NA_real_)
  } else {
    squeezed_var <- limma::squeezeVar(contrast_df$sigma^2, df = contrast_df[[df]], robust = robust)
  }
  if (all(is.infinite(squeezed_var$df.prior))) {
    # count only rows with a variance; NA variances (zero residual variance, failed fits) are not in the prior
    squeezed_var$df.prior <- mean(contrast_df[[df]], na.rm = TRUE) * sum(!is.na(contrast_df$sigma)) / 10
  }
  if (!is.null(variance_floor)) {
    bounded <- !is.na(squeezed_var$var.post) & squeezed_var$var.post < variance_floor
    squeezed_var$var.post[bounded] <- variance_floor
    squeezed_var$df.prior <- rep_len(squeezed_var$df.prior, length(bounded))
    squeezed_var$df.prior[bounded] <- Inf
  }
  squeezed_var <- tibble::as_tibble(squeezed_var)
  squeezed_var <- setNames(squeezed_var, paste0("moderated.", names(squeezed_var)))
  contrast_df <- dplyr::bind_cols(contrast_df, squeezed_var) |>
    dplyr::mutate(
      moderated.std.error = .data$std.error * sqrt(.data$moderated.var.post) / .data$sigma,
      moderated.statistic = .data$statistic * .data$sigma / sqrt(.data$moderated.var.post),
      moderated.df.total = .data[[df]] + .data$moderated.df.prior,
      moderated.p.value = 2 * pt(abs(.data$moderated.statistic), df = .data$moderated.df.total, lower.tail = FALSE)
    )

  conf_quantile <- -qt((1 - confint) / 2, df = contrast_df$moderated.df.total)
  contrast_df$moderated.conf.low <- contrast_df[[estimate]] - conf_quantile * contrast_df$moderated.std.error
  contrast_df$moderated.conf.high <- contrast_df[[estimate]] + conf_quantile * contrast_df$moderated.std.error
  dplyr::ungroup(contrast_df)
}

#' Finalize moderated contrast columns (prune raw stats, promote moderated, FDR).
#'
#' Shared post-moderation column surgery used by \code{ContrastsModerated} and
#' \code{ContrastsModeratedDEqMS}. When \code{all = FALSE} the raw (non-moderated)
#' statistics are dropped and the moderated ones renamed into the standard
#' column names; when \code{all = TRUE} both are kept and only the moderated
#' FDR column is added.
#'
#' @param result moderated contrast table (per-group output of \code{moderated_p_limma}).
#' @param p.adjust p-value adjustment function (signature compatible with
#'   \code{adjust_p_values}).
#' @param all keep all columns (default FALSE).
#' @return the finalized contrast table.
#' @keywords internal
#' @noRd
.finalize_moderated_columns <- function(result, p.adjust, all = FALSE) {
  if (all) {
    result <- p.adjust(
      result,
      column = "moderated.p.value",
      group_by_col = "contrast",
      newname = "FDR.moderated"
    )
  } else {
    result <- result |>
      dplyr::select(
        -c(
          "sigma",
          "df",
          "statistic",
          "std.error",
          "p.value",
          "conf.low",
          "conf.high",
          "FDR",
          "moderated.df.prior",
          "moderated.var.prior"
        )
      )
    result <- result |>
      dplyr::mutate(sigma = sqrt(moderated.var.post), .keep = "unused")
    result <- result |>
      dplyr::rename(
        conf.low = "moderated.conf.low",
        conf.high = "moderated.conf.high",
        statistic = "moderated.statistic",
        std.error = "moderated.std.error",
        df = "moderated.df.total",
        p.value = "moderated.p.value"
      )
    result <- p.adjust(
      result,
      column = "p.value",
      group_by_col = "contrast",
      newname = "FDR"
    )
  }
  result
}


#' Adjust p-values for multiple testing
#'
#' Benjamini-Hochberg false discovery rate, computed separately within each contrast (or other group).
#'
#' @return `contrast_df`, ungrouped, with the Benjamini-Hochberg adjusted p-values in the column `newname`,
#'   adjusted within each group of `group_by_col`.
#' @export
#' @param contrast_df data.frame with p-values to adjust
#' @param column name of column containing p-values
#' @param group_by_col column(s) to group by before adjusting (e.g. contrast), or NULL for no grouping
#' @param newname name of the new column with adjusted p-values
#' @examples
#'
#' bb <- c(runif(1000), rexp(1500,rate=5))
#' length(bb)
#' bb <- bb[bb < 1]
#' length(bb)
#' bb <- bb[1:2000]
#' hist(bb)
#' data <- data.frame(contrast = rep(LETTERS[1:5],400), p.value = bb)
#'
#' dataX <- adjust_p_values(data)
#' Adata <- dataX |> dplyr::filter(contrast == "A")
#' stopifnot(all.equal(Adata$FDR, p.adjust(Adata$p.value, method="BH")))
#' data2 <- adjust_p_values(data, group_by_col = NULL)
#' stopifnot(all.equal(data2$FDR, p.adjust(data2$p.value, method="BH")))
#'
#'
adjust_p_values <- function(
  contrast_df,
  column = "p.value",
  group_by_col = "contrast",
  newname = "FDR"
) {
  grouped_df <- contrast_df |>
    dplyr::group_by(across(all_of(group_by_col)))
  adjusted_df <- dplyr::mutate(grouped_df, !!newname := p.adjust(.data[[column]], method = "BH"))
  dplyr::ungroup(adjusted_df)
}


# ROPECA ----

#' Protein p-value from the p-value of the median fold-change peptide
#' @param median.p.value median of the peptide p-values of a protein
#' @param n.obs number of peptides
#' @param max.n limit number of peptides per protein.
#' @return a numeric vector of p-values, the beta-distribution probability that the median of
#'   `min(n.obs, max.n)` uniform p-values is at most `median.p.value`
#' @noRd
.get_p_values_pbeta <- function(median.p.value, n.obs, max.n = 10) {
  # beta distribution of the median of n.obs uniform p-values; both shapes equal n.obs / 2 + 0.5
  shape <- pmin(n.obs, max.n) / 2 + 0.5
  pbeta(median.p.value, shape1 = shape, shape2 = shape)
}


#' Compute protein-level fold changes and p-values using the beta distribution
#'
#' takes p-value of the scaled p-value
#'
#' @param contrasts_data data frame
#' @param contrast name of column with contrast identifier
#' @param subject_id name of column with typically protein Id
#' @param estimate name of column with effect size estimate
#' @param statistic statistic name of column with statistic (typically t-statistics)
#' @param p.value name of column with moderated.p.value
#' @param max.n used to limit the number of peptides in probability computation.
#' @return data.frame with columns
#'
#' @srrstats {G5.4a} ROPECA (ropeca_nested) is a derivation of PECA/ROPECA: per-peptide linear models, the median
#'   peptide fold change and a beta-distribution p-value of the median scaled p-value. test-Contrasts.R compares it
#'   with an independent per-peptide `stats::lm()`, median and `pbeta()` calculation and pins its conservative null.
#' @noRd
.summary_ropeca_median_p_scaled <- function(
  contrasts_data,
  contrast = "contrast",
  subject_id = "protein_Id",
  estimate = "diff",
  statistic = "statistic",
  p.value = "moderated.p.value",
  max.n = 10
) {
  contrasts_data <- contrasts_data |>
    dplyr::mutate(
      scaled.p = ifelse(.data[[estimate]] > 0, 1 - .data[[p.value]], .data[[p.value]] - 1)
    )

  summarized.protein <- contrasts_data |>
    group_by(across(all_of(c(subject_id, contrast)))) |>
    dplyr::summarize(
      n = dplyr::n(),
      n_not_na = sum(!is.na(.data$scaled.p)),
      mad.estimate = mad(.data[[estimate]], na.rm = TRUE),
      estimate = median(.data[[estimate]], na.rm = TRUE),
      statistic = median(.data[[statistic]], na.rm = TRUE),
      median.p.scaled = median(.data$scaled.p, na.rm = TRUE),
      avgAbd = median(.data$avgAbd, na.rm = TRUE)
    )

  summarized.protein <- summarized.protein |>
    dplyr::mutate(
      median.p.value = 1 - abs(.data$median.p.scaled),
      beta.based.significance = .get_p_values_pbeta(.data$median.p.value, .data$n_not_na, max.n = max.n),
      n.beta = pmin(.data$n_not_na, max.n),
      isSingular = FALSE
    )
  ungroup(summarized.protein)
}
