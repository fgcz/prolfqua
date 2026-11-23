# Contrast linear function machinery ----

#' Split two-sided contrast specifications into left/right group columns.
#'
#' Shared implementation behind the \code{get_contrast_sides()} methods of the
#' contrast classes. Keeps only contrasts containing a \code{-}, strips
#' backticks/spaces, and separates the right-hand side on \code{-}.
#'
#' @param contrasts named character vector of contrast specifications.
#' @return a tibble with columns \code{contrast}, \code{group_1}, \code{group_2}.
#' @note Preserves the historical parser behaviour: it does not handle \code{+}
#'   contrasts or group names that themselves contain \code{-}.
#' @keywords internal
#' @noRd
parse_contrast_sides <- function(contrasts) {
  tt <- contrasts[grep("-", contrasts, fixed = TRUE)]
  tt <- tibble::tibble(contrast = names(tt), rhs = tt)
  tt |>
    dplyr::mutate(rhs = gsub("[` ]", "", rhs)) |>
    tidyr::separate_wider_delim(
      "rhs",
      delim = "-",
      names = c("group_1", "group_2"),
      too_few = "align_start",
      too_many = "drop"
    )
}

#' Represent observed factor-level combinations as coefficient weights
#' @param m fitted regression model
#' @noRd
.model_coeff_matrix <- function(m) {
  if (inherits(m, "lmerModLmerTest")) {
    data <- m@frame
    coeffs <- coefficients(summary(m))[, "Estimate"]
  } else {
    # lm-like fits (lm, rlm, rfit, logistf, ...)
    data <- m$model
    coeffs <- coef(m)
  }
  interaction_columns <- intersect(attributes(terms(m))$term.labels, colnames(data))
  data <- .make_interaction_column(data, interaction_columns, sep = ":")

  inter <- unique(data$interaction)
  coeff_matrix <- matrix(0, nrow = length(inter), ncol = length(coeffs))
  rownames(coeff_matrix) <- inter
  colnames(coeff_matrix) <- names(coeffs)
  coeff_matrix[, 1] <- 1
  non_intercept_coeffs <- coeffs[-1]
  for (i in seq_along(non_intercept_coeffs)) {
    # the grep is needed to extract coefficients of interaction terms belonging to a factor
    # I am using wor boundaries "\\b" to allow for factor levels that are substrings.
    position_index <- grep(paste0("\\b", names(non_intercept_coeffs)[i], "\\b"), inter)
    coeff_matrix[position_index, i + 1] <- 1
  }
  list(coeff_matrix = coeff_matrix, coeffs = coeffs)
}


#' Average combination weights for each individual factor level
#' @param coeff_matrix factor-combination coefficient weights
#' @noRd
.coeff_weights_factor_levels <- function(coeff_matrix) {
  row_name_parts <- stringi::stri_split_fixed(rownames(coeff_matrix), ":")
  get_coeffs <- function(factor_level) {
    idx <- which(vapply(row_name_parts, function(x) factor_level %in% x, logical(1)))
    x <- as.list(colMeans(coeff_matrix[idx, , drop = FALSE]))
    x <- tibble::as_tibble(x)
    tibble::add_column(x, "factor_level" = factor_level, .before = 1)
  }
  purrr::list_rbind(purrr::map(unique(unlist(row_name_parts)), get_coeffs))
}

#' Get group-level linear functions from a fitted model
#'
#' @details The fitted model supplies its coefficient names and observed
#'   categorical factor combinations. Each combination is expressed as weights
#'   on the fitted coefficients; weights for an individual factor level are
#'   averaged equally over the observed combinations containing that level.
#'   Thus \code{linfct_factors} describes marginal factor-level means and
#'   \code{linfct_interactions} describes cell means. With treatment coding,
#'   a reference-level mean includes the intercept; a non-reference mean also
#'   includes its treatment coefficient. The \emph{Under the hood} article
#'   shows a complete \code{~ A * B} example. These helpers assume categorical
#'   predictors and the usual intercept-based coding; they are not a general
#'   marginal-means implementation for continuous predictors or arbitrary coding.
#' @param as_list if TRUE, return the factor and interaction matrices as a list; if FALSE, as one matrix
#' @param m linear model
#' @return if `as_list = TRUE`, a list with the matrices `linfct_factors` (one row per factor level) and
#'   `linfct_interactions` (one row per observed factor combination), each with one column per model coefficient;
#'   otherwise the two matrices bound by rows into a single matrix
#' @srrstats {RE1.1} The article and linfct details show treatment-coded designs, cell means and tested interaction
#'   weights.
#' @noRd
.linfct_from_model <- function(m, as_list = TRUE) {
  coeff_result <- .model_coeff_matrix(m)
  sorted_coeff_matrix <- coeff_result$coeff_matrix[order(rownames(coeff_result$coeff_matrix)), ]

  l_factors <- .coeff_weights_factor_levels(sorted_coeff_matrix)
  linfct_factors <- l_factors |>
    dplyr::select(-dplyr::all_of("factor_level")) |>
    data.matrix()

  rownames(linfct_factors) <- l_factors$factor_level
  linfct_factors <- linfct_factors[order(rownames(linfct_factors)), ]
  res <- list(linfct_factors = linfct_factors, linfct_interactions = sorted_coeff_matrix)

  if (as_list) {
    res
  } else {
    do.call(rbind, res)
  }
}


#' Convert named group-mean contrasts to coefficient weights
#'
#' @details Rows of \code{linfct} represent group or interaction means and
#'   columns represent fitted coefficients. Contrast expressions combine these
#'   rows, so \code{"group_A - group_Ctrl"} becomes their row difference.
#'   Quote interaction names containing colons with backticks. Expressions may
#'   reference earlier named contrasts. A row naming a group unavailable in a
#'   particular model is dropped with a warning. The returned matrix has one
#'   row per estimable contrast and coefficient columns in model order.
#'
#' @param linfct linear functions as created by linfct_from_model
#' @param contrasts named character vector of contrasts to determine linear functions for
#' @param p.message print messages default FALSE
#' @return a matrix with one row per computable contrast and one column per model coefficient
#'
#' @noRd
.linfct_matrix_contrasts <- function(linfct, contrasts, p.message = FALSE) {
  # one coefficient vector per group; computed contrasts join them so later contrasts can use earlier ones
  terms <- stats::setNames(lapply(seq_len(nrow(linfct)), function(i) linfct[i, ]), rownames(linfct))
  computed <- character(0)
  failures <- list()
  for (i in seq_along(contrasts)) {
    contrast_name <- names(contrasts)[i]
    if (is.null(contrast_name) || !nzchar(contrast_name)) {
      contrast_name <- paste0("contrast_", i)
    }
    if (p.message) {
      message(contrast_name, "=", contrasts[i], "\n")
    }
    value <- .eval_contrast(contrasts[[i]], terms)
    if (is.character(value)) {
      failures[[length(failures) + 1]] <- list(contrast = contrast_name, message = value)
    } else {
      terms[[contrast_name]] <- value
      computed <- c(computed, contrast_name)
    }
  }
  if (length(failures) > 0) {
    failure_df <- dplyr::bind_rows(failures)
    msg <- sprintf(
      "linfct_matrix_contrasts: computed %d/%d contrasts; failed %d: %s. %s",
      length(computed),
      length(contrasts),
      nrow(failure_df),
      toString(failure_df$contrast),
      paste(utils::head(unique(failure_df$message), 3), collapse = "; ")
    )
    warning(msg, call. = FALSE)
  }
  matrix(
    as.numeric(unlist(terms[computed], use.names = FALSE)),
    nrow = length(computed),
    ncol = ncol(linfct),
    byrow = TRUE,
    dimnames = list(computed, colnames(linfct))
  )
}

#' Evaluate one contrast expression on the group coefficient vectors
#'
#' A contrast is arithmetic (`+ - * /`, parentheses, numbers) on group and earlier contrast names and names at
#' least one of them, so its value has one weight per coefficient. Any other function is refused before evaluation.
#' @param contrast contrast expression, e.g. "group_A - group_B"
#' @param terms named list of coefficient vectors (groups and earlier contrasts), all of the same length
#' @return the coefficient vector of the contrast, or the reason it could not be computed (character)
#' @noRd
.eval_contrast <- function(contrast, terms) {
  expr <- tryCatch(rlang::parse_expr(contrast), error = function(e) e)
  if (inherits(expr, "error")) {
    return(conditionMessage(expr))
  }
  functions <- setdiff(all.names(expr), all.vars(expr))
  refused <- setdiff(functions, c("+", "-", "*", "/", "("))
  if (length(refused) > 0) {
    abort_bad_argument(
      "contrasts",
      "use only + - * / and parentheses on group names",
      not = sprintf("`%s`, which calls %s", contrast, toString(paste0("`", refused, "`")))
    )
  }
  if (length(all.vars(expr)) == 0) {
    abort_bad_argument("contrasts", "name at least one group", not = sprintf("`%s`", contrast))
  }
  # a group missing from this model (partial design) is reported without raising an error
  unknown <- setdiff(all.vars(expr), names(terms))
  if (length(unknown) > 0) {
    return(paste0("not in the model: ", toString(unknown)))
  }
  operators <- list2env(list(`+` = `+`, `-` = `-`, `*` = `*`, `/` = `/`, `(` = `(`), parent = emptyenv())
  value <- tryCatch(eval(expr, terms, operators), error = function(e) e)
  if (inherits(value, "error")) {
    return(conditionMessage(value))
  }
  value
}

# Computing contrasts helpers -----

#' Variance of each contrast row, l' V l over the coefficients the row uses
#'
#' Coefficients with zero weight are left out, so an infinite or NA variance of an unused coefficient (rfit on
#' small groups) does not turn the row into NaN through 0 * Inf.
#' @param linfct contrast matrix, one row per contrast
#' @param Sigma.hat variance-covariance matrix of the coefficients in the columns of `linfct`
#' @noRd
.contrast_variance <- function(linfct, Sigma.hat) {
  vapply(
    seq_len(nrow(linfct)),
    function(i) {
      used <- !.is_zero(linfct[i, ])
      l <- linfct[i, used]
      sum(l * (Sigma.hat[used, used, drop = FALSE] %*% l))
    },
    numeric(1)
  )
}

#' compute contrasts for full models
#'
#' @param m linear model generated using lm
#' @param linfct linear function
#' @param strategy optional strategy for df and sigma computation
#' @param coef coefficients vector, default from model
#' @param Sigma.hat variance-covariance matrix, default from model
#' @param confint which confidence interval to determine
#' @noRd
.compute_contrast <- function(m, linfct, strategy = NULL, coef = coefficients(m), Sigma.hat = vcov(m), confint = 0.95) {
  if (is.null(strategy)) {
    df <- df.residual(m)
    sigma <- sigma(m)
  } else {
    df <- strategy$df_residual(m)
    sigma <- strategy$sigma(m)
  }

  estimate <- linfct %*% t(t(coef))

  if (df > 0) {
    std.error <- sqrt(.contrast_variance(linfct, Sigma.hat))
    statistic <- estimate / std.error

    p.value <- pt(abs(statistic), df = df, lower.tail = FALSE) * 2
    prqt <- -qt((1 - confint) / 2, df = df)
    conf.low <- estimate - prqt * std.error
    conf.high <- estimate + prqt * std.error
  } else {
    std.error <- NA
    statistic <- NA
    p.value <- NA
    conf.low <- NA
    conf.high <- NA
  }

  res <- data.frame(
    lhs = rownames(linfct),
    sigma = sigma,
    df = df,
    estimate = estimate,
    std.error = std.error,
    statistic = statistic,
    p.value = p.value,
    conf.low = conf.low,
    conf.high = conf.high,
    stringsAsFactors = FALSE
  )
  res
}

#' Compute contrasts for possibly incomplete models
#'
#' only keeps non NA coefficients.
#'
#' @param m linear model generated using lm
#' @param linfct linear function
#' @param confint confidence interval default 0.95
#' @return a data.frame with one row per row of `linfct` and columns `lhs`, `sigma`, `df`, `estimate`,
#'   `std.error`, `statistic`, `p.value`, `conf.low` and `conf.high`; the estimates are `NA` for rows that use a
#'   coefficient the model could not estimate
#'
#' @srrstats {RE7.1} Tests recover exact group differences and exact paired-design means.
#' @noRd
.compute_estimable_contrast <- function(m, linfct, confint = 0.95) {
  if (nrow(linfct) == 0L) {
    return(dplyr::bind_rows(list()))
  }
  zero_rows <- rowSums(!.is_zero(linfct)) == 0
  if (any(zero_rows)) {
    abort_bad_argument(
      "linfct",
      "have a non-zero coefficient in every row",
      not = paste0("all-zero rows ", toString(rownames(linfct)[zero_rows]))
    )
  }
  coef <- stats::na.omit(stats::coefficients(m))
  estimable_coef <- intersect(colnames(linfct), names(coef))
  # a row is estimable when it puts no weight on a coefficient the model could not estimate
  other <- linfct[, setdiff(colnames(linfct), estimable_coef), drop = FALSE]
  estimable <- rowSums(!.is_zero(other)) == 0
  res <- .compute_contrast(
    m,
    linfct[, estimable_coef, drop = FALSE],
    coef = coef[estimable_coef],
    Sigma.hat = stats::vcov(m)[estimable_coef, estimable_coef, drop = FALSE],
    confint = confint
  )
  if (!all(estimable)) {
    res[!estimable, c("estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")] <- NA
  }
  res
}

#' Compute contrasts of a mixed-effects model with lmerTest
#' @param model mixed effects model
#' @param linfct linear function
#' @param ddf method to determine denominator degrees of freedom
#' @return a tibble with one row per row of `linfct` and columns `lhs`, `estimate`, `std.error`, `df`,
#'   `statistic`, `conf.low`, `conf.high`, `p.value` and `sigma`; `NA` (with a warning) if `linfct` has no rows or
#'   the model is rank-deficient
#' @noRd
.compute_lmer_contrast <- function(model, linfct, ddf = c("Satterthwaite", "Kenward-Roger")) {
  ddf <- match.arg(ddf)
  if (nrow(linfct) == 0L) {
    # `linfct_matrix_contrasts()` drops a contrast whose coefficients this
    # model does not carry, and for some proteins that leaves nothing at all.
    # `lmerTest::contest()` indexes `1:nrow(L)`, so an empty matrix reaches
    # `L[1, ]` and aborts the whole run with "subscript out of bounds".
    # Nothing is estimable here, so report it like a rank-deficient model.
    warning("No estimable contrast for this model.")
    return(NA)
  }
  if (length(lme4::fixef(model)) != ncol(linfct)) {
    warning("The model is rank-deficient.")
    return(NA) # catch rank defficient
  } else {
    res <- lmerTest::contest(model, linfct, joint = FALSE, confint = TRUE, ddf = ddf)
  }
  res <- tibble::as_tibble(res, rownames = "lhs")
  res$sigma <- sigma(model)
  res <- res |>
    dplyr::rename(
      estimate = "Estimate",
      std.error = "Std. Error",
      statistic = "t value",
      p.value = "Pr(>|t|)",
      conf.low = "lower",
      conf.high = "upper"
    )
  res
}

# computing contrast ----
#' Pivot long-format model contrasts to wide format
#'
#' Turns a contrast table with one row per subject and contrast into one row per subject, as used by the `to_wide()`
#' method of the contrast classes.
#' @return A tibble with one row per subject and one column per value column and contrast, named
#'   `<column>.<contrast>`.
#' @export
#' @family modelling
#' @param model_interaction_contrasts data.frame with contrast results in long format
#' @param subject_id column name(s) identifying subjects (e.g. protein_Id)
#' @param columns character vector of value columns to pivot wide
#' @param contrast column name containing contrast labels
#' @examples
#'
#' # this function is used by the contrast classes to implement the to wide method
#' contrast_df <- data.frame(
#'   protein_Id = c("P1", "P1", "P2", "P2"),
#'   lhs = c("A_vs_B", "C_vs_D", "A_vs_B", "C_vs_D"),
#'   estimate = c(1.0, 0.5, -0.2, 0.7),
#'   p.value = c(0.01, 0.2, 0.4, 0.05),
#'   p.value.adjusted = c(0.02, 0.25, 0.4, 0.08)
#' )
#' pivot_model_contrasts_to_wide(contrast_df)
#'
pivot_model_contrasts_to_wide <- function(
  model_interaction_contrasts,
  subject_id = "protein_Id",
  columns = c("estimate", "p.value", "p.value.adjusted"),
  contrast = "lhs"
) {
  res <- model_interaction_contrasts |>
    dplyr::select(dplyr::all_of(c(subject_id, contrast, columns))) |>
    tidyr::pivot_wider(
      names_from = dplyr::all_of(contrast),
      values_from = dplyr::all_of(columns),
      names_glue = paste0("{.value}.{", contrast, "}")
    )
  res
}
#' Compute group averages or contrasts for each model
#'
#' used in p2621, p2109
#'
#' @param models model table with one fit per subject in the column `linear_model`, e.g. the `model_df` of
#'   [model_analyse()]
#' @param linfct linear functions: a matrix with one row per contrast and one column per model coefficient, or
#'   a list with one such matrix per model
#' @param subject_id columns identifying a subject
#' @param contrastfun function computing the contrasts of one fit, e.g. [compute_lmer_contrast()]
#' @return a tibble with one row per subject and row of `linfct`, holding the `subject_id` columns, the columns
#'   returned by `contrastfun`, and `isSingular`, `sigma.model` and `df.residual.model`; subjects whose contrasts
#'   failed are dropped with a warning
#' @noRd
.contrasts_linfct <- function(models, linfct, subject_id = "protein_Id", contrastfun = .compute_lmer_contrast) {
  modelcol <- "linear_model"

  # Normalize: if linfct is a matrix, replicate it for each model
  if (inherits(linfct, "matrix")) {
    linfct <- rep(list(linfct), nrow(models))
  }
  if (!is.list(linfct) || length(linfct) != nrow(models)) {
    abort_prolfqua("`linfct` must be a matrix or a list with one element per model.")
  }

  pb <- .make_progress(nrow(models), label = "contrasts")
  computed <- purrr::map2(models[[modelcol]], linfct, function(model, l) {
    res <- .fit_capturing_warnings(contrastfun, model, linfct = l)
    pb$tick()
    res
  })
  interaction_models <- purrr::map(computed, "fit")
  contrast_warnings <- as.character(unlist(purrr::map(computed, "warnings")))
  interaction_model_matrix <- models
  interaction_model_matrix$contrast <- interaction_models

  failed <- purrr::map_lgl(interaction_models, is.logical)
  n_failed <- sum(failed)
  if (n_failed > 0) {
    failed_ids <- interaction_model_matrix[[subject_id[1]]][failed]
    reasons <- table(contrast_warnings)
    rlang::warn(
      paste0(
        "contrasts_linfct: dropped ",
        n_failed,
        " of ",
        nrow(interaction_model_matrix),
        " proteins with failed contrasts: ",
        toString(failed_ids),
        if (length(reasons) > 0L) paste0(". Reasons: ", toString(paste0(names(reasons), " (", reasons, ")")))
      ),
      class = "prolfqua_fit_failed"
    )
  }

  contrasts <- interaction_model_matrix[!failed, ] |>
    dplyr::select(all_of(c(subject_id, "contrast"))) |>
    tidyr::unnest(cols = "contrast")

  # take sigma and df from somewhere else.
  model_infos <- models |>
    dplyr::select(all_of(c(subject_id, "isSingular", "sigma.model" = "sigma", "df.residual.model" = "df.residual"))) |>

    dplyr::distinct()
  contrasts <- dplyr::inner_join(contrasts, model_infos, by = subject_id, relationship = "many-to-one")
  ungroup(contrasts)
}
