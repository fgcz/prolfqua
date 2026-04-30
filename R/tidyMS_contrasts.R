# Contrast linear function machinery ----

.model_coeff_matrix <- function(m) {
  data <- NULL
  if ("lmerModLmerTest" %in% class(m)) {
    data <- m@frame
    coeffs <- coefficients(summary(m))[, "Estimate"]
  } else {
    # for "lmerModLmerTest"
    data <- m$model
    coeffs <- coef(m)
  }
  interaction_columns <- intersect(attributes(terms(m))$term.labels, colnames(data))
  data <- make_interaction_column(data, interaction_columns, sep = ":")

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
  return(list(coeff_matrix = coeff_matrix, coeffs = coeffs))
}


.get_match_idx <- function(coeff_matrix, factor_level) {
  row_name_parts <- names_to_matrix(rownames(coeff_matrix), split = ":")
  factor_match <- apply(
    row_name_parts,
    2,
    function(x, factor_level) {
      x %in% factor_level
    },
    factor_level
  )
  idx <- which(apply(factor_match, 1, sum) > 0)
  return(idx)
}

.coeff_weights_factor_levels <- function(coeff_matrix) {
  get_coeffs <- function(factor_level, coeff_matrix) {
    idx <- .get_match_idx(coeff_matrix, factor_level)
    x <- as.list(apply(coeff_matrix[idx, , drop = FALSE], 2, mean))
    x <- tibble::as_tibble(x)
    tibble::add_column(x, "factor_level" = factor_level, .before = 1)
  }
  factor_levels <- unique(unlist(stringi::stri_split_fixed(rownames(coeff_matrix), ":")))
  weights_by_factor <- purrr::map_df(factor_levels, get_coeffs, coeff_matrix)
  return(weights_by_factor)
}

#' get linfct from model
#' @param m linear model
#' @family modelling
#' @export
#' @keywords internal
#' @examples
#'
#' m <- sim_make_model_lm("factors")
#' linfct <- linfct_from_model(m, as_list = TRUE)
#'
#' linfct$linfct_factors
#' linfct$linfct_interactions
#' lf <- matrix(
#' c(1, 1, 1, 1, 0.5, 0.5, 0, 1, 0, 1, 0.5, 0.5),
#' nrow = 4,
#' byrow = FALSE,
#' dimnames = list(c("BackgroundX", "BackgroundZ", "TreatmentA", "TreatmentB"),
#'                 c("(Intercept)", "TreatmentB", "BackgroundZ"))
#' )
#' stopifnot(lf == linfct$linfct_factors)
#' m <- sim_make_model_lm("interaction")
#' linfct <- linfct_from_model(m)
#'
#' m <- lm(Petal.Width ~ Species, data = iris)
#' linfct_from_model(m)
#' xx <- data.frame( Y = 1:10 , Condition = c(rep("a",5), rep("b",5)) )
#' m <- lm(Y ~ Condition, data = xx)
#' linfct_from_model(m)
#' xx <- data.frame( Y = 1:10 , Condition = c(rep("a",5), rep("b.b",5)) )
#' m <- lm(Y ~ Condition, data = xx)
#' linfct_from_model(m)
#' xx <- data.frame( Y = 1:10 , Condition = c(rep("a",5), rep("ab",5)) )
#' m <- lm(Y ~ Condition, data = xx)
#' linfct_from_model(m)
#'
linfct_from_model <- function(m, as_list = TRUE) {
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
    return(res)
  } else {
    do.call(rbind, res)
  }
}


#' linfct_matrix_contrasts
#'
#' When \code{options(prolfqua.vectorize = TRUE)} is set, dispatches to a
#' vectorized implementation that batch-evaluates all contrast expressions in a
#' single \code{dplyr::mutate()} call instead of looping per expression.
#' Set \code{options(prolfqua.vectorize = FALSE)} (the default) to use the
#' original per-expression loop.
#'
#' @export
#' @param linfct linear functions as created by linfct_from_model
#' @param contrasts named character vector of contrasts to determine linear functions for
#' @param p.message print messages default FALSE
#'
#' @family modelling
#' @keywords internal
#' @examples
#'
#' m <- sim_make_model_lm( "factors")
#' Contr <- c("TreatmentA_vs_B" = "TreatmentA - TreatmentB",
#'     "BackgroundX_vs_Z" = "BackgroundX - BackgroundZ",
#'     "IntoflintoA" = "`TreatmentA:BackgroundX` - `TreatmentA:BackgroundZ`",
#'     "IntoflintoB" = "`TreatmentB:BackgroundX` - `TreatmentB:BackgroundZ`",
#'     "IntoflintoX" = "`TreatmentA:BackgroundX` - `TreatmentB:BackgroundX`",
#'     "IntoflintoZ" = "`TreatmentA:BackgroundZ` - `TreatmentB:BackgroundZ`",
#'     "interactXZ" = "IntoflintoX - IntoflintoZ",
#'     "interactAB" = "IntoflintoA - IntoflintoB"
#'      )
#' linfct <- linfct_from_model(m, as_list = FALSE)
#' x <- linfct_matrix_contrasts(linfct, Contr )
#' stopifnot(sum(x["interactXZ",]) == 0 )
#' stopifnot(sum(x["interactAB",]) == 0 )
#'
#' m <- sim_make_model_lm( "interaction")
#' linfct <- linfct_from_model(m, as_list = FALSE)
#' x <- linfct_matrix_contrasts(linfct, Contr )
#' stopifnot(sum(x["interactXZ",]) ==1 )
#' stopifnot(sum(x["interactAB",]) ==1 )
#'
linfct_matrix_contrasts <- function(linfct, contrasts, p.message = FALSE) {
  if (isTRUE(getOption("prolfqua.vectorize"))) {
    return(linfct_matrix_contrasts_vectorized(linfct, contrasts, p.message = p.message))
  }
  linfct <- t(linfct)
  df <- tibble::as_tibble(linfct, rownames = "interaction")
  make_contrasts <- function(data, contrasts) {
    cnams <- base::setdiff(colnames(data), "interaction")
    failures <- list()
    for (i in seq_along(contrasts)) {
      if (p.message) {
        message(names(contrasts)[i], "=", contrasts[i], "\n")
      }
      contrast_name <- names(contrasts)[i]
      if (is.null(contrast_name) || !nzchar(contrast_name)) {
        contrast_name <- paste0("contrast_", i)
      }
      err <- tryCatch(
        {
          data <- dplyr::mutate(data, !!names(contrasts)[i] := !!rlang::parse_expr(contrasts[i]))
          NULL
        },
        error = function(e) {
          e
        }
      )
      if (inherits(err, "error")) {
        failures[[length(failures) + 1]] <- list(
          contrast = contrast_name,
          message = conditionMessage(err)
        )
      }
    }
    res <- data |> dplyr::select(-all_of(cnams))
    if (length(failures) > 0) {
      failure_df <- dplyr::bind_rows(failures)
      failure_names <- paste(failure_df$contrast, collapse = ", ")
      failure_messages <- unique(failure_df$message)
      failure_summary <- paste(utils::head(failure_messages, 3), collapse = "; ")
      msg <- sprintf(
        "linfct_matrix_contrasts: computed %d/%d contrasts; failed %d: %s. %s",
        ncol(res) - 1,
        length(contrasts),
        nrow(failure_df),
        failure_names,
        failure_summary
      )
      warning(msg, call. = FALSE)
    }
    return(res)
  }

  res <- make_contrasts(df, contrasts)
  res <- tibble::column_to_rownames(res, "interaction")
  res <- t(res)
  return(res)
}

#' Vectorized version of \code{\link{linfct_matrix_contrasts}}
#'
#' Same semantics but uses a single \code{dplyr::mutate(data, !!!parsed)} call
#' instead of one mutate per contrast. Falls back to per-expression evaluation
#' on error so that granular failure reporting is preserved.
#'
#' @inheritParams linfct_matrix_contrasts
#' @keywords internal
linfct_matrix_contrasts_vectorized <- function(linfct, contrasts, p.message = FALSE) {
  linfct <- t(linfct)
  df <- tibble::as_tibble(linfct, rownames = "interaction")
  make_contrasts <- function(data, contrasts) {
    cnams <- base::setdiff(colnames(data), "interaction")

    # Ensure all contrasts have names
    for (i in seq_along(contrasts)) {
      contrast_name <- names(contrasts)[i]
      if (is.null(contrast_name) || !nzchar(contrast_name)) {
        names(contrasts)[i] <- paste0("contrast_", i)
      }
    }

    # Pre-parse all contrast expressions
    parsed <- lapply(contrasts, rlang::parse_expr)
    names(parsed) <- names(contrasts)

    # Fast path: single mutate with !!! splicing
    err <- tryCatch(
      {
        data <- dplyr::mutate(data, !!!parsed)
        NULL
      },
      error = function(e) {
        e
      }
    )

    if (is.null(err)) {
      res <- data |> dplyr::select(-dplyr::all_of(cnams))
      return(res)
    }

    # Fallback: per-expression evaluation for granular error reporting
    failures <- list()
    for (i in seq_along(contrasts)) {
      if (p.message) {
        message(names(contrasts)[i], "=", contrasts[i], "\n")
      }
      err <- tryCatch(
        {
          data <- dplyr::mutate(data, !!names(contrasts)[i] := !!parsed[[i]])
          NULL
        },
        error = function(e) {
          e
        }
      )
      if (inherits(err, "error")) {
        failures[[length(failures) + 1]] <- list(
          contrast = names(contrasts)[i],
          message = conditionMessage(err)
        )
      }
    }
    res <- data |> dplyr::select(-dplyr::all_of(cnams))
    if (length(failures) > 0) {
      failure_df <- dplyr::bind_rows(failures)
      failure_names <- paste(failure_df$contrast, collapse = ", ")
      failure_messages <- unique(failure_df$message)
      failure_summary <- paste(utils::head(failure_messages, 3), collapse = "; ")
      msg <- sprintf(
        "linfct_matrix_contrasts: computed %d/%d contrasts; failed %d: %s. %s",
        ncol(res) - 1,
        length(contrasts),
        nrow(failure_df),
        failure_names,
        failure_summary
      )
      warning(msg, call. = FALSE)
    }
    return(res)
  }

  res <- make_contrasts(df, contrasts)
  res <- tibble::column_to_rownames(res, "interaction")
  res <- t(res)
  return(res)
}

#' create all possible contrasts
#' @export
#' @keywords internal
#' @family modelling
#' @examples
#' m <- sim_make_model_lm( "interaction")
#' linfct <- linfct_from_model(m)
#' xl <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_factors)
#' xx <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_interactions)
#' m <- sim_make_model_lm( "factor")
#' linfct <- linfct_from_model(m)
#' xl <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_factors)
#' xx <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_interactions)
#' m <- sim_make_model_lm( "parallel2")
#' linfct <- linfct_from_model(m)
#' xl <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_factors)
#' stopifnot(all(xl == c(0,-1)))
#' xx <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_interactions)
#' stopifnot(all(xx == c(0,-1)))
#'
#' m <- sim_make_model_lm( "parallel3")
#' linfct <- linfct_from_model(m)
#' xl <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_factors)
#' stopifnot(all(xl == c(0,0,0,-1,0,1,0,-1,-1)))
#' xx <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_interactions)
#' stopifnot(all(xl == c(0,0,0,-1,0,1,0,-1,-1)))
linfct_all_possible_contrasts <- function(lin_int) {
  combs <- combn(nrow(lin_int), 2)
  names <- rownames(lin_int)
  newnames <- rep("", ncol(combs))
  new_lin_fct <- matrix(NA, nrow = ncol(combs), ncol = ncol(lin_int))
  for (i in seq_len(ncol(combs))) {
    newnames[i] <- paste(names[combs[, i]], collapse = " - ")
    new_lin_fct[i, ] <- lin_int[combs[1, i], ] - lin_int[combs[2, i], ]
  }
  rownames(new_lin_fct) <- newnames
  colnames(new_lin_fct) <- colnames(lin_int)
  return(new_lin_fct)
}
#' create contrasts between factor levels
#'
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#'
#' m <- sim_make_model_lm( "interaction")
#' xl <- linfct_factors_contrasts(m)
#' m <- lm(Petal.Width ~ Species, data = iris)
#' linfct_factors_contrasts(m)
linfct_factors_contrasts <- function(m) {
  ffac <- attributes(terms(m))$term.labels
  ffac <- ffac[!grepl(":", ffac)] # remove interactions
  linfct_factors <- linfct_from_model(m)$linfct_factors

  factor_depths <- rownames(linfct_factors)
  res <- vector(length(ffac), mode = "list")
  for (i in seq_along(ffac)) {
    fac <- ffac[i]
    idx <- grep(fac, factor_depths)
    linfct_m <- linfct_factors[idx, ]
    res[[i]] <- linfct_all_possible_contrasts(linfct_m)
  }
  res <- do.call(rbind, res)
  return(res)
}

# Computing contrasts helpers -----

# compute contrasts for full models
# @param m linear model generated using lm
# @param linfct linear function
# @param strategy optional strategy for df and sigma computation
# @param coef coefficients vector, default from model
# @param Sigma.hat variance-covariance matrix, default from model
# @param confint which confidence interval to determine
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
    std.error <- sqrt(diag(linfct %*% Sigma.hat %*% t(linfct)))
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
  return(res)
}

#' handles incomplete models
#'
#' only keeps non NA coefficients.
#'
#' When \code{options(prolfqua.vectorize = TRUE)} is set, dispatches to a
#' vectorized implementation that uses matrix multiplication instead of a
#' per-row loop. Set \code{options(prolfqua.vectorize = FALSE)} (the default)
#' to use the original loop.
#'
#' @param m linear model generated using lm
#' @param linfct linear function
#' @param confint confidence interval default 0.95
#'
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#' m <- sim_make_model_lm( "factors")
#' linfct <- linfct_from_model(m)$linfct_factors
#' compute_contrast(m, linfct, confint = 0.95)
#' compute_contrast(m, linfct, confint = 0.99)
#'
compute_contrast <- function(m, linfct, confint = 0.95) {
  if (isTRUE(getOption("prolfqua.vectorize"))) {
    return(compute_contrast_vectorized(m, linfct, confint = confint))
  }
  Sigma.hat <- vcov(m)

  coef <- na.omit(coefficients(m))

  res <- vector(nrow(linfct), mode = "list")
  for (i in seq_len(nrow(linfct))) {
    linfct_v <- linfct[i, , drop = FALSE]
    idx <- which(linfct_v != 0)
    nam <- colnames(linfct_v)[idx]

    if (all(nam %in% names(coef))) {
      linfct_v_red <- linfct_v[, nam, drop = FALSE]
      Sigma.hat_red <- Sigma.hat[nam, nam, drop = FALSE]
      coef_red <- coef[nam]
      stopifnot(all.equal(colnames(linfct_v_red), colnames(Sigma.hat_red)))
      stopifnot(all.equal(colnames(linfct_v_red), names(coef_red)))
      res[[i]] <- .compute_contrast(m, linfct_v_red, coef = coef_red, Sigma.hat = Sigma.hat_red, confint = confint)
    } else {
      res[[i]] <- data.frame(
        lhs = rownames(linfct_v),
        sigma = sigma(m),
        df = df.residual(m),
        estimate = NA,
        std.error = NA,
        statistic = NA,
        p.value = NA,
        conf.low = NA,
        conf.high = NA,
        stringsAsFactors = FALSE
      )
    }
  }
  res <- dplyr::bind_rows(res)
  return(res)
}

#' Vectorized version of \code{\link{compute_contrast}}
#'
#' Same semantics but uses vectorized matrix multiplication instead of a
#' per-row loop. NAs in \code{coefficients(m)} propagate naturally via
#' \code{linfct \%*\% coef} for rows that reference missing coefficients.
#'
#' @inheritParams compute_contrast
#' @keywords internal
compute_contrast_vectorized <- function(m, linfct, confint = 0.95) {
  n <- nrow(linfct)
  if (n == 0) {
    return(data.frame(
      lhs = character(0),
      sigma = numeric(0),
      df = numeric(0),
      estimate = numeric(0),
      std.error = numeric(0),
      statistic = numeric(0),
      p.value = numeric(0),
      conf.low = numeric(0),
      conf.high = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  coef_full <- coefficients(m)
  available <- names(na.omit(coef_full))
  Sigma.hat <- vcov(m)
  df <- df.residual(m)
  sig <- sigma(m)

  # Align coefficients with linfct columns by name
  coef_aligned <- coef_full[colnames(linfct)]
  na_coefs <- is.na(coef_aligned)

  # Use zero-filled coefficients for multiplication (0 * NA = NA in R, but 0 * 0 = 0)
  coef_zero <- coef_aligned
  coef_zero[na_coefs] <- 0
  estimate <- as.vector(linfct %*% coef_zero)

  # Mark rows invalid if any non-zero weight touches an NA coefficient
  if (any(na_coefs)) {
    invalid <- as.logical(linfct[, na_coefs, drop = FALSE] %*% rep(1, sum(na_coefs)) != 0)
    estimate[invalid] <- NA
  }

  if (df > 0) {
    # std.error: vcov only covers available (non-NA) coefficients
    available_cols <- intersect(colnames(linfct), available)
    linfct_avail <- linfct[, available_cols, drop = FALSE]
    Sigma.hat_avail <- Sigma.hat[available_cols, available_cols, drop = FALSE]
    std.error <- sqrt(diag(linfct_avail %*% Sigma.hat_avail %*% t(linfct_avail)))
    std.error[is.na(estimate)] <- NA
    statistic <- estimate / std.error
    p.value <- pt(abs(statistic), df = df, lower.tail = FALSE) * 2
    prqt <- -qt((1 - confint) / 2, df = df)
    conf.low <- estimate - prqt * std.error
    conf.high <- estimate + prqt * std.error
  } else {
    std.error <- rep(NA_real_, n)
    statistic <- rep(NA_real_, n)
    p.value <- rep(NA_real_, n)
    conf.low <- rep(NA_real_, n)
    conf.high <- rep(NA_real_, n)
  }

  res <- data.frame(
    lhs = rownames(linfct),
    sigma = sig,
    df = df,
    estimate = estimate,
    std.error = std.error,
    statistic = statistic,
    p.value = p.value,
    conf.low = conf.low,
    conf.high = conf.high,
    stringsAsFactors = FALSE
  )
  return(res)
}

#' applies contrast computation using lmerTest::contest function
#' @param model mixed effects model
#' @param linfct linear function
#' @param ddf method to determine denominator degrees of freedom
#' @family modelling
#' @export
#' @keywords internal
#' @examples
#'
#' mb <- sim_make_model_lmer("interaction")
#' summary(mb)
#'
#' linfct <- linfct_from_model(mb)
#' names(linfct)
#' compute_lmer_contrast(mb, linfct$linfct_factors)
#' compute_lmer_contrast(mb, linfct$linfct_interactions)
#' length(lme4::fixef(mb))
#' lmerTest::contest(mb, c( 0 ,1 , 0 , 0),joint = FALSE)
#' summary(mb)
#'
compute_lmer_contrast <- function(model, linfct, ddf = c("Satterthwaite", "Kenward-Roger")) {
  ddf <- match.arg(ddf)
  if (length(lme4::fixef(model)) != ncol(linfct)) {
    warning("Model is rank deficient!")
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
  return(res)
}

# computing contrast ----
#' pivot model contrasts matrix to wide format produced by `contrasts_linfct` and ...
#' @return Contrast definitions or contrast results.
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
  return(res)
}
#' compute group averages
#'
#' used in p2621, p2109
#'
#' @export
#' @keywords internal
#' @examples
#' modelSummary_A <- sim_build_models_lm()
#' m <- get_complete_model_fit(modelSummary_A$model_df)
#'
#' factor_contrasts <- linfct_factors_contrasts( m$linear_model[[1]])
#'
#' factor_levelContrasts <- contrasts_linfct( m,
#'         factor_contrasts,
#'         subject_id = "protein_Id",
#'         contrastfun = prolfqua::compute_contrast)
contrasts_linfct <- function(models, linfct, subject_id = "protein_Id", contrastfun = prolfqua::compute_lmer_contrast) {
  message("contrasts_linfct")
  modelcol <- "linear_model"

  # Normalize: if linfct is a matrix, replicate it for each model
  if (inherits(linfct, "matrix")) {
    linfct <- rep(list(linfct), nrow(models))
  }
  if (!is.list(linfct) || length(linfct) != nrow(models)) {
    stop("linfct must be either a matrix or a list of length == nrow(models)")
  }

  interaction_models <- vector(mode = "list", length = nrow(models))
  pb <- progress::progress_bar$new(total = nrow(models))
  for (i in seq_along(models[[modelcol]])) {
    interaction_models[[i]] <- contrastfun(models[[modelcol]][[i]], linfct = linfct[[i]])
    pb$tick()
  }
  interaction_model_matrix <- models
  interaction_model_matrix$contrast <- interaction_models

  mclass <- function(x) {
    class(x)[1]
  }

  interaction_model_matrix <- interaction_model_matrix |>
    dplyr::mutate(classC = purrr::map_chr(.data$contrast, mclass))

  failed_mask <- interaction_model_matrix$classC == "logical"
  n_failed <- sum(failed_mask)
  if (n_failed > 0) {
    failed_ids <- interaction_model_matrix[[subject_id[1]]][failed_mask]
    warning(
      "contrasts_linfct: dropped ",
      n_failed,
      " of ",
      nrow(interaction_model_matrix),
      " proteins with failed contrasts: ",
      paste(failed_ids, collapse = ", "),
      call. = FALSE
    )
  }

  interaction_model_matrix <- interaction_model_matrix |>
    dplyr::filter(.data$classC != "logical")

  contrasts <- interaction_model_matrix |>
    dplyr::select(all_of(c(subject_id, "contrast"))) |>
    tidyr::unnest(cols = c("contrast"))

  # take sigma and df from somewhere else.
  model_infos <- models |>
    dplyr::select(all_of(c(subject_id, "isSingular", "sigma.model" = "sigma", "df.residual.model" = "df.residual"))) |>

    dplyr::distinct()
  contrasts <- dplyr::inner_join(contrasts, model_infos, by = subject_id)
  return(ungroup(contrasts))
}
