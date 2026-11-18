# ContrastsROPECA -----

#' Map ROPECA output to the standard contrast columns: beta-based significance
#'
#' becomes p.value and FDR; std.error = diff / statistic, sigma = mad.estimate,
#' df = n_not_na, and a t-based 95% interval. `all = TRUE` appends the remaining
#' ROPECA columns after the standard ones.
#' @param res contrast result table
#' @param subject_id column names identifying the modelled subjects
#' @param all retain all diagnostic columns
#' @noRd
.ropeca_standard_columns <- function(res, subject_id, all = FALSE) {
  res <- dplyr::rename(res, p.value = "beta.based.significance", FDR = "FDR.beta.based.significance")
  res$std.error <- ifelse(res$statistic != 0, res$diff / res$statistic, NA_real_)
  res$sigma <- res$mad.estimate
  res$df <- res$n_not_na
  prqt <- stats::qt(0.975, df = pmax(res$df, 1))
  res$conf.low <- res$diff - prqt * abs(res$std.error)
  res$conf.high <- res$diff + prqt * abs(res$std.error)
  standard_cols <- c(
    subject_id,
    "modelName",
    "estimate_type",
    "contrast",
    "avgAbd",
    "diff",
    "FDR",
    "statistic",
    "std.error",
    "df",
    "p.value",
    "conf.low",
    "conf.high",
    "sigma"
  )
  keep <- if (all) c(standard_cols, setdiff(colnames(res), standard_cols)) else standard_cols
  res[, keep, drop = FALSE]
}

#'
#' ROPECA reproducibility-optimization method
#'
#' ROPECA optimizes the reproducibility of statistical testing
#'  on peptide-level and aggregates the peptide-level changes
#'   to determine differential protein-level expression.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#'
#' istar <- prolfqua::sim_lfq_data_peptide_config(Nprot=20)
#' istar_data <- istar$data
#' modelFunction <-
#'   strategy_lm("abundance  ~ group_")
#' pepIntensity <- istar_data
#' config <- istar$config$clone(deep = TRUE)
#' config$hierarchy_depth <- 2
#' config$hierarchy_keys_depth()
#'
#' mod <- build_model(
#'  pepIntensity,
#'  modelFunction,
#'  subject_id = config$hierarchy_keys_depth())
#'
#'  Contr <- c("AvsCtrl" = "group_A - group_Ctrl")
#'
#'
#'  contr <- prolfqua::Contrasts$new(mod, Contr)
#'  dim(contr$get_contrasts())
#'  contrM <- prolfqua::ContrastsModerated$new(contr)
#'  dim(contrM$get_contrasts())
#'  contrast <- prolfqua::ContrastsROPECA$new(contrM)
#'  contrast$get_contrasts()
#'  contrast <- prolfqua::ContrastsROPECA$new(contr)
#'  tmp <- contrast$get_contrasts()
#'  dim(tmp)
#'  pl <- contrast$get_Plotter()
#'  contrast$to_wide()
#'  contrast$get_linfct()
#'  contrast$get_contrast_sides()
#'  pl$histogram()
#'  pl$ma_plot()
#'
ContrastsROPECA <- R6::R6Class(
  "ContrastsROPECA",
  inherit = ContrastsInterface,
  public = list(
    #' @field Contrast Contrast
    Contrast = NULL,
    #' @field contrast_result contrast result
    contrast_result = NULL,
    #' @field model_name model name
    model_name = character(),
    #' @field subject_id columns with protein ID's
    subject_id = character(),
    #' @field p.adjust method to use for p.value adjustment
    p.adjust = NULL,
    #' @description
    #' initialize
    #' @param Contrast e.g. instance of Contrasts class, or ContrastsModerated
    #' @param model_name default ROPECA
    #' @param p.adjust function to use for p.value adjustment
    initialize = function(Contrast, model_name = "ROPECA", p.adjust = prolfqua::adjust_p_values) {
      self$Contrast <- Contrast
      if (length(Contrast$subject_id) <= 1) {
        abort_bad_argument("Contrast", "hold peptide-level contrasts, with more than one `subject_id` column")
      }
      self$model_name <- model_name
      self$subject_id <- Contrast$subject_id
      self$p.adjust <- p.adjust
    },
    #' @description
    #' show names of contrasts
    #' @return data.frame
    get_contrast_sides = function() {
      self$Contrast$get_contrast_sides()
    },
    #' @description
    #' get linear function used to determine contrasts
    #' @return data.frame
    get_linfct = function() {
      self$Contrast$get_linfct()
    },
    #' @description
    #' get contrasts
    #' @seealso \code{\link{summary_ROPECA_median_p_scaled}}
    #' @param all should all columns be returned (default FALSE)
    #' @param global use a global linear function (determined by get_linfct)
    #' @return data.frame
    get_contrasts = function(all = FALSE) {
      if (is.null(self$contrast_result)) {
        contrast_result <- self$Contrast$get_contrasts(all = FALSE)

        contrast_result <- summary_ROPECA_median_p_scaled(
          contrast_result,
          contrast = "contrast",
          subject_id = self$subject_id[length(self$subject_id) - 1],
          estimate = "diff",
          statistic = "statistic",
          p.value = "p.value",
          max.n = 10
        )
        contrast_result <- dplyr::rename(contrast_result, diff = "estimate")
        contrast_result <- self$p.adjust(
          contrast_result,
          column = "beta.based.significance",
          group_by_col = "contrast",
          newname = "FDR.beta.based.significance"
        )
        contrast_result <- self$p.adjust(
          contrast_result,
          column = "median.p.value",
          group_by_col = "contrast",
          newname = "FDR.median.p.value"
        )
        self$contrast_result <- .stamp_model_identity(contrast_result, self$model_name, "observed")
      }

      if (all) {
        res <- self$contrast_result
      } else {
        res <- select(
          self$contrast_result,
          -all_of(c(
            "n_not_na",
            "mad.estimate",
            "n.beta",
            "isSingular",
            "median.p.scaled",
            "median.p.value",
            "FDR.median.p.value"
          ))
        )
      }

      return(res)
    },
    #' @description
    #' get \code{\link{ContrastsPlotter}}
    #' @return \code{\link{ContrastsPlotter}}
    #' @param fdr_threshold FDR threshold
    #' @param fc_threshold FC threshold
    get_Plotter = function(fdr_threshold = 0.1, fc_threshold = 2) {
      contrast_result <- self$get_contrasts()
      res <- ContrastsPlotter$new(
        contrast_result,
        subject_id = self$subject_id[1],
        fcthresh = fc_threshold,
        volcano = list(list(score = "FDR.beta.based.significance", thresh = fdr_threshold)),
        histogram = list(
          list(score = "beta.based.significance", breaks = seq(0, 1, by = 0.05)),
          list(score = "FDR.beta.based.significance", breaks = seq(0, 1, by = 0.05))
        )
      )
      return(res)
    },
    #' @description convert to wide format
    #' @param columns value column default beta.based.significance
    #' @return data.frame
    to_wide = function(columns = c("beta.based.significance", "FDR.beta.based.significance")) {
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_to_wide(
        contrast_minimal,
        subject_id = self$subject_id[length(self$subject_id) - 1],
        columns = c("diff", columns),
        contrast = "contrast"
      )
      return(contrasts_wide)
    }
  )
)

# ROPECA V2 -----

#' Gauss-Hermite rule for E[f(W)], W ~ N(0, 1), by the Golub-Welsch eigen
#'
#' decomposition of the Jacobi matrix of the probabilists' Hermite polynomials.
#' @param n number of Gauss-Hermite quadrature nodes
#' @noRd
.gauss_hermite_normal <- function(n = 40) {
  jacobi <- matrix(0, n, n)
  off_diagonal <- sqrt(seq_len(n - 1))
  jacobi[cbind(seq_len(n - 1), 2:n)] <- off_diagonal
  jacobi[cbind(2:n, seq_len(n - 1))] <- off_diagonal
  decomposition <- eigen(jacobi, symmetric = TRUE)
  list(nodes = decomposition$values, weights = decomposition$vectors[1, ]^2)
}

#' P(Med_n > u) for the median Med_n of n independent N(0, 1) values.
#'
#' Odd n: Med_n is the order statistic X_(k), k = (n + 1) / 2, so Phi(Med_n) ~ Beta(k, k).
#' Even n: Med_n = (X_(k) + X_(k+1)) / 2, k = n / 2. Given X_(k) = a, X_(k+1) is the
#' smallest of n - k values above a, hence
#'   P(Med_n > u) = P(X_(k) > u) + C * int_{-Inf}^{u} Phi(a)^(k-1) (1 - Phi(2u - a))^(n-k) phi(a) da
#' with C = n! / ((k - 1)! (n - k)!). The integral is tabulated on a u grid and
#' interpolated on the log scale.
#' @param u normal-score threshold
#' @param n number of independent normal observations in each median
#' @noRd
.median_normal_sf <- function(u, n) {
  if (n %% 2 == 1) {
    k <- (n + 1) / 2
    return(stats::pbeta(stats::pnorm(-u), k, k))
  }
  # the u grid is every second point of the a grid, so each upper limit u is a node
  step <- 0.01
  grid_a <- seq(-10, 8, by = step)
  u_index <- seq(201, length(grid_a), by = 2)
  grid_u <- grid_a[u_index]
  k <- n / 2
  log_c <- lgamma(n + 1) - lgamma(k) - lgamma(n - k + 1)
  log_integrand <- outer(
    grid_u,
    grid_a,
    function(u, a) {
      log_c +
        (k - 1) * stats::pnorm(a, log.p = TRUE) +
        (n - k) * stats::pnorm(2 * u - a, lower.tail = FALSE, log.p = TRUE) +
        stats::dnorm(a, log = TRUE)
    }
  )
  integrand <- exp(log_integrand) * outer(u_index, seq_along(grid_a), ">=")
  # trapezoid rule from the first a node to u
  integral <- step * (rowSums(integrand) - integrand[, 1] / 2 - integrand[cbind(seq_along(u_index), u_index)] / 2)
  upper_order_stat <- stats::pbeta(stats::pnorm(-grid_u), n - k + 1, k)
  log_sf <- log(upper_order_stat + integral)
  exp(stats::approx(grid_u, log_sf, xout = u, rule = 2)$y)
}

#' Standard normal score of a t-statistic with df degrees of freedom: the
#'
#' z with the same one-sided tail probability, computed on the log scale so that
#' extreme t stay finite.
#' @param t peptide t-statistics
#' @param df peptide residual degrees of freedom
#' @noRd
.peptide_z <- function(t, df) {
  sign(t) * -stats::qnorm(stats::pt(-abs(t), df, log.p = TRUE), log.p = TRUE)
}

#' Two-sided p-value of the median of n peptide z-scores whose pairwise correlation
#'
#' is rho: M = sqrt(rho) W + sqrt(1 - rho) Med_n, W ~ N(0, 1), p = 2 P(M >= |m|).
#' @param median_z median peptide normal score
#' @param n number of peptide z-scores contributing to each element of median_z
#' @param rho within-protein peptide correlation
#' @noRd
#' @srrstats {G5.4a} Tests include exact fold changes and one-peptide/independent-null ROPECA reductions.
#' @srrstats {G5.6a} Numerical recovery and null-calibration tests use explicit tolerances.
.ropeca_v2_pvalue <- function(median_z, n, rho) {
  rule <- .gauss_hermite_normal()
  p <- rep(NA_real_, length(median_z))
  for (size in unique(n[!is.na(median_z)])) {
    selected <- which(n == size & !is.na(median_z))
    u <- outer(abs(median_z[selected]), sqrt(rho) * rule$nodes, "-") / sqrt(1 - rho)
    tail <- matrix(.median_normal_sf(as.vector(u), size), nrow = length(selected))
    p[selected] <- pmin(1, 2 * as.vector(tail %*% rule$weights))
  }
  p
}

#' Residuals of a fitted peptide model, the samples they belong to, and the
#'
#' residual projection I - Q Q' of the peptide's design over its observed samples.
#' @param model fitted peptide regression model
#' @param data long-format measurement data frame
#' @param sample_col sample identifier column
#' @noRd
.peptide_residual_projection <- function(model, data, sample_col) {
  residual <- stats::residuals(model)
  q <- qr.Q(model$qr)[, seq_len(model$rank), drop = FALSE]
  list(
    sample = as.character(data[[sample_col]][as.integer(names(residual))]),
    residual = unname(residual),
    projection = diag(length(residual)) - q %*% t(q)
  )
}

#' Moment estimates of the variance components of one protein's peptides. With
#'
#' y_j = b + e_j (b shared by the peptides), the residuals r_j = P_j y_j satisfy
#' E<r_j, r_l> = Var(b) * tr(P_j P_l) over the samples both peptides were observed
#' in, and E|r_j|^2 = (Var(b) + Var(e)) * tr(P_j). Using the projections keeps the
#' estimates unbiased when peptides miss different samples. Returns the shared
#' variance Var(b) and the total variance, NA without two overlapping peptides.
#' @param peptides peptide residuals and projection matrices for one protein
#' @noRd
.protein_variance_components <- function(peptides) {
  if (length(peptides) < 2) {
    return(c(shared = NA_real_, total = NA_real_))
  }
  samples <- unique(unlist(lapply(peptides, `[[`, "sample")))
  n_samples <- length(samples)
  residuals <- matrix(0, length(peptides), n_samples)
  projections <- matrix(0, length(peptides), n_samples^2)
  degrees <- numeric(length(peptides))
  for (j in seq_along(peptides)) {
    index <- match(peptides[[j]]$sample, samples)
    residuals[j, index] <- peptides[[j]]$residual
    embedded <- matrix(0, n_samples, n_samples)
    embedded[index, index] <- peptides[[j]]$projection
    projections[j, ] <- as.vector(embedded)
    degrees[j] <- sum(diag(peptides[[j]]$projection))
  }
  pairs <- upper.tri(diag(length(peptides)))
  expected <- sum((projections %*% t(projections))[pairs])
  if (sum(degrees) <= 0 || expected <= 0) {
    return(c(shared = NA_real_, total = NA_real_))
  }
  c(
    shared = sum((residuals %*% t(residuals))[pairs]) / expected,
    total = sum(residuals^2) / sum(degrees)
  )
}

#' Correlation of peptides within a protein, estimated once per dataset from the
#'
#' residuals of the peptide models, pooled over proteins with at least two
#' peptides, clamped to [0, 0.95].
#' @param model_df table with one fitted model and diagnostics per subject
#' @param subject_id column names identifying the modelled subjects
#' @param sample_col sample identifier column
#' @noRd
.ropeca_v2_rho <- function(model_df, subject_id, sample_col) {
  protein_cols <- subject_id[-length(subject_id)]
  fitted <- model_df[model_df$has_model_fit, , drop = FALSE]
  peptides <- purrr::map2(
    fitted$linear_model,
    fitted$data,
    .peptide_residual_projection,
    sample_col = sample_col
  )
  protein_key <- interaction(fitted[protein_cols], drop = TRUE)
  components <- vapply(split(peptides, protein_key), .protein_variance_components, numeric(2))
  usable <- stats::complete.cases(t(components))
  # ratio of the pooled components; a mean of per-protein ratios is biased low,
  # because each protein's total variance shares its sample effect with the numerator
  rho <- sum(components["shared", usable]) / sum(components["total", usable])
  min(0.95, max(0, rho))
}

#' Require one peptide correlation in the interval [0, 1)
#' @param rho within-protein peptide correlation
#' @noRd
.assert_peptide_correlation <- function(rho) {
  in_range <- is.numeric(rho) && length(rho) == 1 && isTRUE(rho >= 0 && rho < 1)
  if (!in_range) {
    abort_bad_argument("rho", "be a single number in [0, 1)")
  }
}

#' Protein-level ROPECA V2 contrasts from peptide-level contrasts with columns diff,
#'
#' statistic, df and avgAbd.
#' @param peptide_contrasts peptide-level contrast table
#' @param protein_cols column names identifying parent proteins
#' @param rho within-protein peptide correlation
#' @noRd
.ropeca_v2_protein_contrasts <- function(peptide_contrasts, protein_cols, rho) {
  # rows with NA (or NaN) diff, statistic or df are dropped here, so mad() and median() below see no NA
  peptides <- peptide_contrasts |>
    dplyr::filter(!is.na(.data$diff), !is.na(.data$statistic), !is.na(.data$df)) |>
    dplyr::mutate(z = .peptide_z(.data$statistic, .data$df))
  # sigma is computed before diff is replaced by its median
  res <- peptides |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(protein_cols, "contrast")))) |>
    dplyr::summarise(
      avgAbd = stats::median(.data$avgAbd, na.rm = TRUE),
      sigma = stats::mad(.data$diff),
      diff = stats::median(.data$diff),
      median_z = stats::median(.data$z),
      n_peptides = dplyr::n(),
      .groups = "drop"
    )
  res$p.value <- .ropeca_v2_pvalue(res$median_z, res$n_peptides, rho)
  res$statistic <- sign(res$median_z) * stats::qnorm(res$p.value / 2, lower.tail = FALSE)
  res$std.error <- ifelse(res$statistic != 0, abs(res$diff / res$statistic), NA_real_)
  res$df <- Inf
  res$conf.low <- res$diff - stats::qnorm(0.975) * res$std.error
  res$conf.high <- res$diff + stats::qnorm(0.975) * res$std.error
  res$rho <- rho
  res
}

#' ROPECA V2: median peptide z-score with a correlation-aware null
#'
#' Summarises peptide-level contrasts per protein like \code{\link{ContrastsROPECA}},
#' but tests the median peptide z-score against its exact null distribution: the
#' median of n standard normal scores with pairwise correlation \code{rho}. The
#' protein p-values are uniform under the null for any number of peptides, and
#' correlated peptides do not make proteins with many peptides more significant.
#' \code{\link{ContrastsROPECAV2NestedFacade}} estimates \code{rho} from the peptide
#' model residuals.
#'
#' Each peptide t-statistic is converted to the standard normal score with the
#' same tail probability. \code{diff} is the median peptide fold change,
#' \code{statistic} the normal score with the protein's two-sided p-value,
#' \code{std.error} = \code{|diff / statistic|}, \code{df} = \code{Inf}, and
#' \code{conf.low}, \code{conf.high} the matching normal 95\% interval;
#' \code{sigma} is the MAD of the peptide fold changes.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_peptide_config(Nprot = 20)
#' config <- istar$config$clone(deep = TRUE)
#' config$hierarchy_depth <- 2
#' mod <- build_model(
#'   istar$data,
#'   strategy_lm("abundance ~ group_"),
#'   subject_id = config$hierarchy_keys_depth()
#' )
#' contrasts <- c("AvsCtrl" = "group_A - group_Ctrl")
#' peptide_contrasts <- ContrastsModerated$new(Contrasts$new(mod, contrasts))
#' ropeca2 <- ContrastsROPECAV2$new(peptide_contrasts, rho = 0.2)
#' head(ropeca2$get_contrasts())
#' ropeca2$to_wide()
#' pl <- ropeca2$get_Plotter()
#' pl$volcano()
ContrastsROPECAV2 <- R6::R6Class(
  "ContrastsROPECAV2",
  inherit = ContrastsInterface,
  public = list(
    #' @field Contrast peptide-level contrasts, e.g. \code{ContrastsModerated}
    Contrast = NULL,
    #' @field rho correlation of the peptides of a protein
    rho = NULL,
    #' @field contrast_result contrast result
    contrast_result = NULL,
    #' @field model_name model name
    model_name = character(),
    #' @field subject_id peptide-level subject columns; the protein columns are all but the last
    subject_id = character(),
    #' @field p.adjust method to use for p.value adjustment
    p.adjust = NULL,
    #' @description
    #' initialize
    #' @param Contrast peptide-level contrasts with columns \code{diff}, \code{statistic} and \code{df}
    #' @param rho correlation of the peptides of a protein, in [0, 1)
    #' @param model_name default ROPECAV2
    #' @param p.adjust function to use for p.value adjustment
    initialize = function(Contrast, rho, model_name = "ROPECAV2", p.adjust = prolfqua::adjust_p_values) {
      if (length(Contrast$subject_id) < 2) {
        abort_bad_argument("Contrast", "have peptide-level subject_id (parent and child keys)")
      }
      .assert_peptide_correlation(rho)
      self$Contrast <- Contrast
      self$rho <- rho
      self$model_name <- model_name
      self$subject_id <- Contrast$subject_id
      self$p.adjust <- p.adjust
    },
    #' @description
    #' show names of contrasts
    #' @return data.frame
    get_contrast_sides = function() {
      self$Contrast$get_contrast_sides()
    },
    #' @description
    #' get linear function used to determine contrasts
    #' @return data.frame
    get_linfct = function() {
      self$Contrast$get_linfct()
    },
    #' @description
    #' get protein-level contrasts
    #' @param all \code{FALSE} (default) returns the standard columns;
    #'   \code{TRUE} appends \code{n_peptides}, \code{median_z} and \code{rho}
    #' @return data.frame
    get_contrasts = function(all = FALSE) {
      protein_cols <- private$protein_cols()
      if (is.null(self$contrast_result)) {
        res <- .ropeca_v2_protein_contrasts(self$Contrast$get_contrasts(), protein_cols, self$rho)
        res <- self$p.adjust(res, column = "p.value", group_by_col = "contrast", newname = "FDR")
        self$contrast_result <- .stamp_model_identity(dplyr::ungroup(res), self$model_name, "observed")
      }
      standard <- c(
        protein_cols,
        "modelName",
        "estimate_type",
        "contrast",
        "avgAbd",
        "diff",
        "FDR",
        "statistic",
        "std.error",
        "df",
        "p.value",
        "conf.low",
        "conf.high",
        "sigma"
      )
      if (all) {
        return(dplyr::relocate(self$contrast_result, dplyr::all_of(standard)))
      }
      self$contrast_result[, standard]
    },
    #' @description
    #' get \code{\link{ContrastsPlotter}}
    #' @return \code{\link{ContrastsPlotter}}
    #' @param fdr_threshold FDR threshold
    #' @param fc_threshold FC threshold
    get_Plotter = function(fdr_threshold = 0.1, fc_threshold = 2) {
      ContrastsPlotter$new(
        self$get_contrasts(),
        subject_id = private$protein_cols(),
        fcthresh = fc_threshold,
        volcano = list(list(score = "FDR", thresh = fdr_threshold)),
        histogram = list(
          list(score = "p.value", breaks = seq(0, 1, by = 0.05)),
          list(score = "FDR", breaks = seq(0, 1, by = 0.05))
        )
      )
    },
    #' @description convert to wide format
    #' @param columns value columns, default \code{c("p.value", "FDR")}
    #' @return data.frame
    to_wide = function(columns = c("p.value", "FDR")) {
      pivot_model_contrasts_to_wide(
        self$get_contrasts(),
        subject_id = private$protein_cols(),
        columns = c("diff", columns),
        contrast = "contrast"
      )
    }
  ),
  private = list(
    protein_cols = function() {
      self$subject_id[-length(self$subject_id)]
    }
  )
)
