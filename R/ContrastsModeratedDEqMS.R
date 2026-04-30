# ContrastsModeratedDEqMS -----

#' Find prior degrees of freedom for DEqMS
#'
#' Uses trigammaInverse to estimate d0 from the mean residual variance
#' after removing the count-dependent trend.
#'
#' @param mean_myfct mean of squared residuals minus trigamma correction
#' @return scalar prior degrees of freedom (d0)
#' @keywords internal
find_d0_deqms <- function(mean_myfct) {
  if (is.na(mean_myfct) || mean_myfct <= 0) {
    return(Inf)
  }
  d0 <- 2 * trigammaInverse(mean_myfct)
  return(max(d0, 0.1))
}


#' DEqMS-style count-dependent variance moderation
#'
#' Applies count-dependent empirical Bayes variance shrinkage to a contrast
#' result table. Proteins quantified from many peptides get less shrinkage;
#' proteins from few peptides get more.
#'
#' @param mm data.frame from one contrast group with columns: sigma, df,
#'   statistic, std.error, and the estimate column
#' @param count_col name of column with peptide/PSM count per protein
#' @param df name of the degrees of freedom column
#' @param estimate name of the fold change column
#' @param loess_span span parameter for LOESS fit (default 0.75)
#' @param confint confidence level for intervals (default 0.95)
#' @return data.frame with added moderated.* columns
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#' mm <- data.frame(
#'   sigma = c(0.25, 0.32, 0.28, 0.40, 0.35),
#'   df = rep(6, 5),
#'   statistic = c(2.1, -1.8, 0.5, 3.0, -2.2),
#'   diff = c(0.8, -0.6, 0.2, 1.2, -0.9),
#'   count = c(2, 3, 4, 6, 8)
#' )
#' res <- moderated_p_deqms(mm, count_col = "count")
#' "moderated.p.value" %in% colnames(res)
moderated_p_deqms <- function(mm, count_col, df = "df", estimate = "diff", loess_span = 0.75, confint = 0.95) {
  logvar <- log(mm$sigma^2)
  log2count <- log2(mm[[count_col]])
  df_res <- mm[[df]]

  # Handle edge cases: remove NAs/Infs for fitting
  ok <- is.finite(logvar) & is.finite(log2count) & is.finite(df_res) & df_res > 0
  if (sum(ok) < 4) {
    # Too few observations for LOESS — fall back to standard limma moderation
    return(moderated_p_limma(mm, df = df, estimate = estimate, confint = confint))
  }

  # Fit variance ~ count curve using LOESS
  loess_fit <- stats::loess(logvar[ok] ~ log2count[ok], span = loess_span)

  # Predict for all proteins (including those not used in fit)
  fitted_logvar <- rep(NA_real_, length(logvar))
  fitted_logvar[ok] <- stats::fitted(loess_fit)
  # For non-ok, predict from the fit if possible
  not_ok <- !ok & is.finite(log2count)
  if (any(not_ok)) {
    fitted_logvar[not_ok] <- predict(
      loess_fit,
      newdata = data.frame(
        `log2count[ok]` = log2count[not_ok]
      )
    )
  }
  # For any remaining NAs, use the mean fitted value
  fitted_logvar[is.na(fitted_logvar)] <- mean(fitted_logvar, na.rm = TRUE)

  # Estimate d0 via trigamma approach
  eg <- logvar - digamma(df_res / 2) + log(df_res / 2)
  egpred <- fitted_logvar - digamma(df_res / 2) + log(df_res / 2)
  myfct <- (eg - egpred)^2 - trigamma(df_res / 2)
  mean_myfct <- mean(myfct[ok], na.rm = TRUE)

  d0 <- find_d0_deqms(mean_myfct)

  # Count-specific prior variance
  if (is.finite(d0)) {
    s02 <- exp(egpred + digamma(d0 / 2) - log(d0 / 2))
  } else {
    s02 <- exp(fitted_logvar)
  }

  # Posterior variance
  if (is.finite(d0)) {
    var_post <- (d0 * s02 + df_res * mm$sigma^2) / (d0 + df_res)
    df_total <- d0 + df_res
  } else {
    var_post <- s02
    df_total <- df_res
  }

  # Moderated statistics
  moderated_t <- mm$statistic * mm$sigma / sqrt(var_post)
  moderated_p <- 2 * pt(abs(moderated_t), df = df_total, lower.tail = FALSE)

  # Confidence intervals (same pattern as moderated_p_limma)
  prqt <- -qt((1 - confint) / 2, df = df_total)
  conf_low <- mm[[estimate]] - prqt * sqrt(var_post)
  conf_high <- mm[[estimate]] + prqt * sqrt(var_post)

  # Attach results (same naming convention as moderated_p_limma)
  mm$moderated.var.prior <- s02
  mm$moderated.df.prior <- d0
  mm$moderated.var.post <- var_post
  mm$moderated.statistic <- moderated_t
  mm$moderated.df.total <- df_total
  mm$moderated.p.value <- moderated_p
  mm$moderated.conf.low <- conf_low
  mm$moderated.conf.high <- conf_high

  return(mm)
}


#' DEqMS-style moderation for long contrast table
#'
#' Splits by contrast group and applies \code{\link{moderated_p_deqms}} to
#' each group independently.
#'
#' @param mm contrast result data.frame (long format)
#' @param count_col name of column with peptide/PSM count per protein
#' @param group_by_col column to group by (default "contrast")
#' @param estimate name of the fold change column
#' @param loess_span span parameter for LOESS fit (default 0.75)
#' @return data.frame with moderated columns added
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#' mm <- data.frame(
#'   contrast = rep(c("A_vs_B", "C_vs_D"), each = 5),
#'   sigma = rep(c(0.25, 0.32, 0.28, 0.40, 0.35), 2),
#'   df = rep(6, 10),
#'   statistic = rep(c(2.1, -1.8, 0.5, 3.0, -2.2), 2),
#'   diff = rep(c(0.8, -0.6, 0.2, 1.2, -0.9), 2),
#'   count = rep(c(2, 3, 4, 6, 8), 2)
#' )
#' res <- moderated_p_deqms_long(mm, count_col = "count")
#' table(res$contrast)
moderated_p_deqms_long <- function(mm, count_col, group_by_col = "contrast", estimate = "diff", loess_span = 0.75) {
  dfg <- mm |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_by_col))) |>
    tidyr::nest()

  results <- vector("list", nrow(dfg))
  warning_rows <- list()

  for (i in seq_len(nrow(dfg))) {
    warning_env <- new.env(parent = emptyenv())
    warning_env$messages <- character()
    group_data <- dfg[i, setdiff(names(dfg), "data"), drop = FALSE]
    result_i <- withCallingHandlers(
      moderated_p_deqms(
        dfg$data[[i]],
        count_col = count_col,
        estimate = estimate,
        loess_span = loess_span
      ),
      warning = function(w) {
        warning_env$messages <- c(warning_env$messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    warnings_i <- warning_env$messages
    results[[i]] <- dplyr::bind_cols(
      group_data[rep(1, nrow(result_i)), , drop = FALSE],
      result_i
    )

    if (length(warnings_i) > 0) {
      group_label <- paste(
        sprintf(
          "%s=%s",
          names(group_data),
          unlist(group_data[1, ], use.names = FALSE)
        ),
        collapse = ", "
      )
      if (!nzchar(group_label)) {
        group_label <- paste0("group_", i)
      }
      warning_rows[[length(warning_rows) + 1]] <- tibble::tibble(
        group = group_label,
        messages = paste(unique(warnings_i), collapse = "; ")
      )
    }
  }

  xx <- dplyr::bind_rows(results)

  if (length(warning_rows) > 0) {
    warning_df <- dplyr::bind_rows(warning_rows)
    n_show <- min(nrow(warning_df), 3)
    warning_examples <- paste(
      vapply(
        seq_len(n_show),
        function(i) {
          paste0(warning_df$group[[i]], " (", warning_df$messages[[i]], ")")
        },
        character(1)
      ),
      collapse = "; "
    )
    msg <- sprintf(
      "moderated_p_deqms_long: condition messages in %d/%d groups. %s",
      nrow(warning_df),
      nrow(dfg),
      warning_examples
    )
    warning(msg, call. = FALSE)
  }

  return(xx)
}


#' DEqMS count-dependent moderated contrasts
#'
#' Decorator that wraps any Contrasts object and applies count-dependent
#' empirical Bayes variance shrinkage. Similar to \code{\link{ContrastsModerated}}
#' but the prior variance depends on the number of quantified peptides/PSMs
#' per protein: proteins with many peptides get less shrinkage, proteins with
#' few peptides get more.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#'
#' istar <- sim_lfq_data_protein_config(Nprot = 50)
#' protIntensity <- istar$data
#' config <- istar$config
#'
#' lProt <- LFQData$new(protIntensity, config)
#' lProt$rename_response("transformedIntensity")
#' modelFunction <-
#'   strategy_lm("transformedIntensity ~ group_")
#' mod <- build_model(
#'  lProt,
#'  modelFunction)
#'
#' Contr <- c("dil.b_vs_a" = "group_A - group_Ctrl")
#' contrast <- prolfqua::Contrasts$new(mod, Contr)
#'
#' # Build count_df from config
#' count_df <- dplyr::select(protIntensity,
#'   dplyr::all_of(c(config$hierarchy_keys_depth(), "nr_peptides"))) |>
#'   dplyr::distinct()
#'
#' deqms <- ContrastsModeratedDEqMS$new(contrast,
#'   count_df = count_df,
#'   count_column = "nr_peptides")
#'
#' bb <- deqms$get_contrasts()
#' stopifnot(all(c("diff", "p.value", "FDR", "sigma") %in% colnames(bb)))
#'
#' # Merge with ContrastsMissing
#' csi <- ContrastsMissing$new(lProt, contrasts = Contr)
#' merged <- merge_contrasts_results(deqms, csi)
#'
#' cs <- deqms$get_contrast_sides()
#' cslf <- deqms$get_linfct()
#' ctrwide <- deqms$to_wide()
#' cp <- deqms$get_Plotter()
#' cp$volcano()
#'
ContrastsModeratedDEqMS <- R6::R6Class(
  classname = "ContrastsModeratedDEqMS",
  inherit = ContrastsInterface,
  public = list(
    #' @field Contrast Class implementing the Contrast interface
    Contrast = NULL,
    #' @field count_df data.frame with subject_id + count column
    count_df = NULL,
    #' @field count_column name of the count column in count_df
    count_column = character(),
    #' @field loess_span span parameter for LOESS fit
    loess_span = numeric(),
    #' @field model_name name of model
    model_name = character(),
    #' @field subject_id columns with subject_id (proteinID)
    subject_id = character(),
    #' @field p.adjust function to adjust p-values
    p.adjust = NULL,
    #' @description
    #' initialize
    #' @param Contrast class implementing the ContrastInterface
    #' @param count_df data.frame with subject_id columns and a count column
    #' @param count_column name of the count column in count_df
    #' @param loess_span span for LOESS variance fit (default 0.75)
    #' @param model_name name of the model
    #' @param p.adjust function to adjust p-values - default BH
    initialize = function(
      Contrast,
      count_df,
      count_column,
      loess_span = 0.75,
      model_name = paste0(Contrast$model_name, "_DEqMS"),
      p.adjust = prolfqua::adjust_p_values
    ) {
      self$Contrast <- Contrast
      self$subject_id <- Contrast$subject_id
      # Aggregate count_df to one row per subject_id (max count)
      count_df <- count_df |>
        dplyr::group_by(dplyr::across(dplyr::all_of(Contrast$subject_id))) |>
        dplyr::summarise(!!count_column := max(!!rlang::sym(count_column), na.rm = TRUE), .groups = "drop")
      # Replace -Inf (from all-NA groups) with NA
      count_df[[count_column]][is.infinite(count_df[[count_column]])] <- NA
      self$count_df <- count_df
      self$count_column <- count_column
      self$loess_span <- loess_span
      self$model_name <- model_name
      self$p.adjust <- p.adjust
    },
    #' @description
    #' get both sides of contrasts
    get_contrast_sides = function() {
      self$Contrast$get_contrast_sides()
    },
    #' @description
    #' get linear functions from contrasts
    #' @param global logical TRUE - get linear functions for all models
    get_linfct = function(global = TRUE) {
      self$Contrast$get_linfct()
    },
    #' @description
    #' applies DEqMS-style count-dependent moderation
    #' @seealso \code{\link{moderated_p_deqms_long}}
    #' @param all should all columns be returned (default FALSE)
    get_contrasts = function(all = FALSE) {
      contrast_result <- self$Contrast$get_contrasts(all = FALSE)

      # Join count_df to get per-protein peptide counts
      contrast_result <- dplyr::inner_join(
        contrast_result,
        self$count_df,
        by = self$subject_id
      )

      # Ensure count >= 1 to avoid log2(0)
      contrast_result[[self$count_column]] <- pmax(
        contrast_result[[self$count_column]],
        1
      )

      contrast_result <- moderated_p_deqms_long(
        contrast_result,
        count_col = self$count_column,
        group_by_col = "contrast",
        estimate = "diff",
        loess_span = self$loess_span
      )

      # Remove the count column from output (it was joined for computation)
      contrast_result[[self$count_column]] <- NULL

      if (!all) {
        contrast_result <- contrast_result |>
          dplyr::select(
            -c(
              "sigma",
              "df",
              "statistic",
              "p.value",
              "conf.low",
              "conf.high",
              "FDR",
              "moderated.df.prior",
              "moderated.var.prior"
            )
          )
        contrast_result <- contrast_result |>
          dplyr::mutate(sigma = sqrt(moderated.var.post), .keep = "unused")
        contrast_result <- contrast_result |>
          dplyr::rename(
            conf.low = "moderated.conf.low",
            conf.high = "moderated.conf.high",
            statistic = "moderated.statistic",
            df = "moderated.df.total",
            p.value = "moderated.p.value"
          )
        contrast_result <- self$p.adjust(
          contrast_result,
          column = "p.value",
          group_by_col = "contrast",
          newname = "FDR"
        )
      } else {
        contrast_result <- self$p.adjust(
          contrast_result,
          column = "moderated.p.value",
          group_by_col = "contrast",
          newname = "FDR.moderated"
        )
      }

      contrast_result <- dplyr::ungroup(contrast_result)

      if (inherits(contrast_result$modelName, "factor")) {
        mname <- factor(
          paste0(contrast_result$modelName, "_DEqMS"),
          levels = paste0(levels(contrast_result$modelName), "_DEqMS")
        )
      } else {
        mname <- paste0(contrast_result$modelName, "_DEqMS")
      }
      contrast_result$modelName <- mname

      stopifnot(all(
        super$column_description()$column_name %in% colnames(contrast_result)
      ))

      return(contrast_result)
    },
    #' @description
    #' get \code{\link{ContrastsPlotter}}
    #' @param fc_threshold fold change threshold to show in plots
    #' @param fdr_threshold FDR threshold to show in plots
    get_Plotter = function(
      fc_threshold = 1,
      fdr_threshold = 0.1
    ) {
      contrast_result <- self$get_contrasts()
      res <- ContrastsPlotter$new(
        contrast_result,
        subject_id = self$subject_id,
        fcthresh = fc_threshold,
        volcano = list(list(score = "FDR", thresh = fdr_threshold)),
        histogram = list(
          list(score = "p.value", xlim = c(0, 1, 0.05)),
          list(score = "FDR", xlim = c(0, 1, 0.05))
        ),
        score = list(list(score = "statistic", thresh = 5)),
        modelName = "modelName",
        diff = "diff",
        contrast = "contrast"
      )
      return(res)
    },
    #' @description
    #' convert to wide format
    #' @param columns value column default p.value, FDR, statistic
    #' @return data.frame
    to_wide = function(columns = c("p.value", "FDR", "statistic")) {
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_to_wide(
        contrast_minimal,
        subject_id = self$subject_id,
        columns = c("diff", columns),
        contrast = "contrast"
      )
      return(contrasts_wide)
    }
  )
)
