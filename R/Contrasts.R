.linfct <- function(model, contrast, avg = TRUE) {
  linfct <- linfct_from_model(model, as_list = FALSE)
  linfct <- unique(linfct) # needed for single factor models
  if (avg) {
    namtmp <- paste0("avg_", names(contrast))
    cntr_avg <- paste0("(", gsub(" - ", " + ", contrast), ")/2")
    names(cntr_avg) <- namtmp
    contrast <- c(contrast, cntr_avg)
  }
  linfct_a <- linfct_matrix_contrasts(linfct, contrast)
  return(linfct_a)
}

.linfct_partial_model <- function(model, contrast, avg = TRUE) {
  withCallingHandlers(
    .linfct(model, contrast = contrast, avg = avg),
    warning = function(w) invokeRestart("muffleWarning")
  )
}

# One linfct per model. With `reuse_complete`, models with the maximal number of estimated coefficients share the
# linfct of the complete model fit; all other models get their own from `linfct_fun`.
.linfct_per_model <- function(model_df, contrasts, avg, label, linfct_fun = .linfct, reuse_complete = TRUE) {
  res <- vector(mode = "list", nrow(model_df))
  pb <- .make_progress(nrow(model_df), label = label)
  if (reuse_complete && nrow(model_df) > 0) {
    compmodel <- .linfct(get_complete_model_fit(model_df)$linear_model[[1]], contrasts, avg = avg)
    max_coef <- max(model_df$nr_coef_not_NA)
  }
  for (i in seq_along(model_df$linear_model)) {
    pb$tick()
    res[[i]] <- if (reuse_complete && model_df$nr_coef_not_NA[[i]] == max_coef) {
      compmodel
    } else {
      linfct_fun(model_df$linear_model[[i]], contrast = contrasts, avg = avg)
    }
  }
  res
}

# Wald-test finish shared by Contrasts and ContrastsFirth: name the contrast columns, join the `avg_` rows back as
# avgAbd, adjust p-values and keep the unmoderated standard error and df.
.finalize_wald_contrasts <- function(contrast_result, contrasts, subject_id, p.adjust) {
  contrast_result <- dplyr::rename(ungroup(contrast_result), contrast = "lhs", diff = "estimate")

  differences <- contrast_result |>
    dplyr::filter(.data$contrast %in% names(contrasts))

  avg_abd <- contrast_result |>
    dplyr::select(dplyr::all_of(c(subject_id, "contrast", "diff"))) |>
    dplyr::filter(startsWith(.data$contrast, "avg_"))

  avg_abd$contrast <- gsub("^avg_", "", avg_abd$contrast)
  avg_abd <- avg_abd |> dplyr::rename(avgAbd = "diff")
  contrast_result <- left_join(differences, avg_abd)

  contrast_result <- p.adjust(contrast_result, column = "p.value", group_by_col = "contrast")
  contrast_result <- contrast_result |> relocate("FDR", .after = "diff")
  contrast_result |>
    dplyr::mutate(
      std.error.unmoderated = .data$std.error,
      df.unmoderated = .data$df
    )
}

.select_wald_columns <- function(contrast_result, all) {
  if (all) {
    return(contrast_result)
  }
  dplyr::select(contrast_result, -dplyr::all_of(c("sigma.model", "df.residual.model", "isSingular")))
}

# Contrasts -----

#' Estimate contrasts using Wald Test
#'
#' The per-protein contrast computation uses \code{\link{compute_contrast}} and
#' \code{\link{linfct_matrix_contrasts}} internally.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @examples
#'
#' # Fitting mixed effects model to peptide data
#' istar <- prolfqua::sim_lfq_data_peptide_config()
#'
#' modelFunction <-
#'   strategy_lmer("abundance  ~ group_ + (1 | peptide_Id) + (1 | sample)")
#'
#' config <- istar$config
#' config$hierarchy_keys_depth()
#'
#' mod <- build_model(
#'   istar$data,
#'   modelFunction,
#'   subject_id = config$hierarchy_keys_depth()
#' )
#'
#' ref_lfc <- data.frame(
#'   `(Intercept)` = c(0, 0, 0),
#'   group_B = c(0, 1, 0),
#'   group_Ctrl = c(-1, -1, 1),
#'   row.names = c("groupA_vs_Ctrl", "dil.e_vs_b", "dil.ctrl_vs_b")
#' )
#' prolfqua::model_summary(mod)
#' Contr <- c(
#'   "groupA_vs_Ctrl" = "group_A - group_Ctrl",
#'   "dil.e_vs_b" = "group_B - group_Ctrl",
#'   "dil.ctrl_vs_b" = "group_Ctrl - group_A"
#' )
#' contrastX <- prolfqua::Contrasts$new(mod, Contr)
#' y <- contrastX$get_linfct(avg = FALSE)
#' stopifnot(all(ref_lfc == y))
#' t <- contrastX$get_linfct(global = FALSE)
#'
#' x <- contrastX$get_contrasts()
#' stopifnot(all(x$p.value < 1 & x$p.value > 0))
#' stopifnot(all(x$avgAbd > 0))
#' stopifnot(all(x$FDR > 0 & x$FDR < 1))
#'
#' x <- contrastX$get_contrast_sides()
#' xd <- contrastX$column_description()
#' modelFunction <-
#'   strategy_lm("abundance  ~ group_")
#' mod <- build_model(
#'   istar$data,
#'   modelFunction,
#'   subject_id = config$hierarchy_keys_depth()
#' )
#' contrastX <- prolfqua::Contrasts$new(mod, Contr)
#' y <- contrastX$get_linfct(avg = FALSE)
#' stopifnot(all(ref_lfc == y))
#'
Contrasts <- R6::R6Class(
  "Contrast",
  inherit = ContrastsInterface,
  public = list(
    #' @field models Model
    models = NULL,
    #' @field contrasts character with contrasts
    contrasts = character(),
    #' @field contrastfun function to compute contrasts
    contrastfun = NULL,
    #' @field model_name model name
    model_name = character(),
    #' @field subject_id name of column containing e.g., protein Id's
    subject_id = character(),
    #' @field p.adjust function to adjust p-values (default prolfqua::adjust_p_values)
    p.adjust = NULL,
    #' @field contrast_result data frame containing results of contrast computation
    contrast_result = NULL,
    #' @field global use a global linear function (determined by get_linfct)
    global = TRUE,
    #' @description
    #' initialize
    #' create Contrast
    #' @param model a dataframe with a structure similar to that generated by \code{\link{build_model}}
    #' @param contrasts a character vector with contrast specificiation
    #' @param p.adjust function to adjust the p-values
    #' @param global development/internal argument (if FALSE determine linfct for each model.)
    #' @param model_name name of contrast method, default WaldTest
    initialize = function(
      model,
      contrasts,
      p.adjust = prolfqua::adjust_p_values,
      global = FALSE,
      model_name = "WaldTest"
    ) {
      self$models <- model$model_df |> dplyr::filter(has_model_fit == TRUE)
      self$contrasts <- contrasts
      self$contrastfun <- model$model_strategy$contrast_fun
      self$model_name <- model_name
      self$subject_id <- model$subject_id
      self$p.adjust <- p.adjust
      self$global <- global
    },
    #' @description
    #' get both sides of contrasts
    get_contrast_sides = function() {
      parse_contrast_sides(self$contrasts)
    },
    #' @description
    #' get linear functions from contrasts
    #' @param global logical TRUE - get the a linear functions for all models, FALSE - linear function for each model
    #' @param avg logical TRUE - get also linfct for averages
    get_linfct = function(global = TRUE, avg = TRUE) {
      if (global) {
        model <- get_complete_model_fit(self$models)$linear_model[[1]]
        res <- .linfct(model, self$contrasts, avg = avg)
        return(res)
      }
      .linfct_per_model(self$models, self$contrasts, avg, "linfct", linfct_fun = .linfct_partial_model)
    },
    #' @description
    #' get table with contrast estimates
    #' @param all should all columns be returned (default FALSE)
    #' @return data.frame with contrasts
    #'
    get_contrasts = function(all = FALSE) {
      if (is.null(self$contrast_result)) {
        message("determine linear functions:")
        linfct <- self$get_linfct(global = self$global)
        message("get_contrasts -> contrasts_linfct")
        # TODO (goes into calling code)
        contrast_result <- contrasts_linfct(
          self$models,
          linfct,
          subject_id = self$subject_id,
          contrastfun = self$contrastfun
        )
        contrast_result <- .finalize_wald_contrasts(contrast_result, self$contrasts, self$subject_id, self$p.adjust)
        # Stamp modelName uniformly with the model identity. Rescue/imputation
        # state lives in a separate `estimate_type` column: rows refit from an
        # LOD-imputed model (flagged by impute_refit_singular) get
        # "lod_imputed", everything else "observed".
        estimate_type <- "observed"
        if ("imputed" %in% colnames(self$models)) {
          imputed_lookup <- self$models |>
            dplyr::select(dplyr::all_of(c(self$subject_id, "imputed")))
          contrast_result <- dplyr::left_join(
            contrast_result,
            imputed_lookup,
            by = self$subject_id
          )
          estimate_type <- ifelse(
            !is.na(contrast_result$imputed) & contrast_result$imputed,
            "lod_imputed",
            "observed"
          )
          contrast_result$imputed <- NULL
        }
        self$contrast_result <- .stamp_model_identity(contrast_result, self$model_name, estimate_type)
      }
      res <- .select_wald_columns(self$contrast_result, all)
      stopifnot(all(super$column_description()$column_name %in% colnames(res)))
      return(res)
    }
  )
)
