#' compute contrasts
#'
#'
#' @export
#' @keywords internal
#' @examples
#'
#' mod3 <- sim_build_models_logistf(model = "parallel3", weight_missing = 1, peptide=TRUE)
#' contrasts <- c(Avs = "group_A - group_B", AvsCtrl = "group_A - group_Ctrl")
#' ctrpep <- ContrastsFirth$new(mod3,contrasts)
#' ctrpep$get_contrast_sides()
#'
#' linfct_models <- ctrpep$get_linfct()
#' tmp1 <- contrasts_linfct_firth(linfct_models$models1)
#' tmp2 <- contrasts_linfct_firth(linfct_models$models2)
#' stopifnot(all(dim(tmp1) > 10))
#' stopifnot(all(dim(tmp2) > 10))
contrasts_linfct_firth <- function(models, subject_id = "protein_Id") {
  model_df <- models$model_df |>
    dplyr::filter(
      .data$has_model_fit,
      !is.na(.data$nr_coef_not_NA),
      !purrr::map_lgl(.data$linear_model, is.character)
    )
  #computeGroupAverages
  message("contrasts_linfct_firth")
  if (nrow(model_df) == 0) {
    return(tibble::tibble())
  }
  modelcol <- "linear_model"

  interaction_models <- vector(mode = "list", length = nrow(model_df))
  pb <- .make_progress(length(model_df[[modelcol]]), label = "firth contrasts")

  for (i in seq_along(model_df[[modelcol]])) {
    # nolint start: object_usage_linter
    interaction_models[[i]] <- tryCatch(
      .compute_contrast(
        model_df[[modelcol]][[i]],
        linfct = model_df$linfct[[i]],
        strategy = models$strategy
      ),
      error = function(e) FALSE
    )
    # nolint end
    pb$tick()
  }

  interaction_model_matrix <- model_df
  interaction_model_matrix$contrast <- interaction_models

  mclass <- function(x) {
    class(x)[1]
  }

  interaction_model_matrix <- interaction_model_matrix |>
    dplyr::mutate(classC = purrr::map_chr(.data$contrast, mclass))

  n_failed <- sum(interaction_model_matrix$classC == "logical")
  if (n_failed > 0) {
    message(
      "contrasts_linfct_logistf: dropped ",
      n_failed,
      " of ",
      nrow(interaction_model_matrix),
      " proteins with failed contrasts."
    )
  }

  interaction_model_matrix <- interaction_model_matrix |>
    dplyr::filter(.data$classC != "logical")

  contrasts <- interaction_model_matrix |>
    dplyr::select(all_of(c(subject_id, "contrast"))) |>
    tidyr::unnest(cols = c("contrast"))

  # take sigma and df from somewhere else.
  model_infos <- model_df |>
    dplyr::select(all_of(c(subject_id, "isSingular", "sigma.model" = "sigma", "df.residual.model" = "df.residual"))) |>

    dplyr::distinct()
  contrasts <- dplyr::inner_join(contrasts, model_infos, by = subject_id)
  return(ungroup(contrasts))
}


.prepare_logistf_lfqdata <- function(lfqdata) {
  stopifnot("LFQData" %in% class(lfqdata))
  lfq_missing <- lfqdata$get_copy()
  lfq_missing$complete_cases()
  lfq_missing$set_data(prolfqua::encode_bin_resp(lfq_missing))
  lfq_missing$set_config_value("bin_resp", "bin_resp")
  lfq_missing
}


#' Build Firth logistic model for aggregated LFQData
#'
#' Encodes missingness as a binary response and fits the Firth logistic backend
#' used by the missingness model path in \code{prolfquapp}.
#'
#' @param lfqdata aggregated \code{\link{LFQData}} object
#' @param modelstr model formula string without the response variable
#'   (e.g. \code{"~ group_"})
#' @return a \code{\link{ModelFirth}} object
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_protein_config(
#'   Nprot = 10, with_missing = TRUE, weight_missing = 0.5, seed = 3)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' mod <- build_model_glm_protein(lfqdata, "~ group_")
#' head(mod$get_coefficients())
build_model_glm_protein <- function(lfqdata, modelstr) {
  .assert_aggregated_facade_input(lfqdata, "build_model_glm_protein")
  lfq_missing <- .prepare_logistf_lfqdata(lfqdata)
  formula <- paste(lfq_missing$get_config()$bin_resp, modelstr)
  build_model_logistf(lfq_missing, formula)
}


#' Build Firth logistic model for nested LFQData
#'
#' Encodes missingness as a binary response and fits the peptide-aware Firth
#' logistic backend. Proteins with multiple child features are fitted with the
#' lowest hierarchy key appended to the formula.
#'
#' @param lfqdata nested \code{\link{LFQData}} object
#' @param modelstr model formula string without the response variable
#'   (e.g. \code{"~ group_"})
#' @return a \code{\link{ModelFirth}} object
#' @export
#' @family modelling
#' @examples
#' istar <- sim_lfq_data_peptide_config(
#'   Nprot = 10, with_missing = TRUE, weight_missing = 0.5, seed = 3)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' mod <- build_model_glm_peptide(lfqdata, "~ group_")
#' head(mod$get_coefficients())
build_model_glm_peptide <- function(lfqdata, modelstr) {
  .assert_nested_facade_input(lfqdata, "build_model_glm_peptide")
  lfq_missing <- .prepare_logistf_lfqdata(lfqdata)
  formula <- paste(lfq_missing$get_config()$bin_resp, modelstr)
  build_model_logistf(lfq_missing, formula)
}


#' build_model_logistf
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#' istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 10, with_missing = TRUE,
#'   weight_missing = 0.5, seed = 3)
#' istar$data <- prolfqua::encode_bin_resp(LFQData$new(istar$data, istar$config))
#' istar$config$bin_resp <- "bin_resp"
#' tmp <- LFQData$new(istar$data, istar$config)
#' formula <- paste0(tmp$get_config()$bin_resp , "~ group_")
#' xx2 <- build_model_logistf(tmp, formula)
#'
#' istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 10, with_missing = TRUE,
#'   weight_missing = 0.5, seed = 3)
#' istar$data <- prolfqua::encode_bin_resp(LFQData$new(istar$data, istar$config))
#' istar$config$bin_resp <- "bin_resp"
#' tmp <- LFQData$new(istar$data, istar$config)
#' formula <- paste0(tmp$get_config()$bin_resp , "~ group_")
#' xx <- build_model_logistf(tmp, formula)
#'
#'
#'
#' m <- xx$models$models1$model_df$linear_model[[1]]
#' linfct <- linfct_from_model(m)
#' linfct_all_possible_contrasts(linfct$linfct_factors)
#' x <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_interactions)
#' linfct <- linfct_factors_contrasts(m)
#'
#' m <- xx2$models$models2$model_df$linear_model[[1]]
#' linfct <- linfct_from_model(m)
#' x <- linfct_all_possible_contrasts(linfct$linfct_factors)
#' x <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_interactions)
#' linfct <- linfct_factors_contrasts(m)
#'
#'
#'
build_model_logistf <- function(data, formula) {
  pep <- data
  df <- pep$summarize_hierarchy()
  df2 <- df[df[[ncol(df)]] > 1, ]

  models2 <- NULL
  hkey <- NULL
  if (nrow(df2) > 0) {
    hkey <- tail(pep$get_config()$hierarchy_keys(), n = 1)
    lfq2 <- pep$get_subset(df2)
    formula2 <- paste0(formula, "+", hkey)
    model_strategy2 <- prolfqua::strategy_logistf(formula2)
    models2 <- model_analyse(
      lfq2$data_long(),
      model_strategy2,
      model_name = "logistf_2",
      label = "firth multi-peptide",
      subject_id = lfq2$subject_id()
    )
    models2$strategy <- model_strategy2
  }

  df1 <- df[df[[ncol(df)]] == 1, ]
  models1 <- NULL
  if (nrow(df1) > 0) {
    lfq1 <- pep$get_subset(df1)
    model_strategy1 <- prolfqua::strategy_logistf(formula)
    models1 <- model_analyse(
      lfq1$data_long(),
      model_strategy1,
      model_name = "logistf_1",
      label = "firth single-peptide",
      subject_id = lfq1$subject_id()
    )
    models1$strategy <- model_strategy1
  }
  res <- ModelFirth$new(list(models2 = models2, models1 = models1, hkey = hkey))
  return(res)
}


#' build dataframe with models for testing
#' @family modelling
#' @export
#' @keywords internal
#' @examples
#' modi <- sim_build_models_logistf(model = "interaction", weight_missing = 1)
#' stopifnot(dim(modi$model_df) == c(10,9))
#' mod2 <- sim_build_models_logistf(model = "parallel2", weight_missing = 1)
#' mod2$model_df$linear_model[[1]]
#' mod3 <- sim_build_models_logistf(model = "parallel3", weight_missing = 1)
#' modf <- sim_build_models_logistf(model = "factors", weight_missing = 1)
#'
#' mod3 <- sim_build_models_logistf(model = "parallel3", weight_missing = 1, peptide=TRUE)
#' modf <- sim_build_models_logistf(model = "factors", weight_missing = 1, peptide=TRUE)

sim_build_models_logistf <- function(
  model = c("parallel2", "parallel3", "factors", "interaction"),
  Nprot = 10,
  with_missing = TRUE,
  weight_missing = 1,
  peptide = FALSE
) {
  model <- match.arg(model)
  if (!peptide) {
    if (model != "parallel3") {
      istar <- prolfqua::sim_lfq_data_2factor_config(
        Nprot = Nprot,
        with_missing = with_missing,
        weight_missing = weight_missing
      )
      istar$data <- encode_bin_resp(LFQData$new(istar$data, istar$config))
      istar$config$bin_resp <- "bin_resp"
    } else {
      istar <- prolfqua::sim_lfq_data_protein_config(
        Nprot = Nprot,
        with_missing = with_missing,
        weight_missing = weight_missing
      )
      istar$data <- encode_bin_resp(LFQData$new(istar$data, istar$config))
      istar$config$bin_resp <- "bin_resp"
    }
    istar <- prolfqua::LFQData$new(istar$data, istar$config)
  } else {
    if (model != "parallel3") {
      istar <- prolfqua::sim_lfq_data_2factor_config(
        Nprot = Nprot,
        with_missing = with_missing,
        weight_missing = weight_missing,
        PEPTIDE = TRUE
      )
      istar$data <- encode_bin_resp(LFQData$new(istar$data, istar$config))
      istar$config$bin_resp <- "bin_resp"
    } else {
      istar <- prolfqua::sim_lfq_data_peptide_config(
        Nprot = Nprot,
        with_missing = with_missing,
        weight_missing = weight_missing
      )
      istar$data <- encode_bin_resp(LFQData$new(istar$data, istar$config))
      istar$config$bin_resp <- "bin_resp"
    }
    istar <- prolfqua::LFQData$new(istar$data, istar$config)
  }

  model <- if (model == "factors") {
    "~ Treatment + Background"
  } else if (model == "interaction") {
    "~ Treatment * Background"
  } else if (model == "parallel2") {
    "~ Treatment"
  } else if (model == "parallel3") {
    "~ group_"
  } else {
    NULL
  }
  model_function <- paste0(istar$get_config()$bin_resp, model)
  mod <- build_model_logistf(
    istar,
    model_function
  )
  return(mod)
}


#' Firth's logistic regression strategy (R6 class)
#'
#' Encapsulates everything needed to fit per-protein Firth's bias-reduced
#' logistic regression via \code{\link[logistf]{logistf}} and extract contrasts.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @importFrom logistf logistf
#' @examples
#' strat <- StrategyLogistf$new("bin_resp ~ condition")
#' strat$model_fun(get_formula = TRUE)
StrategyLogistf <- R6::R6Class(
  "StrategyLogistf",
  public = list(
    #' @field formula model formula
    formula = NULL,
    #' @field model_name name of model
    model_name = NULL,
    #' @field report_columns columns to report
    report_columns = NULL,
    #' @field is_mixed always FALSE for logistf
    is_mixed = FALSE,
    #' @field anova_df list with anova function and column names
    anova_df = NULL,

    #' @description Create a new StrategyLogistf
    #' @param modelstr model formula string
    #' @param model_name name of model
    #' @param report_columns columns to report
    #' @param test type of test statistic to use (e.g. "Chisq")
    initialize = function(
      modelstr,
      model_name = "logistf",
      report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value", "moderated.p.value.adjusted"),
      test = "Chisq"
    ) {
      self$formula <- as.formula(modelstr)
      self$model_name <- model_name
      self$report_columns <- report_columns
      self$anova_df <- get_anova_df(test = test)
    },

    #' @description Fit logistf to one protein's data
    #' @param x data.frame for one protein
    #' @param pb optional progress bar
    #' @param get_formula if TRUE, return formula instead of fitting
    model_fun = function(x, pb, get_formula = FALSE) {
      if (get_formula) {
        return(self$formula)
      }
      if (!missing(pb)) {
        pb$tick()
      }
      predictor_vars <- all.vars(update(self$formula, . ~ .))
      DFT <- x |>
        dplyr::group_by(dplyr::across(dplyr::all_of(predictor_vars))) |>
        dplyr::summarize(Freq = dplyr::n(), .groups = "drop")
      tryCatch(logistf::logistf(self$formula, data = DFT, weights = Freq, pl = FALSE), error = .error_handler)
    },

    #' @description Check if model is singular (NA coefficients or df < 2)
    #' @param model fitted model
    isSingular = function(model) {
      if (any(is.na(coefficients(model)))) {
        return(TRUE)
      }
      if (self$df_residual(model) >= 2) {
        return(FALSE)
      }
      TRUE
    },

    #' @description Compute contrasts from fitted model
    #' @param ... passed to \code{\link{compute_contrast}}
    contrast_fun = function(...) compute_contrast(...),

    #' @description Get residual degrees of freedom
    #' @param model fitted model
    df_residual = function(model) {
      model$n - length(model$coefficients)
    },

    #' @description Get residual standard error (always 1 for logistic)
    #' @param model fitted model
    sigma = function(model) 1
  )
)


#' Create Firth's logistic regression strategy
#'
#' Convenience wrapper that creates a \code{\link{StrategyLogistf}} object.
#' @export
#' @rdname strategy
#' @param modelstr model formula
#' @param model_name name of model
#' @param report_columns columns to report
#' @param test type of test statistic to use (e.g. "Chisq")
#' @family modelling
#' @return a \code{\link{StrategyLogistf}} object
#' @examples
#' tmp <- strategy_logistf("bin_resp ~ condition", model_name = "parallel design")
#' tmp$model_fun(get_formula = TRUE)
#'
#' istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 10, with_missing = TRUE,
#'   weight_missing = 0.5, seed = 3)
#' istar$data <- encode_bin_resp(LFQData$new(istar$data, istar$config))
#' istar$config$bin_resp <- "bin_resp"
#' istar <- LFQData$new(istar$data, istar$config)
#' df <- istar$summarize_hierarchy()
#' df2 <- df[df[[ncol(df)]] > 1, ]
#' istar2 <- istar$get_subset(df2)
#' istar2$data_long() |>
#'   dplyr::group_by(protein_Id) |>
#'   tidyr::nest() -> nestProtein
#' modelFunction <- strategy_logistf("bin_resp ~ group_ + peptide_Id",
#'   model_name = "random_example")
#' modelFunction$model_fun(nestProtein$data[[1]])
#' modelFunction$model_fun(nestProtein$data[[4]])
strategy_logistf <- function(
  modelstr,
  model_name = "logistf",
  report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value", "moderated.p.value.adjusted"),
  test = "Chisq"
) {
  StrategyLogistf$new(modelstr, model_name, report_columns, test)
}
