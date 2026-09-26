.prepare_detection_lfqdata <- function(lfqdata) {
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
  lfq_missing <- .prepare_detection_lfqdata(lfqdata)
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
  lfq_missing <- .prepare_detection_lfqdata(lfqdata)
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
#' contrasts <- c(AvsB = "group_A - group_B")
#' m <- xx$models$models1$model_df$linear_model[[1]]
#' linfct_matrix_contrasts(linfct_from_model(m)$linfct_factors, contrasts)
#' m <- xx2$models$models2$model_df$linear_model[[1]]
#' linfct_matrix_contrasts(linfct_from_model(m)$linfct_factors, contrasts)
build_model_logistf <- function(data, formula) {
  df <- data$summarize_hierarchy()
  n_children <- df[[ncol(df)]]
  hkey <- if (any(n_children > 1)) tail(data$get_config()$hierarchy_keys(), n = 1) else NULL
  fit_subset <- function(rows, formula, model_name, label) {
    if (nrow(rows) == 0) {
      return(NULL)
    }
    lfq <- data$get_subset(rows)
    strategy <- strategy_logistf(formula)
    models <- model_analyse(
      lfq$data_long(),
      strategy,
      model_name = model_name,
      label = label,
      subject_id = lfq$subject_id()
    )
    models$strategy <- strategy
    models
  }
  models2 <- fit_subset(df[n_children > 1, ], paste0(formula, "+", hkey), "logistf_2", "firth multi-peptide")
  models1 <- fit_subset(df[n_children == 1, ], formula, "logistf_1", "firth single-peptide")
  ModelFirth$new(list(models2 = models2, models1 = models1, hkey = hkey))
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
#' strat$formula
StrategyLogistf <- R6::R6Class(
  "StrategyLogistf",
  inherit = StrategyBase,
  public = list(
    #' @description Create a new StrategyLogistf
    #' @param modelstr model formula string
    #' @param model_name name of model
    initialize = function(modelstr, model_name = "logistf") {
      super$initialize(modelstr, model_name)
    },

    #' @description Get residual degrees of freedom
    #' @param model fitted model
    df_residual = function(model) model$n - length(model$coefficients),

    #' @description Get residual standard error (always 1 for logistic)
    #' @param model fitted model
    sigma = function(model) 1
  ),
  private = list(
    # Collapse identical predictor rows into frequency weights.
    prepare = function(x) {
      predictor_vars <- all.vars(update(self$formula, . ~ .))
      x |>
        dplyr::group_by(dplyr::across(dplyr::all_of(predictor_vars))) |>
        dplyr::summarize(Freq = dplyr::n(), .groups = "drop")
    },
    fit = function(DFT) logistf::logistf(self$formula, data = DFT, weights = Freq, pl = FALSE)
  )
)


#' Create Firth's logistic regression strategy
#'
#' Convenience wrapper that creates a \code{\link{StrategyLogistf}} object.
#' @export
#' @rdname strategy
#' @family modelling
#' @return a \code{\link{StrategyLogistf}} object
#' @examples
#' tmp <- strategy_logistf("bin_resp ~ condition", model_name = "parallel design")
#' tmp$formula
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
strategy_logistf <- function(modelstr, model_name = "logistf") {
  StrategyLogistf$new(modelstr, model_name)
}
