#' compute contrasts
#'
#'
#' @export
#' @keywords internal
contrasts_linfct_firth <- function(models,
                             linfct,
                             subject_Id = "protein_Id" ,
                             contrastfun = prolfqua::my_contest){
  #computeGroupAverages
  message("contrasts_linfct_firth")
  modelcol <- "linear_model"
  # TODO (goes into calling code)
  # models <- models |> dplyr::filter(.data$exists_lmer == TRUE)

  interaction_models <- vector(mode = "list", length = nrow(models))

  if ("matrix" %in% class(linfct)) {
    pb <- progress::progress_bar$new(total = length(models[[modelcol]]))
    for (i in seq_along(models[[modelcol]])) {
      interaction_models[[i]] <- contrastfun(models[[modelcol]][[i]],
                                             linfct = linfct)

      pb$tick()
    }
    interaction_model_matrix <- models
    interaction_model_matrix$contrast <- interaction_models
  } else if (("list" %in% class(linfct)) && (length(linfct) == nrow(models))) {
    pb <- progress::progress_bar$new(total = length(models[[modelcol]]))

    for (i in seq_along(models[[modelcol]])) {

      interaction_models[[i]] <- contrastfun(models[[modelcol]][[i]],
                                             linfct = linfct[[i]]
      )
      pb$tick()
    }
    interaction_model_matrix <- models
    interaction_model_matrix$contrast <- interaction_models
  } else {
    stop("linct must be either a matrix or a list of length == nrow models")
  }

  #interaction_model_matrix <- models |>
  #  dplyr::mutate("contrast" := purrr::map(!!sym(modelcol) , contrastfun , linfct = linfct ))

  mclass <- function(x){
    class(x)[1]
  }

  interaction_model_matrix <-  interaction_model_matrix |>
    dplyr::mutate(classC = purrr::map_chr(.data$contrast, mclass)) |>
    dplyr::filter(.data$classC != "logical")

  contrasts <- interaction_model_matrix |>
    dplyr::select_at( c(subject_Id, "contrast") ) |>
    tidyr::unnest_legacy()

  # take sigma and df from somewhere else.
  modelInfos <- models |>
    dplyr::select_at(c(subject_Id,
                       "isSingular",
                       "sigma.model" = "sigma",
                       "df.residual.model" = "df.residual" )) |>

    dplyr::distinct()
  contrasts <- dplyr::inner_join(contrasts, modelInfos, by = subject_Id)
  return(ungroup(contrasts))
}


#' build_model_logistf
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#' istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 10, with_missing = TRUE, weight_missing = 0.5, seed = 3)
#' istar$data <- prolfqua::encode_bin_resp(istar$data, istar$config)
#' tmp <- LFQData$new(istar$data, istar$config)
#' formula <- paste0(tmp$config$table$bin_resp , "~ group_")
#' xx <- build_model_logistf(tmp, formula)
#'
#' istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 10, with_missing = TRUE, weight_missing = 0.5, seed = 3)
#' istar$data <- prolfqua::encode_bin_resp(istar$data, istar$config)
#' tmp <- LFQData$new(istar$data, istar$config)
#' formula <- paste0(tmp$config$table$bin_resp , "~ group_")
#' xx <- build_model_logistf(tmp, formula)
#'
#'
#' m <- xx$models$models1$modelDF$linear_model[[1]]
#' linfct <- linfct_from_model(m)
#' linfct_all_possible_contrasts(linfct$linfct_factors)
#' prolfqua::linfct_all_possible_contrasts(linfct$linfct_interactions)
#' linfct <- linfct_factors_contrasts(m)
#'
#' #contrasts_linfct_firth(xx$models$models1$modelDF[1,, drop=FALSE],linfct)
#'
#' contrasts <- c(AvsB = "group_A - group_B", AvsCtrl = "group_A - group_Ctrl")
#' prolfqua:::.linfct(m, contrasts)
#' # Contrasts$debug("get_contrasts")
#'
#' #xx <- Contrasts$new(xx, contrasts)
#' #xx$get_contrast_sides()
#' #xx$get_linfct()
#' #xx$get_contrasts()
#'
build_model_logistf <- function(data,
                                formula){
  pep <- data
  df <- pep$summarize_hierarchy()
  df2 <- df[df[[ncol(df)]] > 1,  ]

  models2 <- NULL
  hkey <- NULL
  if (nrow(df2) > 0) {
    hkey <- tail(pep$config$table$hierarchy_keys(), n = 1)
    lfq2 <- pep$get_subset(df2)
    formula2 <- paste0(formula, "+", hkey)
    model_strategy2 <- prolfqua::strategy_logistf(formula2)
    models2 <- model_analyse(lfq2$data,
                             model_strategy2,
                             modelName = "logistf_2",
                             subject_Id = lfq2$subject_Id())
    models2$strategy = model_strategy2
  }

  df1 <- df[df[[ncol(df)]] == 1,]
  models1 <- NULL
  if (nrow(df1) > 0) {
    lfq1 <- pep$get_subset(df1)
    model_strategy1 <- prolfqua::strategy_logistf(formula)
    models1 <- model_analyse(lfq1$data,
                             model_strategy1,
                             modelName = "logistf_1",
                             subject_Id = lfq1$subject_Id())
    models1$strategy = model_strategy1
  }
  res <- ModelFirth$new(list(models2 = models2, models1 = models1, hkey = hkey))
  return(res)
}


#' build dataframe with models for testing
#' @family modelling
#' @export
#' @keywords internal
#' @examples
#' # debug(sim_build_models_logistf)
#' modi <- sim_build_models_logistf(model = "interaction", weight_missing = 1)
#' stopifnot(dim(modi$modelDF) == c(10,9))
#' mod2 <- sim_build_models_logistf(model = "parallel2", weight_missing = 1)
#' mod2$modelDF$linear_model[[1]]
#' mod3 <- sim_build_models_logistf(model = "parallel3", weight_missing = 1)
#' modf <- sim_build_models_logistf(model = "factors", weight_missing = 1)
#'
sim_build_models_logistf <- function(model = c("parallel2","parallel3","factors", "interaction"),
                                     Nprot = 10,
                                     with_missing = TRUE,
                                     weight_missing = 1) {
  model <- match.arg(model)
  if (model != "parallel3") {
    istar <- prolfqua::sim_lfq_data_2Factor_config(
      Nprot = Nprot,
      with_missing = with_missing,
      weight_missing = weight_missing)
    istar$data <- encode_bin_resp(istar$data, istar$config)
  } else {
    istar <- prolfqua::sim_lfq_data_protein_config(
      Nprot = Nprot,
      with_missing = with_missing,
      weight_missing = weight_missing)
    istar$data <- encode_bin_resp(istar$data, istar$config)
  }
  istar <- prolfqua::LFQData$new(istar$data,istar$config)

  model <- if (model == "factors") {
    "~ Treatment + Background"
  } else if (model == "interaction") {
    "~ Treatment * Background"
  } else if (model == "parallel2") {
    "~ Treatment"
  } else if (model == "parallel3") {
    "~ group_"
  } else {NULL}
  modelFunction <- paste0( istar$config$table$bin_resp, model)
  mod <- build_model_logistf(
    istar,
    modelFunction)
  return(mod)
}


#' Firth's Bias-Reduced Logistic Regression (logistf)
#' @export
#' @rdname strategy
#' @param modelstr model formula
#' @param model_name name of model
#' @param report_columns columns to report
#' @family modelling
#' @examples
#' library(tidyverse)
#' tmp <- strategy_logistf("bin_resp ~ condition", model_name = "parallel design")
#' tmp$model_fun(get_formula = TRUE)
#' tmp$isSingular
#'
#' istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 10, with_missing = TRUE, weight_missing = 0.5, seed = 3)
#' istar$data <- encode_bin_resp(istar$data, istar$config)
#' istar <- LFQData$new(istar$data, istar$config)
#' df <- istar$summarize_hierarchy()
#' df2 <- df[df[[ncol(df)]] > 1,  ]
#' istar2 <- istar$get_subset(df2)
#' istar2$data |>
#' dplyr::group_by(protein_Id) |>
#'  tidyr::nest() -> nestProtein
#' modelFunction <- strategy_logistf("bin_resp ~ group_ + peptide_Id", model_name = "random_example")
#' modelFunction$model_fun(nestProtein$data[[1]])
#' modelFunction$model_fun(nestProtein$data[[4]])
#'
strategy_logistf <- function(
    modelstr,
    model_name = "logistf",
    report_columns = c("statistic",
                       "p.value",
                       "p.value.adjusted",
                       "moderated.p.value",
                       "moderated.p.value.adjusted"),
    test = "Chisq"
) {
  formula <- as.formula(modelstr)
  model_fun <- function(x, pb, get_formula = FALSE){
    if (get_formula) {
      return(formula)
    }
    if (!missing(pb)) {
      pb$tick()
    }
    # tt <- ftable(formula, x)
    # DFT <- as.data.frame(tt)
    predictor_vars <- all.vars(update(formula, . ~ .))
    DFT <- x %>%
      group_by(across(predictor_vars)) %>%
      summarize(Freq = n(), .groups = "drop")

    modelTest <- tryCatch(logistf::logistf( formula ,
                                            data = DFT,
                                            weights = Freq ),
                          error = .ehandler)
    return(modelTest)
  }
  df_residual_logistf = function(m) {
    n <- m$n                          # Number of observations
    p <- length(m$coefficients)       # Number of estimated parameters
    df_residual <- n - p
    return(df_residual)
  }
  isSingular_logistf = function(m)
  {
    anyNA <- any(is.na(coefficients(m)))
    if (anyNA) {
      return(TRUE)
    } else {
      if (df_residual_logistf(m) >= 2) {
        return(FALSE)
      }
      return(TRUE)
    }

  }
  sigma_logistf = function(m){
    return(1)
  }
  res <- list(model_fun = model_fun,
              isSingular = isSingular_logistf,
              contrast_fun = my_contrast_V2,
              model_name = model_name,
              report_columns = report_columns,
              anova_df = get_anova_df(test = test),
              is_mixed = FALSE,
              df_residual = df_residual_logistf,
              sigma = sigma_logistf)
  return(res)
}

#'
#'
#'
.logistf_coeff_matrix <- function(m){
  data <- NULL
  data <- m$model
  interactionColumns <- intersect(attributes(terms(m))$term.labels,colnames(data))
  data <- make_interaction_column(data, interactionColumns, sep = ":")
  coeffs <- coef(m)
  inter <- unique(data$interaction)
  mm <- matrix(0, nrow = length(inter), ncol = length(coeffs))
  rownames(mm) <- inter
  colnames(mm) <- names(coeffs)
  mm[,1] <- 1
  coefi <- coeffs[-1]
  for (i in seq_along(coefi)) {
    # the grep is needed to extract coefficients of interaction terms belonging to a factor
    # I am using wor boundaries "\\b" to allow for factor levels that are substrings.
    positionIDX <- grep(paste0("\\b",names(coefi)[i],"\\b"), inter)
    mm[positionIDX, i + 1 ] <- 1
  }
  return(list(mm = mm, coeffs = coeffs))
}

