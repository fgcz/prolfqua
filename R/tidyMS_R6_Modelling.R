#' build dataframe with models for testing
#' @family modelling
#' @export
#' @keywords internal
#' @examples
#' mod <- build_models(model = "interaction", weight_missing = 1)
#' stopifnot(dim(mod$modelDF) == c(10,9))
#'
build_models <- function(model = c("factors", "interaction"), Nprot = 10, with_missing = TRUE, weight_missing = 1) {
  model <- match.arg(model)
  model <- if (model == "factors") {
    "~ Treatment + Background"
  } else {
    "~ Treatment * Background"
  }
  istar <- prolfqua::sim_lfq_data_2Factor_config(Nprot = Nprot, with_missing = with_missing, weight_missing = weight_missing)
  istar <- prolfqua::LFQData$new(istar$data,istar$config)
  modelFunction <- strategy_lm(paste0(istar$response(), model))
  mod <- build_model(
    istar,
    modelFunction)
  return(mod)
}

#' make interaction model for examples
#' @family modelling
#' @export
#' @keywords internal
#' @examples
#' m <- make_model()
make_model <- function(model = c("factors", "interaction")){
  mod <- build_models(model = model, Nprot = 1, with_missing = FALSE)
  return(mod$modelDF$linear_model[[1]])
}


# Creating models from configuration ----

.ehandler = function(e){
  warning("WARN :", e)
  # return string here
  as.character(e)
}

#' Create custom lmer model
#' @rdname strategy
#' @param modelstr model formula
#' @param model_name name of model
#' @param report_columns columns to report
#' @export
#' @family modelling
#' @examples
#'
#' istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 10, with_missing = FALSE)
#' istar <- prolfqua::LFQData$new(istar$data,istar$config)
#' istar$data <- istar$data |> dplyr::group_by(protein_Id) |>
#' dplyr::mutate(abundanceC = abundance - mean(abundance)) |> dplyr::ungroup()
#' istar$factors()
#' modelFunction <- strategy_lmer("abundanceC ~ group_ + (1|peptide_Id) ", model_name = "random_example")
#' mod <- build_model(
#'  istar,
#'  modelFunction)
#' sum(mod$modelDF$exists_lmer)
#' sum(mod$modelDF$isSingular, na.rm=TRUE)
#'
strategy_lmer <- function(modelstr,
                          model_name = "Model",
                          report_columns = c("statistic",
                                             "p.value",
                                             "p.value.adjusted",
                                             "moderated.p.value",
                                             "moderated.p.value.adjusted")
) {
  formula <- as.formula(modelstr)
  model_fun <- function(x, pb , get_formula = FALSE){
    if (get_formula) {
      return(formula)
    }
    if (!missing(pb)) {
      pb$tick()
    }
    modelTest <- tryCatch(lmerTest::lmer(formula , data = x ),
                          error = .ehandler)
    return(modelTest)
  }
  res <- list(model_fun = model_fun,
              isSingular = lme4::isSingular,
              contrast_fun = my_contest,
              model_name = model_name,
              report_columns = report_columns,
              anova_df = get_anova_df(test = "F"),
              is_mixed = TRUE)
  return(res)
}

#' Create linear model
#'
#' The strategy contains functions to fit the model but also compute the contrasts etc.
#' @rdname strategy
#' @export
#' @param modelstr model formula
#' @param model_name name of model
#' @param report_columns columns to report
#' @family modelling
#' @return list with model function, contrast computation function etc.
#' @examples
#'
#'
#' tmp <- strategy_lm("Intensity ~ condition", model_name = "parallel design")
#' tmp$model_fun(get_formula = TRUE)
#' tmp$isSingular
strategy_lm <- function(modelstr,
                        model_name = "Model",
                        report_columns = c("statistic",
                                           "p.value",
                                           "p.value.adjusted",
                                           "moderated.p.value",
                                           "moderated.p.value.adjusted")
) {
  formula <- as.formula(modelstr)
  model_fun <- function(x, pb, get_formula = FALSE){
    if (get_formula) {
      return(formula)
    }
    if (!missing(pb)) {
      pb$tick()
    }
    modelTest <- tryCatch(lm( formula , data = x ),
                          error = .ehandler)
    return(modelTest)
  }
  res <- list(model_fun = model_fun,
              isSingular = isSingular_lm,
              contrast_fun = my_contrast_V2,
              model_name = model_name,
              report_columns = report_columns,
              anova_df = get_anova_df(test = "F"),
              is_mixed = FALSE)
  return(res)
}


#' Create robust linear regression model
#'
#' The strategy contains functions to fit the model but also compute the contrasts etc.
#'
#' @rdname strategy
#' @export
#' @param modelstr model formula
#' @param model_name name of model
#' @param report_columns columns to report
#' @family modelling
#' @return list with model function, contrast computation function etc.
#' @examples
#' tmp <- strategy_rlm("Intensity ~ condition", model_name = "parallel design")
#' tmp$model_fun(get_formula = TRUE)
#' tmp$isSingular
strategy_rlm <- function(modelstr,
                         model_name = "Model",
                         report_columns = c("statistic",
                                            "p.value",
                                            "p.value.adjusted",
                                            "moderated.p.value",
                                            "moderated.p.value.adjusted")
) {
  formula <- as.formula(modelstr)
  model_fun <- function(x, pb, get_formula = FALSE){
    if (get_formula) {
      return(formula)
    }
    if (!missing(pb)) {
      pb$tick()
    }
    modelTest <- tryCatch(MASS::rlm( formula , data = x , method = "M" ),
                          error = .ehandler)
    return(modelTest)
  }
  res <- list(model_fun = model_fun,
              isSingular = isSingular_lm,
              contrast_fun = my_contrast_V2,
              model_name = model_name,
              report_columns = report_columns,
              anova_df = get_anova_df(test = "F"),
              is_mixed = FALSE)
  return(res)
}

#' Create quasibinomial glm model
#' @export
#' @rdname strategy
#' @param modelstr model formula
#' @param model_name name of model
#' @param family either binomial or quasibinomial
#' @param multiplier for tuning default is 1.
#' @param report_columns columns to report
#' @family modelling
#' @examples
#' tmp <- strategy_glm("Intensity ~ condition", model_name = "parallel design")
#' tmp$model_fun(get_formula = TRUE)
#' tmp$isSingular
strategy_glm <- function(modelstr,
                         model_name = "Model",
                         test = "Chisq",
                         family = stats::binomial,
                         multiplier = 1,
                         offset = 1,
                         report_columns = c("statistic",
                                            "p.value",
                                            "p.value.adjusted",
                                            "moderated.p.value",
                                            "moderated.p.value.adjusted")
) {
  formula <- as.formula(modelstr)
  model_fun <- function(x, pb, get_formula = FALSE){
    if (get_formula) {
      return(formula)
    }
    if (!missing(pb)) {
      pb$tick()
    }
    # to avoid perfect separation (hack)
    tt <- ftable(formula, x)
    tt <- tt * multiplier + offset
    DFT <- as.data.frame(tt)
    modelTest <- tryCatch(glm( formula ,
                               data = DFT ,
                               weights = Freq,
                               family = family),
                          error = .ehandler)
    return(modelTest)
  }
  res <- list(model_fun = model_fun,
              isSingular = isSingular_lm,
              contrast_fun = my_contrast_V2,
              model_name = model_name,
              report_columns = report_columns,
              anova_df = get_anova_df(test = test),
              is_mixed = FALSE)
  return(res)
}

#' anova returning dataframe
#' @keywords internal
#' @family modelling
#' @export
#' @examples
#' x <- get_anova_df(test = "F")
#' x <- get_anova_df(test = "Chisq")
get_anova_df <- function(test = "F"){
  res <- function(x){
    x <- anova(x, test = test)
    colnames(x) <- make.names(colnames(x))
    x <- data.frame(factor = rownames(x), x)
    return(x)
  }
  return(list(fun = res,
              col_pval = paste0("Pr..",substr(test,1,3),"."),
              col_fdr = paste0("FDR.Pr..",substr(test,1,3),".")))

}



.likelihood_ratio_test <- function(modelNO, model) {
  res <- tryCatch(  anova(modelNO,model), error = function(x) NULL)
  if (!is.null(res)) {
    res <- suppressWarnings(broom::tidy(res))[2,"p.value"]
    return(as.numeric(res))

  }
  else{
    return(NA)
  }
}


# Fit the models to data ----


#' check if lm model is singular
#' @keywords internal
#' @family modelling
#' @export
#'
isSingular_lm <- function(m){
  anyNA <- any(is.na(coefficients(m)))
  if (anyNA) {
    return(TRUE)
  } else {
    if (df.residual(m) >= 2) {
      return(FALSE)
    }
    return(TRUE)
  }
}

#' retrieve complete model.
#' @keywords internal
#' @family modelling
#' @export
#'
get_complete_model_fit <- function(modelProteinF){
  modelProteinF <- modelProteinF |> dplyr::filter(.data$exists_lmer == TRUE)
  modelProteinF <- modelProteinF |> dplyr::filter(.data$nrcoeff_not_NA == max(.data$nrcoeff_not_NA)) |>
    dplyr::arrange(dplyr::desc(.data$nrcoeff_not_NA))
  modelProteinF <- modelProteinF |> dplyr::filter(df.residual > 1)
  return(modelProteinF)
}

#' analyses lmer4 and lm models created using help function `strategy_lm` or `strategy_lmer`
#'
#' used in project p2901
#'
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#'
#'
#' ionstar <- prolfqua_data('data_ionstar')$normalized()
#' ionstar$config <- old2new(ionstar$config)
#'
#' ionstar$data <- ionstar$data |> dplyr::filter(protein_Id %in% sample(protein_Id,10))
#' prolfqua::table_factors(ionstar$data, ionstar$config)
#' formula_randomPeptide <-
#'   strategy_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id)")
#' mr <- model_analyse( ionstar$data,
#'  formula_randomPeptide,
#'  subject_Id = ionstar$config$table$hierarchy_keys_depth())
#' get_complete_model_fit(mr$modelProtein)
model_analyse <- function(pepIntensity,
                          model_strategy,
                          subject_Id = "protein_Id",
                          modelName = "Model")
{
  pepIntensity |>
    dplyr::group_by(!!!syms(subject_Id)) |>
    tidyr::nest() -> nestProtein

  lmermodel <- "linear_model"

  pb <- progress::progress_bar$new(total = nrow(nestProtein))
  modelProtein <- nestProtein |>
    dplyr::mutate(!!lmermodel := purrr::map(data, model_strategy$model_fun, pb = pb))

  modelProtein <- modelProtein |>
    dplyr::mutate(!!"exists_lmer" := purrr::map_lgl(!!sym(lmermodel), function(x){!is.character(x)}))

  modelProteinF <- modelProtein |>
    dplyr::filter( !!sym("exists_lmer") == TRUE)
  modelProteinF <- modelProteinF |>
    dplyr::mutate(!!"isSingular" := purrr::map_lgl(!!sym(lmermodel), model_strategy$isSingular ))
  modelProteinF <- modelProteinF |>
    dplyr::mutate(!!"df.residual" := purrr::map_dbl(!!sym(lmermodel), df.residual ))
  modelProteinF <- modelProteinF |>
    dplyr::mutate(!!"sigma" := purrr::map_dbl( !!sym(lmermodel) , sigma))

  nrcoeff <- function(x) {
    cc <- coefficients(x)
    if (class(cc) == "numeric") {
      return(length(cc))
    }else{
      return(ncol(cc[[1]]))
    }
  }

  nrcoeff_not_NA <- function(x){
    cc <- coefficients(x)
    if (class(cc) == "numeric") {
      return(sum(!is.na(cc)))
    }else{
      return(ncol(cc[[1]]))
    }
  }

  modelProteinF <- modelProteinF |> dplyr::mutate(nrcoef = purrr::map_int(!!sym(lmermodel), nrcoeff))
  modelProteinF <- modelProteinF |> dplyr::mutate(nrcoeff_not_NA = purrr::map_int(!!sym(lmermodel), nrcoeff_not_NA))

  #return(list(modelProtein = modelProtein, modelProteinF = modelProteinF))
  modelProteinF <- modelProteinF |>
    dplyr::select_at(c(subject_Id,"isSingular", "df.residual","sigma" ,"nrcoef", "nrcoeff_not_NA") )
  modelProtein <- dplyr::left_join(modelProtein, modelProteinF)

  return(list(modelProtein = modelProtein,
              modelName = modelName
  ))
}



# visualize lmer modelling results ----

#' Plot prdictions
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#' m <- prolfqua_data('data_interactionModel_p1807')
#' plot_lmer_peptide_predictions(m)
plot_lmer_peptide_predictions <- function(m){
  data <- m@frame
  data$prediction <- predict(m)
  interactionColumns <- intersect(attributes(terms(m))$term.labels,colnames(data))
  data <- make_interaction_column(data, interactionColumns, sep = ":")
  gg <- ggplot(data, aes(x = .data$interaction ,
                         y = .data$transformedIntensity)) + geom_point()
  gg <- gg + geom_point(aes(x = .data$interaction,
                            y = .data$prediction), color = 2) + facet_wrap(~peptide_Id)
  gg <- gg + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  return(gg)
}


#' plot peptide intensities per interaction with random effects removed
#'
#' @param m model
#' @param legend.position none
#' @family modelling
#' @export
#' @keywords internal
#' @examples
#'
#'
#' m <- prolfqua_data('data_basicModel_p1807')
#' plot_lmer_peptide_noRandom(m)
#'
#' m <- prolfqua_data('data_interactionModel_p1807')
#' plot_lmer_peptide_noRandom(m)
plot_lmer_peptide_noRandom <- function(m,legend.position = "none"){
  data <- m@frame
  ran <- lme4::ranef(m)[[1]]
  randeffect <- base::setdiff(all.vars( terms(formula(m)) ) , all.vars(terms(m)))
  ran <- tibble::as_tibble(ran,rownames = randeffect)
  colnames(ran) <- gsub("[()]","",colnames(ran))
  ran <- dplyr::inner_join(data, ran, by = randeffect)

  ran <- ran |> dplyr::mutate(int_randcorrected  = .data$transformedIntensity  - .data$Intercept)
  interactionColumns <- intersect(attributes(terms(m))$term.labels,colnames(data))
  ran <- make_interaction_column(ran,interactionColumns, sep = ":" )

  meanx <- function(x){mean(x,na.rm = TRUE)}
  gg <- ggplot(ran,aes(x = .data$interaction,
                       y = .data$int_randcorrected,
                       color = .data$peptide_Id)) +
    geom_point(position = position_jitterdodge())
  gg <- gg + stat_summary(fun = meanx, colour = "black", geom = "point",
                          shape = 12, size = 3,show.legend = FALSE)
  gg <- gg + theme(axis.text.x = element_text(angle = -90, hjust = 0),
                   legend.position = legend.position)
  gg <- gg + geom_boxplot(alpha = 0.1)

  return(gg)
}



#' Add predicted values for each interaction
#' @export
#' @keywords internal
#' @family modelling
#' @examples
#' m <- prolfqua_data('data_interactionModel_p1807')
#' plot_lmer_predicted_interactions(plot_lmer_model_and_data(m,"dumm"),m)
plot_lmer_predicted_interactions <- function(gg, m){
  cm <- .lmer4_coeff_matrix(m)
  xstart_end <- data.frame(xstart = rownames(cm$mm), xend = rownames(cm$mm))
  ystart_end <- data.frame(xend = rownames(cm$mm), ystart = rep(0, nrow(cm$mm)),
                           yend = cm$mm %*% cm$coeffs)
  segments <- dplyr::inner_join(xstart_end, ystart_end, by = "xend")
  gg <- gg + geom_segment(aes(x = .data$xstart,
                              y = .data$ystart,
                              xend = .data$xend,
                              yend = .data$yend),
                          data = segments, color = "blue", arrow = arrow())
  return(gg)
}

#' Make model plot with title - protein Name.
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#'
#' m <- prolfqua_data('data_interactionModel_p1807')
#' plot_lmer_model_and_data(m,"dumm")
#'
plot_lmer_model_and_data <- function(m, proteinID, legend.position = "none"){
  gg <- plot_lmer_peptide_noRandom(m,legend.position = legend.position)
  gg <- plot_lmer_predicted_interactions(gg, m)
  gg <- gg + ggtitle(proteinID)
  gg
}





# Generate linear functions -----

# get matrix of indicator coefficients for each interaction
.lmer4_coeff_matrix <- function(m){
  data <- NULL
  if ("lm" %in% class(m)) {
    data <- m$model
  }else{
    # for "lmerModLmerTest"
    data <- m@frame
  }
  interactionColumns <- intersect(attributes(terms(m))$term.labels,colnames(data))
  data <- make_interaction_column(data, interactionColumns, sep = ":")

  if("rlm" %in% class(m)){
    coeffs <- coefficients(summary(m))[,'Value']
  } else {
    coeffs <- coefficients(summary(m))[,'Estimate']
  }
  inter <- unique(data$interaction)
  mm <- matrix(0, nrow = length(inter), ncol = length(coeffs))
  rownames(mm) <- inter
  colnames(mm) <- names(coeffs)
  mm[,1] <- 1
  coefi <- coeffs[-1]
  for (i in seq_along(coefi)) {
    positionIDX <- grep(names(coefi)[i], inter)
    mm[positionIDX, i + 1 ] <- 1
  }
  return(list(mm = mm, coeffs = coeffs))
}


.get_match_idx <- function(mm, factor_level){
  ddd <- names_to_matrix(rownames(mm), split = ":")
  xd <- apply(ddd, 2, function(x, factor_level){x %in% factor_level}, factor_level)
  idx <- which(apply(xd,1, sum) > 0)
  return(idx)
}

.coeff_weights_factor_levels <- function(mm){
  getCoeffs <- function(factor_level, mm){
    idx <- .get_match_idx(mm, factor_level)
    x <- as.list(apply(mm[idx,, drop = FALSE],2,mean) )
    x <- tibble::as_tibble(x)
    tibble::add_column(x, "factor_level" = factor_level,.before = 1)
  }
  factor_levels <- unique(unlist(stringr::str_split(rownames(mm), ":")))
  xx <- purrr::map_df(factor_levels, getCoeffs, mm)
  return(xx)
}

#' get linfct from model
#' @param m linear model
#' @family modelling
#' @export
#' @keywords internal
#' @examples
#'
#' m <- prolfqua_data('data_basicModel_p1807')
#' # debug(linfct_from_model)
#' linfct <- linfct_from_model(m)
#'
#' linfct$linfct_factors
#' linfct$linfct_interactions
#'
#' m <- prolfqua_data('data_interactionModel_p1807')
#' # debug(.coeff_weights_factor_levels)
#' linfct <- linfct_from_model(m)
#'
#' all.equal(linfct$linfct_factors["CelltypeCMP/MEP",] ,
#'  apply(linfct$linfct_interactions[grep("CelltypeCMP/MEP", rownames(linfct$linfct_interactions)),],2, mean))
#' linfct$linfct_interactions
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
linfct_from_model <- function(m, as_list = TRUE){

  cm <- .lmer4_coeff_matrix(m)
  cm_mm <- cm$mm[order(rownames(cm$mm)),]

  l_factors <- .coeff_weights_factor_levels(cm_mm)
  linfct_factors <- l_factors |>
    dplyr::select(-.data$factor_level) |>
    data.matrix()

  rownames(linfct_factors) <- l_factors$factor_level
  linfct_factors <- linfct_factors[order(rownames(linfct_factors)),]
  res <- list(linfct_factors = linfct_factors , linfct_interactions = cm_mm)

  if (as_list) {
    return(res)
  }else{
    do.call( rbind, res)
  }
}




#' linfct_matrix_contrasts
#' @export
#' @param linfct linear functions as created by linfct_from_model
#' @param contrast contrasts to determine linear functions for
#' @param p.message print messages default FALSE
#'
#' @family modelling
#' @keywords internal
#' @examples
#'
#' m <- make_model( "factors")
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
#' x<- linfct_matrix_contrasts(linfct, Contr )
#' stopifnot(sum(x["interactXZ",]) ==0 )
#' stopifnot(sum(x["interactAB",]) ==0 )
#'
#' m <- make_model( "interaction")
#' linfct <- linfct_from_model(m, as_list = FALSE)
#' x<- linfct_matrix_contrasts(linfct, Contr )
#' stopifnot(sum(x["interactXZ",]) ==1 )
#' stopifnot(sum(x["interactAB",]) ==1 )
linfct_matrix_contrasts <- function(linfct , contrasts, p.message = FALSE){
  linfct <- t(linfct)
  df <- tibble::as_tibble(linfct, rownames = "interaction")
  make_contrasts <- function(data,
                             contrasts)
  {
    cnams <- base::setdiff(colnames(data),"interaction")
    for (i in seq_along(contrasts)) {
      if (p.message) {message(names(contrasts)[i], "=", contrasts[i],"\n")}
      expr_string <- as.character(rlang::parse_expr(contrasts[i]))
      tryCatch({
        data <- dplyr::mutate(data, !!names(contrasts)[i] := !!rlang::parse_expr(contrasts[i]))
      }, error = function(e) {
        warning("Error:", e$message, "\n")
        # Handle the error, e.g., by skipping the current iteration, logging the error, etc.
      })
    }
    res <- data |> dplyr::select(-one_of(cnams))
    return(res)
  }

  res <- make_contrasts(df, contrasts )
  res <- tibble::column_to_rownames(res,"interaction")
  res <- t(res)
  return(res)
}


#' create all possible contrasts
#' @export
#' @keywords internal
#' @family modelling
#' @examples
#' m <- make_model( "interaction")
#' linfct <- linfct_from_model(m)
#' xl <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_factors)
#' xx <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_interactions)
#'
linfct_all_possible_contrasts <- function(lin_int ){
  combs <- combn(nrow(lin_int),2)
  names <- rownames(lin_int)
  newnames <- rep("", ncol(combs))
  new_lin_fct <- matrix(NA,  nrow = ncol(combs), ncol = ncol(lin_int))
  for (i in seq_len(ncol(combs))) {
    newnames[i] <- paste(names[combs[,i]], collapse = " - ")
    new_lin_fct[i,] <- lin_int[combs[1,i],] - lin_int[combs[2,i],]
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
#' m <- make_model( "interaction")
#' xl <- linfct_factors_contrasts(m)
#' m <- lm(Petal.Width ~ Species, data = iris)
#' linfct_factors_contrasts(m)
linfct_factors_contrasts <- function(m){
  ffac <- attributes(terms(m))$term.labels
  ffac <- ffac[!grepl(":",ffac)] # remove interactions
  linfct_factors <- linfct_from_model(m)$linfct_factors

  factorDepths <- rownames(linfct_factors)
  res <- vector(length(ffac), mode = "list")
  for (i in seq_along(ffac)) {
    fac <- ffac[i]
    idx <- grep(fac, factorDepths)
    linfct_m <- linfct_factors[idx,]
    res[[i]] <- linfct_all_possible_contrasts(linfct_m)
  }
  res <- do.call(rbind, res)
  return(res)
}

# Computing contrasts helpers -----

#' apply multcomp::glht method to linfct
#'
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#'
#' mb <- make_model( "interaction")
#' linfct <- linfct_from_model(mb)
#' names(linfct)
#' my_glht(mb, linfct$linfct_factors)
#'
#' m <-  make_model( "factors")
#' linfct <- linfct_from_model(m)$linfct_factors
#' my_glht(m, linfct)
#'
my_glht <- function(model, linfct , sep = TRUE ) {
  if (!class(model) == "lm") # fixes issue of mutlcomp not working on factors of class character
  {
    warning("USE ONLY WITH LM models ", class(model))
    if (length(lme4::fixef(model)) != ncol(linfct)) {
      return(NA) # catch rank defficient
    }
  }else{
    if (isSingular_lm(model)) {
      return(NA)
    }
    model$model <- as.data.frame(unclass(model$model))
  }
  if (sep) {
    res <- list()
    for (i in seq_len(nrow(linfct))) {
      x <- multcomp::glht(model, linfct = linfct[i,,drop = FALSE])
      RHS <- broom::tidy(confint(x)) |> dplyr::select(-.data$estimate)

      RHS$df <- x$df
      RHS$sigma <- sigma(model)

      x <- dplyr::inner_join(broom::tidy(summary(x)), RHS, by = c("contrast")) # |> dplyr::select(-contrast)
      res[[i]] <- x
    }
    res <- dplyr::bind_rows(res)
    return(res)
  }else{
    x <- multcomp::glht(model, linfct = linfct)
    RHS <- broom::tidy(confint(x)) |> dplyr::select(-.data$estimate)
    RHS$df <- x$df
    RHS$sigma <- sigma(model)
    res <- dplyr::inner_join(broom::tidy(summary(x)), RHS, by = c("contrast")) |>
      dplyr::select(-.data$rhs)
    return(res)
  }
}

#' compute contrasts for full models
#' @param m linear model generated using lm
#' @param linfct linear function
#' @param coef use default
#' @param use default
#' @param confint which confidence interval to determine
#'
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#'
#' m <-  make_model( "factors")
#' linfct <- linfct_from_model(m)$linfct_factors
#' my_glht(m, linfct)
#' my_contrast(m, linfct, confint = 0.95)
#' my_contrast(m, linfct, confint = 0.99)
#'
my_contrast <- function(m,
                        linfct,
                        coef = coefficients(m),
                        Sigma.hat = vcov(m),
                        confint = 0.95){

  df <- df.residual(m)
  sigma <- sigma(m)

  estimate <- linfct %*% t(t(coef))

  if (df > 0) {
    std.error <- sqrt(diag(linfct %*% Sigma.hat %*% t(linfct)))
    statistic <- estimate / std.error

    #p.value <- pt(-abs(statistic), df = df) * 2

    p.value <- pt(abs(statistic), df = df, lower.tail = FALSE) * 2
    prqt <- -qt((1 - confint)/2, df = df)
    conf.low <- estimate  - prqt * std.error
    conf.high <- estimate + prqt * std.error

  }else{
    std.error <- NA
    statistic <- NA
    p.value <- NA
    conf.low <- NA
    conf.high <- NA
  }

  res <- data.frame(lhs = rownames(linfct),
                    sigma = sigma,
                    df = df,
                    estimate = estimate,
                    std.error = std.error,
                    statistic = statistic ,
                    p.value = p.value,
                    conf.low = conf.low,
                    conf.high = conf.high,
                    stringsAsFactors = FALSE)
  return(res)
}

#' handles incomplete models by setting coefficients to 0
#' @param m linear model generated using lm
#' @param linfct linear function
#' @param confint confidence interval default 0.95
#'
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#' m <- make_model( "factors")
#' linfct <- linfct_from_model(m)$linfct_factors
#' my_contrast_V1(m, linfct, confint = 0.95)
#' my_contrast_V1(m, linfct, confint = 0.99)
my_contrast_V1 <- function(incomplete, linfct, confint = 0.95){
  Sigma.hat <- vcov(incomplete)
  Sigma.hat[is.na(Sigma.hat)] <- 0
  coef <- coefficients(incomplete)
  coef[is.na(coef)] <- 0
  res <- my_contrast(incomplete, linfct,
                     coef = coef,
                     Sigma.hat = Sigma.hat,
                     confint = confint)
  return(res)
}

#' handles incomplete models
#'
#' only keeps non NA coefficients.
#'
#' @param m linear model generated using lm
#' @param linfct linear function
#' @param confint confidence interval default 0.95
#'
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#' m <- make_model( "factors")
#' linfct <- linfct_from_model(m)$linfct_factors
#' my_contrast_V2(m, linfct, confint = 0.95)
#' my_contrast_V2(m, linfct, confint = 0.99)
#'
my_contrast_V2 <- function(m, linfct,confint = 0.95){
  Sigma.hat <- vcov(m)

  coef <- na.omit(coefficients(m))

  res <- vector(nrow(linfct), mode = "list")
  for (i in seq_len(nrow(linfct))) {
    linfct_v <- linfct[i,,drop = FALSE]
    idx <- which(linfct_v != 0)
    nam <- colnames(linfct_v)[idx]

    if (all(nam %in% names(coef))) {
      linfct_v_red <- linfct_v[, nam, drop = FALSE]
      Sigma.hat_red <- Sigma.hat[nam,nam,drop = FALSE]
      coef_red <- coef[nam]
      stopifnot(all.equal(colnames(linfct_v_red),colnames(Sigma.hat_red)))
      stopifnot(all.equal(colnames(linfct_v_red),names(coef_red)))
      res[[i]] <- my_contrast(m,linfct_v_red,
                              coef = coef_red,
                              Sigma.hat = Sigma.hat_red, confint = confint)
    }else{
      res[[i]] <-  data.frame(lhs = rownames(linfct_v),
                              sigma = sigma(m),
                              df = df.residual(m),
                              estimate = NA,
                              std.error = NA,
                              statistic = NA ,
                              p.value = NA,
                              conf.low = NA,
                              conf.high = NA,
                              stringsAsFactors = FALSE)
    }
  }
  res <- dplyr::bind_rows(res)
  return(res)
}

#' applies contrast computation using lmerTest::contest function
#' @param m mixed effects model
#' @param linfct linear function
#' @param ddf method to determine denominator degrees of freedom
#' @family modelling
#' @export
#' @keywords internal
#' @examples
#'
#' mb <- prolfqua_data('data_basicModel_p1807')
#' summary(mb)
#'
#' linfct <- linfct_from_model(mb)
#' names(linfct)
#' my_contest(mb, linfct$linfct_factors)
#' my_contest(mb, linfct$linfct_interactions)
#'
#' # my_glht(mb, linfct$linfct_factors)
#' # my_glht(mb, linfct$linfct_interactions)
#' lmerTest::contest(mb, c( 0 ,1 , 0 , 0),joint = FALSE)
#' summary(mb)
#'
#'
#' #library(pbkrtest)
#' #(fm1 <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
#' #class(fm1)
#' #pbkrtest::get_ddf_Lb.lmerMod(fm1)
#'
my_contest <- function(model, linfct, ddf = c("Satterthwaite", "Kenward-Roger")){
  ddf <- match.arg(ddf)
  if (length(lme4::fixef(model)) != ncol(linfct) ) {
    warning("Model is rank deficient!")
    return(NA) # catch rank defficient
  }else{
    res <- lmerTest::contest(model,
                             linfct,
                             joint = FALSE,
                             confint = TRUE,
                             ddf = ddf)
  }
  res <- tibble::as_tibble(res, rownames = "lhs")
  res$sigma <- sigma(model)
  res <- res |> dplyr::rename(estimate = "Estimate",
                              std.error = "Std. Error",
                              statistic = "t value",
                              p.value = "Pr(>|t|)",
                              conf.low = "lower",
                              conf.high = "upper")
  return(res)
}

# computing contrast ----
#' pivot model contrasts matrix to wide format produced by `contrasts_linfct` and ...
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#'
#' # this function is used by the contrast classes to implement the to wide method
#'
pivot_model_contrasts_2_Wide <- function(modelWithInteractionsContrasts,
                                         subject_Id = "protein_Id",
                                         columns = c("estimate", "p.value","p.value.adjusted"),
                                         contrast = "lhs"){

  m_spread <- function(longContrasts, subject_Id, column , contrast){
    res <- longContrasts |>
      dplyr::select(all_of(c(subject_Id, contrast,  column)))
    res <- res |> dplyr::mutate(!!contrast := paste0(column, ".", !!sym(contrast)))
    res <- res |> tidyr::pivot_wider(names_from = contrast, values_from = column)
    return(res)
  }
  res <- list()
  for (column in columns) {
    res[[column]] <- m_spread(modelWithInteractionsContrasts, subject_Id,column, contrast)
  }
  res <- res |> reduce(left_join, by = c(subject_Id))
  return(res)
}
#' compute group averages
#'
#' used in p2621, p2109
#'
#' @export
#' @keywords internal
#' @examples
#' modelSummary_A <- build_models()
#' m <- get_complete_model_fit(modelSummary_A$modelDF)
#'
#' factor_contrasts <- linfct_factors_contrasts( m$linear_model[[1]])
#' factor_levelContrasts <- contrasts_linfct( m,
#'         factor_contrasts,
#'         subject_Id = "protein_Id",
#'         contrastfun = prolfqua::my_contrast_V2)
#'
#'
contrasts_linfct <- function(models,
                             linfct,
                             subject_Id = "protein_Id" ,
                             contrastfun = prolfqua::my_contest){
  #computeGroupAverages
  message("computing contrasts.")
  modelcol <- "linear_model"
  # TODO (goes into calling code)
  # models <- models |> dplyr::filter(.data$exists_lmer == TRUE)

  interaction_models <- vector(mode = "list", length = nrow(models))

  if ("matrix" %in% class(linfct)) {
    for (i in seq_along(models[[modelcol]])) {
      interaction_models[[i]] <- contrastfun(models[[modelcol]][[i]], linfct = linfct)
    }
    interaction_model_matrix <- models
    interaction_model_matrix$contrast <- interaction_models
  } else if (("list" %in% class(linfct)) && (length(linfct) == nrow(models))) {
    for (i in seq_along(models[[modelcol]])) {
      interaction_models[[i]] <- contrastfun(models[[modelcol]][[i]], linfct = linfct[[i]])
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



# LIMMA ----

#' Moderate p-values - limma approach
#' @export
#' @family modelling
#' @keywords internal
#'
moderated_p_limma <- function(mm, df = "df", estimate = "diff", robust = FALSE, confint = 0.95){
  sv <- prolfqua::squeezeVarRob(mm$sigma^2, df = mm[[df]],robust = robust)

  # pior degrees of freedom are Inf
  if (all(is.infinite(sv$df.prior))) {
    sv$df.prior <- mean(mm[[df]]) * nrow(mm)/10
  }

  sv <- tibble::as_tibble(sv)
  sv <- setNames(sv, paste0('moderated.', names(sv)))
  mm <- dplyr::bind_cols(mm, sv)
  mm <- mm |> dplyr::mutate(moderated.statistic  =  .data$statistic * .data$sigma /  sqrt(.data$moderated.var.post))
  mm <- mm |> dplyr::mutate(moderated.df.total = !!sym(df) + .data$moderated.df.prior)
  mm <- mm |> dplyr::mutate(moderated.p.value = 2*pt( abs(.data$moderated.statistic),
                                                      df = .data$moderated.df.total, lower.tail = FALSE) )

  prqt <- -qt((1 - confint)/2, df = mm$moderated.df.total)
  mm$moderated.conf.low <- mm[[estimate]]  - prqt * sqrt(mm$moderated.var.post)
  mm$moderated.conf.high <- mm[[estimate]] + prqt * sqrt(mm$moderated.var.post)
  mm <-  dplyr::ungroup(mm)


  return(mm)
}

#' Moderate p-value for long table
#' @param mm result of \code{\link{contrasts_linfct}}
#' @param group_by_col colnames with contrast description - default 'lhs'
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#'
#' mod <- build_models()
#' m <- get_complete_model_fit(mod$modelDF)
#' factor_contrasts <- linfct_factors_contrasts(m$linear_model[[1]])
#' factor_levelContrasts <- contrasts_linfct(
#'   mod$modelDF,
#'   factor_contrasts,
#'   subject_Id = "protein_Id",
#'   contrastfun = my_contrast_V2)
#'
#' mmm <- moderated_p_limma_long(factor_levelContrasts, group_by_col = "lhs")
#'
moderated_p_limma_long <- function(mm ,
                                   group_by_col = "lhs",
                                   estimate = "estimate",
                                   robust = FALSE){
  dfg <- mm |>
    dplyr::group_by_at(group_by_col) |>
    dplyr::group_split()
  xx <- purrr::map_df(dfg, moderated_p_limma, estimate = estimate, robust = robust)
  return(xx)
}


#' adjust columns
#'
#' @export
#' @keywords internal
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
    mm,
    column = "p.value",
    group_by_col = "contrast",
    newname = "FDR"
){
  dfg <- mm |>
    dplyr::group_by_at(group_by_col)
  xx <- dplyr::mutate(dfg, !!newname := p.adjust(!!sym(column),method = "BH"))
  return(xx)
}


# HELPER ----



# ROPECA ----

#' p-value of protein from p.value of the median fold change peptide.
#' @param max.n limit number of peptides per protein.
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#' plot(get_p_values_pbeta(0.1,1:10,10), ylim=c(0,0.1))
#' plot(get_p_values_pbeta(0.1,1:10,3), ylim=c(0,0.1))
#' plot(get_p_values_pbeta(0.3,1:30, 3), ylim=c(0,0.1))
#' abline(h=.05,col = 2)
#' plot(seq(0.0,1.0,length=30),get_p_values_pbeta(seq(0.0,1.0,length=30),rep(10,30)))
#' abline(0,1)
#' plot(seq(0.0,1.0,length=30),get_p_values_pbeta(seq(0.0,1.0,length=30),rep(10,30),3))
#' abline(0,1)
#' testthat::expect_equal(get_p_values_pbeta(0.3,10, 3),0.216, tolerance = 1e-4)
#' testthat::expect_equal(get_p_values_pbeta(0,10, 3),0, tolerance = 1e-4)
#' testthat::expect_equal(get_p_values_pbeta(1,10, 3),1, tolerance = 1e-4)
#' testthat::expect_equal(get_p_values_pbeta(1,10, 3),get_p_values_pbeta(1,3, 10), tolerance = 1e-4)
#'
get_p_values_pbeta <- function(median.p.value,
                               n.obs,
                               max.n = 10){
  n.obs <- pmin(n.obs, max.n)

  shape1 <- (n.obs/2 + 0.5)
  shape2 <- (n.obs - (n.obs/2 + 0.5) + 1)
  # n.obs/2 + 0.5

  stopifnot(shape1 == shape2)
  res.p.value <- pbeta(median.p.value,
                       shape1 = shape1,
                       shape2 = shape2)
  return(res.p.value)
}



#' compute protein level fold changes and p.values (using beta distribution)
#' takes p-value of the scaled p-value
#'
#' @param contrasts_data data frame
#' @param contrast name of column with contrast identifier
#' @param subject_Id name of column with typically protein Id
#' @param estimate name of column with effect size estimate
#' @param statistic statistic name of column with statistic (typically t-statistics)
#' @param p.value name of column with moderated.p.value
#' @param max.n used to limit the number of peptides in probablity computation.
#' @export
#' @family modelling
#' @keywords internal
#' @return data.frame with columns
#'
#'
#' @examples
#'
#' set.seed(10)
#' nrPep <- 10000
#' nrProtein <- 800
#' p.value <- runif(nrPep)
#' estimate <- rnorm(nrPep)
#' avgAbd <- runif(nrPep)
#' protein_Id <- sample(1:800, size = nrPep,
#'   replace = TRUE, prob = dexp(seq(0,5,length = 800)))
#'
#' plot(table(table(protein_Id)))
#'
#' testdata <- data.frame(contrast = "contrast1",
#'   protein_Id = protein_Id,
#'   estimate = estimate,
#'   pseudo_estimate = estimate,
#'   p.value = p.value,
#'   avgAbd = avgAbd )
#'
#' xx30 <- summary_ROPECA_median_p.scaled(testdata,
#'                                     subject_Id = "protein_Id",
#'                                     estimate = "estimate",
#'                                     p.value = "p.value",
#'                                     max.n = 30)
#'
#' xx2 <- summary_ROPECA_median_p.scaled(testdata,
#'                                     subject_Id = "protein_Id",
#'                                     estimate = "estimate",
#'                                     p.value = "p.value",
#'                                     max.n = 1)
#'
#' testthat::expect_equal(mad(xx2$estimate, na.rm = TRUE),0.384409, tolerance = 1e-4)
#' testthat::expect_equal(median(xx2$estimate), -0.006874857, tolerance = 1e-4)
#' testthat::expect_equal(xx2$beta.based.significance[1],0.819, tolerance = 1e-3)
#' testthat::expect_equal(xx2$beta.based.significance[2],0.9234362,tolerance = 1e-3)
#'
#' # Uniform distribution
#' hist(testdata$p.value)
#' hist(xx30$median.p.scaled, breaks = 20)
#' hist(xx2$median.p.scaled, breaks = 20)
#' # shows that beta.based.significance has NO uniform distribution
#' # although H0 is true for all cases.
#'
#' hist(xx30$beta.based.significance, breaks = 20)
#' hist(xx2$beta.based.significance, breaks = 20)
#'
#' hist(xx2$median.p.value, breaks = 20)
#' hist(xx2$beta.based.significance, breaks = 20)
#' hist(estimate)
#'
summary_ROPECA_median_p.scaled <- function(
    contrasts_data,
    contrast = "contrast",
    subject_Id = "protein_Id",
    estimate = "diff",
    statistic = "statistic",
    p.value = "moderated.p.value",
    max.n = 10){

  nrpepsPerProt <- contrasts_data |>
    group_by_at(c(subject_Id, contrast)) |>
    dplyr::summarize(n = dplyr::n() )

  contrasts_data <- contrasts_data |>
    dplyr::mutate(
      scaled.p =
        ifelse(!!sym(estimate) > 0, 1 - !!sym(p.value) , !!sym(p.value) - 1))

  summarized.protein <- contrasts_data |>
    group_by_at(c(subject_Id, contrast)) |>
    dplyr::summarize(
      n_not_na = n(),
      mad.estimate = mad(!!sym(estimate), na.rm = TRUE),
      estimate = median(!!sym(estimate), na.rm = TRUE),
      statistic = median(!!sym(statistic), na.rm = TRUE),
      median.p.scaled = median(.data$scaled.p, na.rm = TRUE),
      avgAbd = median(.data$avgAbd, na.rm = TRUE))

  summarized.protein <- summarized.protein |>
    dplyr::mutate(median.p.value = 1 - abs(.data$median.p.scaled))

  summarized.protein <- summarized.protein |>
    dplyr::mutate(beta.based.significance = get_p_values_pbeta(.data$median.p.value, .data$n_not_na, max.n = max.n))
  summarized.protein <- summarized.protein |>
    dplyr::mutate(n.beta = pmin(.data$n_not_na, max.n))

  summarized.protein <- dplyr::inner_join(nrpepsPerProt,
                                          summarized.protein,
                                          by = c(subject_Id, contrast))

  summarized.protein$isSingular <- FALSE
  # scale it back here.
  return(ungroup( summarized.protein ))
}


#' Fishers exact test on a datframe
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#' Nprot <- 1000
#' condA <- 8
#' condB <- 8
#' observedA <- sample(0:8, Nprot, replace = TRUE)
#' observedB <- sample(0:8, Nprot, replace = TRUE)
#' xb <- data.frame(observedA = observedA, observedB = observedB)
#'
#' xb$samplesA <- condA
#' xb$samplesB <- condB
#' proteinID <- unique(stringi::stri_rand_strings(Nprot + 20,5))[1:Nprot]
#' xb$proteinID <- proteinID
#' xlater <- xb
#' res <- contrasts_fisher_exact(xlater)
#'
contrasts_fisher_exact <- function(
    xb,
    observedA = "observedA",
    observedB = "observedB",
    samplesA = "samplesA",
    samplesB = "samplesB"
) {
  relativeRisk <- function(observedA, observedB, samplesA, samplesB) {
    rr <- (observedA/(observedA + observedB))/(samplesA/(samplesA + samplesB))
    return(rr)
  }
  odsRatio <- function(observedA, observedB, samplesA, samplesB) {
    rr <- (observedA/observedB)/(samplesA/samplesB)
    return(rr)
  }
  apply_fischer <- function(proteinID,observedA, observedB, samplesA, samplesB){
    mat <- matrix(c(observedA, samplesA - observedA,
                    observedB,samplesB - observedB), nrow = 2)
    fisher_result <- fisher.test(mat)
    return(data.frame(proteinID =  proteinID,
                      p_value = fisher_result$p.value,
                      OdsRatio = (fisher_result$estimate),
                      conf.lower = (fisher_result$conf.int[1]),
                      conf.higher = (fisher_result$conf.int[2]))
    )
  }

  xb$OdsRatioM <- odsRatio(
    observedA = xb[["observedA"]],
    observedB = xb[["observedB"]],
    samplesA = xb[["samplesA"]],
    samplesB = xb[["samplesB"]])
  xb$relativeRiskM <- relativeRisk(
    observedA = xb[["observedA"]],
    observedB = xb[["observedB"]],
    samplesA = xb[["samplesA"]],
    samplesB = xb[["samplesB"]])

  res <- vector(mode = "list", length(nrow(xb)))

  for (i in seq(nrow(xb))) {
    res[[i]] <- apply_fischer(
      xb[["proteinID"]][i],
      xb[["observedA"]][i],
      xb[["observedB"]][i],
      xb[["samplesA"]][i],
      xb[["samplesB"]][i] )
  }

  result <- dplyr::bind_rows(res)
  xx <- dplyr::inner_join(xb , result, by = c("proteinID" = "proteinID"))
  return(xx)
}



