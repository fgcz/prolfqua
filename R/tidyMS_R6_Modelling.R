# Creating models from configuration ----

.ehandler = function(e){
  warning("WARN :", e)
  # return string here
  as.character(e)
}


#' Create custom lmer model
#' @param modelstr model formula
#' @param model_name name of model
#' @param report_columns columns to report
#' @export
#' @family modelling
#' @examples
#'
#' tmp <- strategy_lmer("Intensity ~ condition + (1|peptide_Id)", model_name = "random_example")
#' tmp$model_fun(get_formula = TRUE)
#' tmp$isSingular
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
    modelTest <- tryCatch( lmerTest::lmer( formula , data = x ),
                          error = .ehandler)
    return(modelTest)
  }
  res <- list(model_fun = model_fun,
              isSingular = lme4::isSingular,
              contrast_fun = my_contest,
              model_name = model_name,
              report_columns = report_columns,
              is_mixed = TRUE)
  return(res)
}

#' Create custom lm model
#' @export
#' @param modelstr model formula
#' @param model_name name of model
#' @param report_columns columns to report
#' @family modelling
#' @examples
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
              is_mixed = FALSE)
  return(res)
}


#' Create custom quasibinomial glm model
#' @export
#' @param modelstr model formula
#' @param model_name name of model
#' @param report_columns columns to report
#' @family modelling
#' @examples
#' tmp <- strategy_glm("Intensity ~ condition", model_name = "parallel design")
#' tmp$model_fun(get_formula = TRUE)
#' tmp$isSingular
strategy_glm <- function(modelstr,
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
    modelTest <- tryCatch(glm( formula ,
                               data = x ,
                               family = stats::quasibinomial),
                          error = .ehandler)
    return(modelTest)
  }
  res <- list(model_fun = model_fun,
              isSingular = isSingular_lm,
              contrast_fun = my_contrast_V2,
              model_name = model_name,
              report_columns = report_columns,
              is_mixed = FALSE)
  return(res)
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
    if (df.residual(m) > 0) {
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
  modelProteinF <- modelProteinF %>% dplyr::filter(!!sym("exists_lmer") == TRUE)
  modelProteinF <- modelProteinF %>% dplyr::filter(nrcoeff_not_NA == max(nrcoeff_not_NA)) %>%
    dplyr::arrange(dplyr::desc(nrcoeff_not_NA))
  # modelProteinF <- modelProteinF %>% dplyr::filter(df.residual > 0)
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
#' library(tidyverse)
#' library(prolfqua)
#' ionstar <- prolfqua::data_ionstar$normalized()
#' ionstar$data <- ionstar$data %>% dplyr::filter(protein_Id %in% sample(protein_Id,10))
#' prolfqua::table_factors(ionstar$data, ionstar$config)
#' formula_randomPeptide <-
#'   strategy_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id)")
#' mr <- model_analyse( ionstar$data,
#'  formula_randomPeptide,
#'  subject_Id = ionstar$config$table$hkeysDepth())
#' get_complete_model_fit(mr$modelProtein)
model_analyse <- function(pepIntensity,
                          modelFunction,
                          subject_Id = "protein_Id",
                          modelName = "Model")
{
  pepIntensity %>%
    dplyr::group_by(!!!syms(subject_Id)) %>%
    tidyr::nest() -> nestProtein

  lmermodel <- "linear_model"

  pb <- progress::progress_bar$new(total = nrow(nestProtein))
  modelProtein <- nestProtein %>%
    dplyr::mutate(!!lmermodel := purrr::map(data, modelFunction$model_fun, pb = pb))

  modelProtein <- modelProtein %>% dplyr::mutate(!!"exists_lmer" := purrr::map_lgl(!!sym(lmermodel), function(x){!is.character(x)}))

  modelProteinF <- modelProtein %>%
    dplyr::filter( !!sym("exists_lmer") == TRUE)
  modelProteinF <- modelProteinF %>%
    dplyr::mutate(!!"isSingular" := purrr::map_lgl(!!sym(lmermodel), modelFunction$isSingular ))
  modelProteinF <- modelProteinF %>%
    dplyr::mutate(!!"df.residual" := purrr::map_dbl(!!sym(lmermodel), df.residual ))
  modelProteinF <- modelProteinF %>%
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

  modelProteinF <- modelProteinF %>% dplyr::mutate(nrcoef = purrr::map_int(!!sym(lmermodel), nrcoeff))
  modelProteinF <- modelProteinF %>% dplyr::mutate(nrcoeff_not_NA = purrr::map_int(!!sym(lmermodel), nrcoeff_not_NA))

  #return(list(modelProtein = modelProtein, modelProteinF = modelProteinF))
  modelProteinF <- modelProteinF %>%
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
#' m <- prolfqua::data_interactionModel_p1807
#' plot_lmer_peptide_predictions(m)
plot_lmer_peptide_predictions <- function(m){
  data <- m@frame
  data$prediction <- predict(m)
  interactionColumns <- intersect(attributes(terms(m))$term.labels,colnames(data))
  data <- make_interaction_column(data, interactionColumns, sep = ":")
  gg <- ggplot(data, aes(x = interaction , y = transformedIntensity)) + geom_point()
  gg <- gg + geom_point(aes(x = interaction, y = prediction), color = 2) + facet_wrap(~peptide_Id)
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
#' library(tidyverse)
#' library(prolfqua)
#' m <- prolfqua::data_basicModel_p1807
#' plot_lmer_peptide_noRandom(m)
#'
#' m <- prolfqua::data_interactionModel_p1807
#' plot_lmer_peptide_noRandom(m)
plot_lmer_peptide_noRandom <- function(m,legend.position = "none"){
  data <- m@frame
  ran <- lme4::ranef(m)[[1]]
  randeffect <- setdiff(all.vars( terms(formula(m)) ) , all.vars(terms(m)))
  ran <- tibble::as_tibble(ran,rownames = randeffect)
  colnames(ran) <- gsub("[()]","",colnames(ran))
  head(data)
  head(ran)
  ran <- dplyr::inner_join(data, ran, by = randeffect)

  ran <- ran %>% dplyr::mutate(int_randcorrected  = transformedIntensity  - Intercept)
  interactionColumns <- intersect(attributes(terms(m))$term.labels,colnames(data))
  ran <- make_interaction_column(ran,interactionColumns, sep = ":" )

  meanx <- function(x){mean(x,na.rm = TRUE)}
  gg <- ggplot(ran,aes(x = interaction , y= int_randcorrected, color = peptide_Id)) +
    geom_point(position = position_jitterdodge())
  gg <- gg + stat_summary(fun = meanx, colour = "black", geom = "point",
                          shape = 12, size = 3,show.legend = FALSE)
  gg <- gg + theme(axis.text.x = element_text(angle = -90, hjust = 0),
                   legend.position = legend.position)
  gg <- gg + geom_boxplot(alpha = 0.1)

  return(gg)
}

#' plot intensities per interaction with Two
#'  independent random effects removed (1|A) + (1|B)
#'  @param m model
#'  @param legend.position none
#'  @param firstlast todo
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#' #todo look up example
#'
plot_lmer_peptide_noRandom_TWO <- function(m, legend.position = "none", firstlast = TRUE){

  updateDataWithRandom <- function(data, m, i, randeffect){
    rand_i <- randeffect[i]
    ran <- lme4::ranef(m)[[rand_i]]
    name <- paste0(gsub("[()]","",colnames(ran)),"_", rand_i)
    colnames(ran) <- name
    ran <- tibble::as_tibble(ran,rownames = rand_i)
    ran <- dplyr::inner_join(data, ran, by = rand_i)
    ran_res <- ran %>% dplyr::mutate(int_randcorrected  = transformedIntensity  - !!sym(name))
    ran_res
  }

  data <- m@frame
  if (firstlast) {
    i1 <- 1;  i2 <- 2
  } else {
    i1 <- 2;  i2 <- 1
  }
  randeffect <- setdiff(all.vars( terms(formula(m)) ) , all.vars(terms(m)))
  ran_res <- updateDataWithRandom(data, m, i1, randeffect)
  ran_res <- updateDataWithRandom(ran_res, m, i2, randeffect)

  #randeffect <- setdiff(all.vars( terms(formula(m)) ) , all.vars(terms(m)))
  #rand_i <- randeffect[i]
  #ran <- ranef(m)[[rand_i]]
  #name <- paste0(gsub("[()]","",colnames(ran)),"_", rand_i)
  #colnames(ran) <- name
  #ran <- as_tibble(ran,rownames = rand_i)
  #ran_res <- dplyr::inner_join(ran_res, ran, by = rand_i)
  #ran_res <- ran_res %>% dplyr::mutate(int_randcorrected  = int_randcorrected  - !!sym(name))
  interactionColumns <- intersect(attributes(terms(m))$term.labels,colnames(data))
  ran_res <- make_interaction_column(ran_res,interactionColumns, sep = ":" )

  meanx <- function(x){mean(x,na.rm = TRUE)}
  gg <- ggplot(ran_res,aes(x = interaction , y = int_randcorrected,
                           color = !!sym(randeffect[i1]))) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2))
  gg <- gg + stat_summary(fun = meanx, colour = "black", geom = "point",
                          shape = 12, size = 3,show.legend = FALSE)
  gg <- gg + theme(axis.text.x = element_text(angle = -90, hjust = 0),legend.position = legend.position)
  gg <- gg + geom_boxplot(alpha = 0.1)
  gg
  return(gg)
}


#' Add predicted values for each interaction
#' @export
#' @keywords internal
#' @family modelling
#' @examples
#' m <- prolfqua::data_interactionModel_p1807
#' plot_lmer_predicted_interactions(plot_lmer_model_and_data(m,"dumm"),m)
plot_lmer_predicted_interactions <- function(gg, m){
  cm <- .lmer4_coeff_matrix(m)
  xstart_end <- data.frame(xstart = rownames(cm$mm), xend = rownames(cm$mm))
  ystart_end <- data.frame(xend = rownames(cm$mm), ystart = rep(0, nrow(cm$mm)),
                           yend = cm$mm %*% cm$coeffs)
  segments <- dplyr::inner_join(xstart_end, ystart_end, by = "xend")
  gg <- gg + geom_segment(aes(x = xstart, y = ystart , xend = xend, yend =yend),
                          data = segments, color = "blue", arrow = arrow())
  return(gg)
}

#' Make model plot with title - protein Name.
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#' library(prolfqua)
#' m <- prolfqua::data_interactionModel_p1807
#' plot_lmer_model_and_data(m,"dumm")
#'
plot_lmer_model_and_data <- function(m, proteinID, legend.position = "none"){
  gg <- plot_lmer_peptide_noRandom(m,legend.position = legend.position)
  gg <- plot_lmer_predicted_interactions(gg, m)
  gg <- gg + ggtitle(proteinID)
  gg
}


#' Plotting two independent random effects
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#' #todo
#'
plot_lmer_model_and_data_TWO <- function(m,
                                         proteinID,
                                         legend.position = "none" ,
                                         firstlast= TRUE){
  gg <- plot_lmer_peptide_noRandom_TWO(m, legend.position = legend.position, firstlast = firstlast)
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

  coeffs <- coefficients(summary(m))[,'Estimate']

  inter <- unique(data$interaction)
  mm <- matrix(0, nrow = length(inter), ncol = length(coeffs))
  rownames(mm) <- inter
  colnames(mm) <- names(coeffs)
  mm[,1] <- 1
  coefi <- coeffs[-1]
  for (i in 1:length(coefi)) {
    positionIDX <- grep(names(coefi)[i], inter)
    mm[positionIDX, i + 1 ] <- 1
  }
  return(list(mm = mm, coeffs = coeffs))
}


.get_match_idx <- function(mm, factor_level){
  ddd <- split2table(rownames(mm), split = ":")
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
#' #if(FALSE){
#'
#' m <- prolfqua::data_basicModel_p1807
#' m
#' linfct <- linfct_from_model(m)
#'
#' linfct$linfct_factors
#' linfct$linfct_interactions
#'
#' m <- prolfqua::data_interactionModel_p1807
#' linfct <- linfct_from_model(m)
#' all.equal(linfct$linfct_factors["CelltypeCMP/MEP",] ,
#'  apply(linfct$linfct_interactions[grep("CelltypeCMP/MEP", rownames(linfct$linfct_interactions)),],2, mean))
#' linfct$linfct_interactions
#' #}
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
  linfct_factors <- l_factors %>%
    dplyr::select(-factor_level) %>%
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
#' m <- prolfqua::data_basicModel_p1807
#' linfct <- linfct_from_model(m,as_list = FALSE)
#' linfct
#'
#' Contrasts <- c("CMP/MEP - HSC" = "`CelltypeCMP/MEP` - `CelltypeHSC`",
#' "someWeird" = "`class_therapyc.NO:CelltypeCMP/MEP` - `class_therapyp.HU:CelltypeCMP/MEP`")
#' linfct_matrix_contrasts(linfct, Contrasts )
#'
linfct_matrix_contrasts <- function(linfct , contrasts, p.message = FALSE){
  linfct <- t(linfct)
  df <- tibble::as_tibble(linfct, rownames = "interaction")
  make_contrasts <- function(data,
                             contrasts)
  {
    cnams <- setdiff(colnames(data),"interaction")
    for (i in 1:length(contrasts)) {
      if (p.message) {message(names(contrasts)[i], "=", contrasts[i],"\n")}
      data <- dplyr::mutate(data, !!names(contrasts)[i] := !!rlang::parse_expr(contrasts[i]))
    }

    res <- data %>% dplyr::select(-one_of(cnams))
    return(res)
  }

  res <- make_contrasts(df, contrasts )
  res <- column_to_rownames(res,"interaction")
  res <- t(res)
  return(res)
}


#' create all possible contrasts
#' @export
#' @keywords internal
#' @family modelling
#' @examples
#' m <- prolfqua::data_basicModel_p1807
#' m
#' linfct <- linfct_from_model(m)
#'
#' xl <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_factors)
#' xx <- prolfqua::linfct_all_possible_contrasts(linfct$linfct_interactions)
#'
linfct_all_possible_contrasts <- function(lin_int ){
  combs <- combn(nrow(lin_int),2)
  names <- rownames(lin_int)
  newnames <- rep("", ncol(combs))
  new_lin_fct <- matrix(NA,  nrow = ncol(combs), ncol = ncol(lin_int))
  for (i in 1:ncol(combs)) {
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
#' library(prolfqua)
#' m <- prolfqua::data_basicModel_p1807
#' xl <- linfct_factors_contrasts(m)
#' xl
#' m <- lm(Petal.Width ~ Species, data = iris)
#' linfct_factors_contrasts(m)
linfct_factors_contrasts <- function(m){
  ffac <- attributes(terms(m))$term.labels
  ffac <- ffac[!grepl(":",ffac)] # remove interactions
  linfct_factors <- linfct_from_model(m)$linfct_factors

  factorDepths <- rownames(linfct_factors)
  res <- vector(length(ffac), mode = "list")
  for (i in 1:length(ffac)) {
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
#' library(prolfqua)
#' mb <- prolfqua::data_basicModel_p1807
#' linfct <- linfct_from_model(mb)
#' names(linfct)
#' my_glht(mb, linfct$linfct_factors)
#'
#' m <- prolfqua::data_modellingResult_A$modelProtein$linear_model[[1]]
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
    for (i in 1:nrow(linfct)) {
      x <- multcomp::glht(model, linfct = linfct[i,,drop = FALSE])
      RHS <- broom::tidy(confint(x)) %>% dplyr::select(-estimate)

      RHS$df <- x$df
      RHS$sigma <- sigma(model)

      x <- dplyr::inner_join(broom::tidy(summary(x)),RHS,by = c("contrast"))# %>% dplyr::select(-contrast)
      res[[i]] <- x
    }
    res <- dplyr::bind_rows(res)
    return(res)
  }else{
    x <- multcomp::glht(model, linfct = linfct)
    RHS <- broom::tidy(confint(x)) %>% dplyr::select(-estimate)
    RHS$df <- x$df
    RHS$sigma <- sigma(model)
    res <- dplyr::inner_join(broom::tidy(summary(x)), RHS, by = c("contrast")) %>%
      dplyr::select(-rhs)
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
#' m <- prolfqua::data_modellingResult_A$modelProtein$linear_model[[1]]
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
#' m <- prolfqua::data_modellingResult_A$modelProtein$linear_model[[1]]
#' linfct <- linfct_from_model(m)$linfct_factors
#' m
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
#' m <- prolfqua::data_modellingResult_A$modelProtein$linear_model[[1]]
#' linfct <- linfct_from_model(m)$linfct_factors
#' my_contrast_V2(m, linfct, confint = 0.95)
#' my_contrast_V2(m, linfct, confint = 0.99)
#'
my_contrast_V2 <- function(m, linfct,confint = 0.95){
  Sigma.hat <- vcov(m)
  coef <- coefficients(m)
  res <- vector(nrow(linfct), mode = "list")
  for (i in 1:nrow(linfct)) {
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
                              Sigma.hat = Sigma.hat_red,confint = confint)
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
#' library(prolfqua)
#' mb <- prolfqua::data_basicModel_p1807
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
  res <- res %>% dplyr::rename(estimate = Estimate,
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
#' dd <- prolfqua::data_factor_levelContrasts
#' head(dd)
#' tmp <- pivot_model_contrasts_2_Wide(dd, subject_Id = "Compound")
#' tmp
pivot_model_contrasts_2_Wide <- function(modelWithInteractionsContrasts,
                                         subject_Id = "protein_Id",
                                         columns = c("estimate", "p.value","p.value.adjusted"),
                                         contrast = "lhs"){

  m_spread <- function(longContrasts, subject_Id, column , contrast){
    res <- longContrasts %>%
      dplyr::select_at(c(subject_Id, "isSingular", contrast,  column))
    res <- res %>% dplyr::mutate(!!contrast := paste0(column, ".", !!sym(contrast)))
    res <- res %>% tidyr::spread(contrast, !!sym(column) )
    return(res)
  }
  res <- list()
  for (column in columns) {
    res[[column]] <- m_spread(modelWithInteractionsContrasts, subject_Id,column, contrast)
  }
  res <- res %>% reduce(left_join, by = c(subject_Id,"isSingular"))
  return(res)
}
#' compute group averages
#'
#' used in p2621, p2109
#'
#' @export
#' @keywords internal
#' @examples
#'
#' modelSummary_A <- prolfqua::data_modellingResult_A
#' m <- get_complete_model_fit(modelSummary_A$modelProtein)
#'
#' factor_contrasts <- linfct_factors_contrasts( m$linear_model[[1]])
#' factor_contrasts
#'
#' factor_levelContrasts <- contrasts_linfct( m,
#'         factor_contrasts,
#'         subject_Id = "Compound",
#'         contrastfun = prolfqua::my_contrast_V2)
#'
#' #usethis::use_data(factor_levelContrasts, overwrite = TRUE)
#'
#' data_models_interaction <- prolfqua::data_models_interaction
#'
#' m <- get_complete_model_fit(data_models_interaction$modelProtein)
#' m$linear_model[[1]]
#' factor_contrasts <- linfct_factors_contrasts( m$linear_model[[1]])
#' factor_levelContrasts <- contrasts_linfct( m,
#'                            factor_contrasts,
#'                        subject_Id = "protein_Id")
#' head(factor_levelContrasts)
#' m$linear_model[[1]]
#' my_contest(m$linear_model[[1]],factor_contrasts )
#'
#' plot(factor_levelContrasts$df, factor_levelContrasts$df.residual.model )
#' abline(c(0,1))
#' plot(factor_levelContrasts$df.residual.model , factor_levelContrasts$df - factor_levelContrasts$df.residual.model )
#'
contrasts_linfct <- function(models,
                             linfct,
                             subject_Id = "protein_Id" ,
                             contrastfun = prolfqua::my_contest){
  #computeGroupAverages
  message("computing contrasts.")
  modelcol <- "linear_model"
  models <- models %>% dplyr::filter(exists_lmer == TRUE)

  if ("list" %in% class(linfct)) {
    #stopifnot(length(linfct) == length(models[[modelcol]]))
    #interaction_model_matrix <- models %>%
    #  dplyr::mutate(contrast = purrr::map2(!!sym(modelcol) , linfct, function(x, y){ contrastfun(x, linfct = y) } ))

    stopifnot(length(linfct) == length(models[[modelcol]]))
    res <- vector(mode = "list", length = length(linfct))
    for (i in 1:length(linfct)) {
      res[[i]] <- contrast$contrastfun(models$linear_model[[i]], linfct = linfct[[i]])
    }

    rown <- lapply(res, rownames)
    intersect <- Reduce(intersect, rown)
    res <- lapply(res , function(x){ x[rownames(x) %in% intersect,]})
    interaction_model_matrix <- models
    interaction_model_matrix$contrast <- res
  } else {
    interaction_model_matrix <- models %>%
      dplyr::mutate(contrast = purrr::map(!!sym(modelcol) , contrastfun , linfct = linfct ))

  }

  mclass <- function(x){
    class(x)[1]
  }

  interaction_model_matrix <-  interaction_model_matrix %>%
    dplyr::mutate(classC = purrr::map_chr(contrast,mclass)) %>%
    dplyr::filter(classC != "logical")

  contrasts <- interaction_model_matrix %>%
    dplyr::select_at( c(subject_Id, "contrast") ) %>%
    tidyr::unnest_legacy()

  # take sigma and df from somewhere else.
  modelInfos <- models %>%
    dplyr::select_at(c(subject_Id,
                       "isSingular",
                       "sigma.model" = "sigma",
                       "df.residual.model" = "df.residual" )) %>%

    dplyr::distinct()
  contrasts <- dplyr::inner_join(contrasts, modelInfos, by = subject_Id)
  return(contrasts)
}



# LIMMA ----

#' Moderate p-values - limma approach
#' @export
#' @family modelling
#' @keywords internal
#'
moderated_p_limma <- function(mm, df = "df", robust = FALSE, confint = 0.95){
  sv <- limma::squeezeVar(mm$sigma^2, df = mm[[df]],robust = robust)
  sv <- tibble::as_tibble(sv)
  sv <- sv %>% setNames(paste0('moderated.', names(.)))
  mm <- dplyr::bind_cols(mm, sv)
  mm <- mm %>% dplyr::mutate(moderated.statistic  =  .data$statistic * .data$sigma /  sqrt(.data$moderated.var.post))
  mm <- mm %>% dplyr::mutate(moderated.df.total = !!sym(df) + .data$moderated.df.prior)
  mm <- mm %>% dplyr::mutate(moderated.p.value = 2*pt( abs(.data$moderated.statistic),
                                                       df = .data$moderated.df.total, lower.tail = FALSE) )

  prqt <- -qt((1 - confint)/2, df = mm$moderated.df.total)
  mm$moderated.conf.low <- mm$estimate  - prqt * sqrt(mm$moderated.var.post)
  mm$moderated.conf.high <- mm$estimate + prqt * sqrt(mm$moderated.var.post)
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
#' library(prolfqua)
#' modelSummary_A <- prolfqua::data_modellingResult_A
#' m <- get_complete_model_fit(modelSummary_A$modelProtein)
#' factor_contrasts <- linfct_factors_contrasts(m$linear_model[[1]])
#' factor_levelContrasts <- contrasts_linfct( modelSummary_A$modelProtein,
#'                                                    factor_contrasts,
#'                                                    subject_Id = "Compound",
#'                                                    contrastfun = my_contrast_V2)
#'
#' mmm <- moderated_p_limma_long(factor_levelContrasts, group_by_col = "lhs")
#' plot(mmm$p.value, mmm$moderated.p.value, log = "xy")
#' abline(0,1, col = 2)
#'
#' # updating lmer model
#' data_models_interaction <- prolfqua::data_models_interaction
#'
#' m <- get_complete_model_fit(data_models_interaction$modelProtein)
#' factor_contrasts <- linfct_factors_contrasts(m$linear_model[[1]])
#'
#' factor_levelContrasts <- contrasts_linfct(m,
#'                                          factor_contrasts,
#'                                          subject_Id = "protein_Id")
#'
#' mmm <- moderated_p_limma_long(factor_levelContrasts, group_by_col = "lhs")
#' head(mmm)
#' plot(mmm$p.value, mmm$moderated.p.value, log = "xy")
#' abline(0,1, col = 2)
#'
moderated_p_limma_long <- function(mm ,
                                   group_by_col = "lhs",
                                   robust = FALSE){
  dfg <- mm %>%
    dplyr::group_by_at(group_by_col) %>%
    dplyr::group_split()
  xx <- purrr::map_df(dfg, moderated_p_limma, robust = robust)
  return(xx)
}


#' adjust columns
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
#' data <- adjust_p_values(data)
#' Adata <- data %>% dplyr::filter(contrast == "A")
#' stopifnot(all.equal(Adata$p.value.adjusted, p.adjust(Adata$p.value, method="BH")))
#' data2 <- adjust_p_values(data, group_by_col = NULL)
#' stopifnot(all.equal(data2$p.value.adjusted, p.adjust(data2$p.value, method="BH")))
#'
#'
adjust_p_values <- function(mm,
                            column = "p.value",
                            group_by_col = "contrast",
                            newname = NULL
){
  if (is.null(newname)) {
    newname <- paste0(column, ".adjusted")
  }
  dfg <- mm %>%
    dplyr::group_by_at(group_by_col)

  xx <- dplyr::mutate(dfg, !!newname := p.adjust(!!sym(column),method = "BH"))

  return(xx)

}

#' write results of `contrasts_linfct`
#' @keywords internal
#' @export
#'
contrasts_linfct_write <- function(results,
                                   config,
                                   path,
                                   modelName = "Model",
                                   prefix = "Contrasts",
                                   columns = c("estimate", "p.value", "p.value.adjusted")){

  subject_Id <- config$table$hkeysDepth()

  if (!is.null(path)) {
    fileLong <- paste0(prefix, "_", modelName)
    message("Writing: ", fileLong, "\n")
    lfq_write_table(separate_hierarchy(results, config) , path = path , name = fileLong)
    fileWide <- paste0(prefix, "_", modelName, "_PIVOT")
    message("Writing: ", fileWide, "\n")
    resultswide <- pivot_model_contrasts_2_Wide(results,
                                                subject_Id = subject_Id,
                                                columns = columns)
    lfq_write_table(separate_hierarchy(resultswide, config), path = path, name = fileWide)
  }
}

# HELPER ----

#' get coefficients from all models
#' @export
#' @family modelling
#' @keywords internal
get_model_coefficients <- function(modeldata, config){
  l_coeff <- function(m){
    if (!is.null(m)) {
      res <- as.numeric(coefficients(m))
      return(res)
    }
    return(numeric())
  }

  n_coeff <- function(m){
    if (!is.null(m)) {
      res <- names(coefficients(m))
      return(res)
    }
    return(character())
  }

  xx <- modeldata %>%
    dplyr::mutate(coefficients_values = purrr::map(linear_model, l_coeff)) %>%
    dplyr::mutate(coefficients_names = purrr::map(linear_model, n_coeff))

  xxs <- xx %>% dplyr::select( config$table$hkeysDepth(),
                               coefficients_values,
                               coefficients_names)
  xxxn <- xxs %>% unnest_legacy()
  xxcoef <- xxxn %>% spread(coefficients_names,coefficients_values)
  return(xxcoef)
}


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
#' @param max.n used to limit the number of peptides in probablity computation.
#' @export
#' @family modelling
#' @keywords internal
#' @return data.frame with columns
#'
#'
#' @examples
#'
#' library(prolfqua)
#' library(tidyverse)
#' nrPep <- 10000
#' nrProtein <- 800
#' p.value <- runif(nrPep)
#' estimate <- runif(nrPep)
#' protein_Id <- sample(1:800, size = nrPep,
#'   replace = TRUE, prob = dexp(seq(0,5,length = 800)))
#'
#' plot(table(table(protein_Id)))
#'
#' testdata <- data.frame(contrast = "contrast1", protein_Id = protein_Id,
#'   estimate = estimate, pseudo_estimate = estimate, p.value = p.value )
#'   head(testdata)
#' xx30 <- summary_ROPECA_median_p.scaled(testdata,
#'                                     subject_Id = "protein_Id",
#'                                     estimate = "estimate",
#'                                     p.value = "p.value",
#'                                     max.n = 30)
#' xx30$mad.estimate
#' xx2 <- summary_ROPECA_median_p.scaled(testdata,
#'
#'                                     subject_Id = "protein_Id",
#'                                     estimate = "estimate",
#'                                     p.value = "p.value",
#'                                     max.n = 1)
#' mad(xx2$estimate, na.rm=TRUE)
#'
#' hist(testdata$p.value)
#' hist(xx30$median.p.scaled, breaks = 20)
#' hist(xx2$median.p.scaled, breaks = 20)
#' hist(xx30$beta.based.significance, breaks = 20)
#' hist(xx2$beta.based.significance, breaks = 20)
#' hist(xx2$median.p.value, breaks = 20)
#' hist(xx2$beta.based.significance, breaks = 20)
#' hist(xx2$mad.estimate)
#'
#' summary_ROPECA_median_p.scaled(prolfqua::data_exampleForRopeca, contrast = "contrast")
#' xx2$mad.estimate
#'
#'
summary_ROPECA_median_p.scaled <- function(
  contrasts_data,
  contrast = "contrast",
  subject_Id = "protein_Id",
  estimate = "estimate",
  statistic = "statistic",
  p.value = "moderated.p.value",
  max.n = 10){

  nrpepsPerProt <- contrasts_data %>%
    group_by_at(c(subject_Id, contrast)) %>%
    dplyr::summarize(n = n() )

  contrasts_data <- contrasts_data %>%
    # dplyr::filter(!is.na(!!sym(p.value))) %>%
    dplyr::mutate(scaled.p = ifelse(!!sym(estimate) > 0, 1 - !!sym(p.value) , !!sym(p.value) - 1))

  summarized.protein <- contrasts_data %>%
    group_by_at(c(subject_Id, contrast)) %>%
    dplyr::summarize(
      n_not_na = n(),
      mad.estimate = mad(!!sym(estimate), na.rm = TRUE),
      estimate = median(!!sym(estimate), na.rm = TRUE),
      statistic = median(!!sym(statistic), na.rm = TRUE),
      median.p.scaled = median(scaled.p, na.rm = TRUE))


  if (rlang::has_name(contrasts_data, "c1_name")) {
    ccsummary <- contrasts_data %>%
      group_by_at(c(subject_Id, contrast)) %>%
      dplyr::summarize(
        c1_name = unique(c1_name),
        c1 = median(c1, na.rm = TRUE),
        c2_name = unique(c2_name),
        c2 = median(c2, na.rm = TRUE) )
    summarized.protein <- inner_join(summarized.protein, ccsummary, by = c(subject_Id, contrast))
  }

  summarized.protein <- summarized.protein %>%
    dplyr::mutate(median.p.value = 1 - abs(median.p.scaled))

  summarized.protein <- summarized.protein %>%
    dplyr::mutate(beta.based.significance = get_p_values_pbeta(median.p.value, n_not_na, max.n = max.n))
  summarized.protein <- summarized.protein %>%
    dplyr::mutate(n.beta = pmin(n_not_na, max.n))

  summarized.protein <- dplyr::inner_join(nrpepsPerProt,
                                          summarized.protein,
                                          by = c(subject_Id, contrast))

  summarized.protein$isSingular <- FALSE
  # scale it back here.
  return(ungroup( summarized.protein ))
}



