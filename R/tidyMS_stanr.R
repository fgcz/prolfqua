## functions supporting working with stan and brms
#'
#' Compute contrasts an
#' covert brms output
#' to code mcmc.list
#' @param model brms model
#' @param linfct_A linear function
#' @export
#'
#'
ms_mcmc_constrast <- function(model, linfct_A){
  colnames(linfct_A) <- paste0("b_",gsub("[()]","",colnames(linfct_A)))
  my_MCMC_contrast_1 <- function(x,linfct_A){
    x <- as_tibble(x) %>% dplyr::select(starts_with("b_"))
    x <- select_at(x,colnames(linfct_A))
    output <-  coda::mcmc(as.matrix(x) %*% t(linfct_A))
  }
  xx <- brms::as.mcmc(model)
  res <- coda::mcmc.list(lapply(xx, my_MCMC_contrast_1, linfct_A))
  return(res)
}

#' plugin method for MCMCsummary.
#' can be passed to MCMCsummary func and func_name argument.
#' @export
#'
#'
ms_mcmc_checkzero <- function(x = NULL){
  if (!is.null(x)) {
    less0 <- mean(x < 0)
    great0 <- mean(x > 0)
    return(c( less0 = less0,
              great0 = great0,
              minp = min(less0, great0)))
  }else{
    return(c("less0", "great0", "minp"))
  }
}

# do not even try to fit if levels are missing
check_factors_level_coverage <- function(mdata26, fixeff){
  complete <- mdata26 %>%
    dplyr::select_at(fixeff) %>%
    dplyr::distinct()
  omitted <- mdata26 %>%
    na.omit %>%
    dplyr::select_at(fixeff) %>%
    dplyr::distinct()
  return(nrow(complete) == nrow(omitted))
}


### Model with stan
#' fits brms model and computes summary.
#' @param mdata todo
#' @param memodel todo
#' @param fixef todo
#' @param linfct_A todo
#' @param func todo
#' @param summarize todo
#' @param cores todo
#' @export
#' @examples
#' # todo
ms_brms_model <- function(mdata,
                          memodel,
                          fixef,
                          linfct_A,
                          func =ms_mcmc_checkzero,
                          summarize=TRUE,
                          cores = parallel::detectCores()
                          , iter = 2000,
                          control = list(adapt_delta = 0.95)){

  if (!check_factors_level_coverage(mdata, fixef)) {
    return(NULL)
  }

  if (class(memodel) == "character") {
    resultmodel <- tryCatch(
      brms::brm(memodel,
                data = mdata,
                cores = cores,
                refresh = 0,
                iter = iter,
                control = control),
      error = function(x) warning(x))
  }else if (class(memodel) == "brmsfit") {
    resultmodel <- tryCatch(
      update(memodel,
             newdata = mdata,
             cores = cores,
             refresh = 0,
             iter = iter,
             control = control),
      error = function(x) warning(x))
  }
  res <- NULL
  if (!is.null(resultmodel)) {
    if (summarize) {
      res <- ms_mcmc_constrast(resultmodel, linfct_A)
      if (is.null(res)) {return(NULL)}
      res <- as_tibble(MCMCvis::MCMCsummary(res, func = func, func_name = func()), rownames = "contrast")
      return(res)
    } else{
      return(resultmodel)
    }
  }else{
    return(NULL)
  }
}
