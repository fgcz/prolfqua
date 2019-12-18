## functions supporting working with stand
#' Compute contrasts an
#' covert brms output
#' to code mcmc.list
#' @export
#'
ms_mcmc_constrast <- function(model, linfct_A){
  my_MCMC_contrast_1 <- function(x,linfct_A){
    x <- as_tibble(x) %>% dplyr::select(starts_with("b_"))
    linfct_A <- dd$linfct_A
    colnames(linfct_A) <- paste0("b_",gsub("[()]","",colnames(dd$linfct_A)))
    x <- select_at(x,colnames(linfct_A))
    output <-  coda::mcmc(as.matrix(x)%*% t(dd$linfct_A))
  }
  xx <- brms::as.mcmc(model)
  res <- coda::mcmc.list(lapply(xx, my_MCMC_contrast_1, dd$linfct_A))
  return(res)
}

#' plugin method for MCMCsummary. can be passed to MCMCsummary func and func_name argument.
#' @export
#'
#'
ms_mcmc_checkzero <- function(x = NULL){
  if(!is.null(x)){
    less0 <- mean(x < 0)
    great0 <- mean(x > 0)
    return(c( less0 = less0,
              great0 = great0,
              minp = min(less0, great0)))
  }else{
    return(c("less0", "great0", "minp"))
  }
}



### Model with stan
#' fits brms model and computes summary.
#' @export
ms_brms_model<- function( mdata,
                        memodel,
                        linfct_A,
                        func=ms_mcmc_checkzero,
                        summarize=TRUE){
  if(class(memodel) == "character"){
    resultmodel <- tryCatch(brm(memodel, data = mdata, cores= 6), error = function(x) NULL)
  }else if(class(memodel) == "brmsfit"){
    resultmodel <- tryCatch(update(memodel, newdata = mdata, cores= 6), error = function(x) NULL)
  }
  res <- NULL
  if(!is.null(resultmodel)){
    if(summarize){
      res <- ms_mcmc_constrast(resultmodel, linfct_A)
      res <- as_tibble(MCMCsummary(res, func = func, func_name = func()), rownames="contrast")
      return(res)
    } else{
      return(resultmodel)
    }
  }else{
    return(NULL)
  }
}
