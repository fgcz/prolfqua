# function for modelling go here.
# rocs helper function
.rocs <- function(data ,response, predictor){
  responseX <- data %>% dplyr::pull(!!sym(response))
  predictorX <- data %>% dplyr::pull(!!sym(predictor))
  levels = levels(as.factor(responseX))
  if(length(levels) < 2){
    return(NULL)
  }
  res <- list()
  comparisons <- combn(levels, 2)
  for(i in 1:ncol(comparisons)){
    comp <- comparisons[,i]
    res[[i]] <-  tryCatch(pROC::roc(response = responseX,
                                    predictor = predictorX , levels=comp), error = function(x) NULL)
  }
  return(res)
}

#' Apply roc analysis on main factor on lowest level
#' @export
#' @importFrom purrr map
#' @examples
#'
#' library(tidyverse)
#' library(LFQService)
#' config <- spectronautDIAData250_config$clone(deep=TRUE)
#' config$parameter$min_nr_of_notNA  <- 20
#' data <- spectronautDIAData250_analysis
#' x <- sample(data$protein_Id,2)
#' data <- data %>% dplyr::filter(protein_Id %in% x)
#' res <- compute_roc(data, config)
#' dim(res)
#' head(res)
#' i <- 2
#'
#' pROC::plot.roc(res$rocs[[i]], print.auc = TRUE, main = paste(res$protein_Id[[i]], "\n",paste(res$rocs[[i]]$levels, collapse = " vs ")))
#' unique(res$protein_Id)
#'
compute_roc <- function(data, config){
  nested <- data %>% dplyr::group_by(!!sym(config$table$hierarchyKeys()[1]) ,
                                     !!sym(config$table$hierarchyKeys(TRUE)[1])) %>% nest()
  nested <- nested %>% dplyr::mutate(rocs = map(data ,
                                                .rocs, response = config$table$factorKeys()[1],
                                                predictor= config$table$getWorkIntensity() ))

  nested <- nested %>% dplyr::mutate(cls = map_lgl(rocs, is.null))  %>%
    dplyr::filter(cls == FALSE)
  #nested <- nested %>% mutate(names = map(rocs, names))

  dumm <- nested %>% dplyr::select(!!sym(config$table$hierarchyKeys()[1]),
                                   !!sym(config$table$hierarchyKeys(TRUE)[1]),
                                   rocs) %>%  tidyr::unnest()
  dumm <- dumm %>% dplyr::mutate(comparison = purrr::map_chr(rocs, function(x){paste(x$levels, collapse = " ")}))
  dumm <- dumm %>% tidyr::separate(comparison, into = c("response1" , "response2"), sep=" ")
  dumm <- dumm %>% dplyr::mutate(auc = map_dbl(rocs, pROC::auc)) %>%
    arrange(desc(auc))
  return(dumm)
}
