# function for modelling go here.
# rocs helper function
.rocs <- function(data ,response, predictor){
  responseX <- data %>% dplyr::pull(!!sym(response))
  predictorX <- data %>% dplyr::pull(!!sym(predictor))
  levels = levels(as.factor(responseX))
  if (length(levels) < 2) {
    return(NULL)
  }
  res <- list()
  comparisons <- combn(levels, 2)
  for (i in 1:ncol(comparisons)) {
    comp <- comparisons[,i]
    res[[i]] <-  tryCatch(pROC::roc(response = responseX,
                                    predictor = predictorX,
                                    direction = "<",
                                    levels = comp),
                          error = function(e) warning(e))
  }
  return(res)
}

#' Apply roc analysis on main factor on lowest hierarchy level
#'
#' deprecate function.
#' @param data data
#' @param config AnalysisConfiguration
#' @export
#' @keywords internal
#' @family deprecated
#' @examples
#'
#' FIXED <- FALSE
#' if(FIXED){
#'
#' bb <- prolfqua::data_ionstar$normalized()
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#' x <- sample(data$protein_Id,2)
#' data <- data %>% dplyr::filter(protein_Id %in% x)
#' res <- compute_roc(na.omit(data), config)
#' i <- 1
#'
#' pROC::plot.roc(res$rocs[[i]], print.auc = TRUE,
#'  main = paste(res$protein_Id[[i]], "\n",paste(res$rocs[[i]]$levels, collapse = " vs ")))
#' unique(res$protein_Id)
#' }
#'
compute_roc <- function(data, config){
  nested <- data %>% dplyr::group_by(!!!syms(config$table$hierarchyKeys())) %>% nest()
  nested <- nested %>% dplyr::mutate(rocs = map(data ,
                                                .rocs,
                                                response = config$table$fkeysDepth(),
                                                predictor = config$table$getWorkIntensity() ))

  nested <- nested %>% dplyr::mutate(cls = map_lgl(.data$rocs, is.null))  %>%
    dplyr::filter(cls == FALSE)

  #nested <- nested %>% mutate(names = map(rocs, names))

  dumm <- nested %>% dplyr::select(!!!syms(config$table$hierarchyKeys()),
                                   .data$rocs) %>%  tidyr::unnest(cols = c("rocs"))

  dumm <- dumm %>%
    dplyr::mutate(comparison =
                    purrr::map_chr(.data$rocs, function(x){paste(x$levels, collapse = " ")}))
  dumm <- dumm %>% tidyr::separate(.data$comparison, into = c("response1" , "response2"), sep = " ")
  dumm <- dumm %>% dplyr::mutate(auc = map_dbl(.data$rocs, pROC::auc)) %>%
    arrange(desc(.data$auc))
  return(dumm)
}
