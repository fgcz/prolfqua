#' preprocess
#' @export
ms_bench_preprocess <- function(data) {
  tmp <- data %>%
    ungroup() %>%
    mutate(ss  = case_when(
      grepl("HUMAN", protein_Id) ~ "HUMAN",
      grepl("ECOLI", protein_Id) ~ "ECOLI",
      TRUE ~ "OTHER"
    ))
  res <- tmp %>% dplyr::filter(!ss == "OTHER")

  res <- res %>% mutate(TP = (ss == "ECOLI"))
  return(list(data = res , table = table(tmp$ss)))
}


#' adds FDR and TPR plots
#'
#' FDP - false discovery proportion (Q in Benjamini hochber table)
#' FPR - false positive rate
#' TPR - true positive rate
#' TP_hits - true positives
#'
#' @export
ms_bench_add_FPRTPR <- function(data,
                                TP_col = "TP",
                                arrangeby = "estimate",
                                desc = TRUE){
  #data <- est

  data <- if (!desc) {
    data %>% arrange(!!sym(arrangeby))#,desc(!!sym(TPcol)))
  }else{
    data %>% arrange(desc(!!sym(arrangeby)))#,desc(!!sym(TPcol)))
  }
  data <- data %>% select(scorecol = !!sym(arrangeby) ,TP_col)
  data$what <- arrangeby
  data$F_ <- sum(!data$TP)
  data$T_ <- sum(data$TP)

  data <- na.omit(data)
  data <- data %>% mutate(
    R = 1:n()
    ,FDP = cummean(!TP)
    , TP_hits = cumsum(TP)
    , FN_hits = T_ - TP_hits
    , FP_hits = cumsum(!TP)
    , TN_hits = F_ - FP_hits
    , FPR = FP_hits/F_
    , TPR  = TP_hits/T_
    , ACC = (TP_hits + TN_hits)/(T_ + F_)

  ) %>% ungroup
  return(data)
}


#' computes auc and pauc given output from ms_bench_add_FPRTPR
#' @export
ms_bench_auc <- function(FPR, TPR, fpr_threshold = 1){
  # make sure that sorted.
  oFPR <- order(FPR)
  FPR <- FPR[oFPR]
  TPR <- TPR[oFPR]

  idx <- FPR < fpr_threshold
  TPR <- TPR[idx]
  FPR <- FPR[idx]
  #integrate
  res <- 1/2*sum(diff(FPR) * (head(TPR,-1) + tail(TPR,-1)))
  return(res)
}


