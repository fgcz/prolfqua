#' adds FDR and TPR plots
#' @export
ms_bench_add_FPRTPR <- function(tmp, idcol = "protein_Id",
                       TPcol = "TP",
                       groupby = "contrast",
                       arrangby = "beta.based.significance",
                       type = c("probability","foldchange"),desc = FALSE){

  tmp$TP_total <- length(unique(tmp[[idcol]][ tmp[[TPcol]] == TRUE]))
  tmp <- tmp %>% group_by(!!!syms(groupby))
  tmp %>% summarise(n = n())

  tmp <- if(!desc){
    tmp %>% arrange(!!sym(arrangby))
  }else{
    tmp %>% arrange(desc(!!sym(arrangby)))
  }
  tmp <- tmp %>% mutate(FDP = cummean(!TP)
                        , FPR = cumsum(!TP)/sum(!TP)
                        , TPR  = cumsum(TP)/TP_total
                        , TP_hits = cumsum(TP)

  ) %>% ungroup
  res <- tmp %>%
    dplyr::select_at(c(idcol, TPcol, groupby, "FDP", "FPR", "TPR", "TP_hits"))
  res$arrangeby <- arrangby
  res$type <- match.arg(type)
  return(res)
}

#' computes auc and pauc given output from ms_bench_add_FPRTPR
#' @export
ms_bench_auc <- function(FPR, TPR, fpr_threshold = 1){
  idx <- FPR < fpr_threshold
  TPR <- TPR[idx]
  FPR <- FPR[idx]

  res <- 1/2*sum(diff(FPR) * (head(TPR,-1)+tail(TPR,-1)))
  return(res)
}


