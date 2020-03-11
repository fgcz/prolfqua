#' preprocess
#' @export
#'
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
  return(res/fpr_threshold*100)
}


#' scale_probabilities
#'@export
#'
scale_probabilities <- function(est ,estimate = "estimate", toscale){

  addScaledP <- function(data , estimate , scale){
    scaled.p = paste0("scaled.", scale)
    data <- data %>% dplyr::mutate(!!scaled.p := ifelse(!!sym(estimate) > 0,
                                                        1 - !!sym(scale) ,
                                                        !!sym(scale) - 1))
    return(data)
  }

  for (scale in toscale) {
    message(scale)
    est <- addScaledP( est, estimate = estimate , scale = scale)
  }
  return(est)
}


#' do_confusion
#' @export
#'
do_confusion <- function(data, arrangeby = c("pseudo_estimate", "estimate",  "statistic",  "scaled.p", "scaled.moderated.p" )){
  # TODO add to LFQService
  est <- data %>% dplyr::select_at(c("TP",
                                     "contrast",
                                     arrangeby))
  res <- list()
  for (arrange in arrangeby) {
    message(arrange)
    res[[arrange]] <- ms_bench_add_FPRTPR(est,TP_col = "TP", arrangeby = arrange, desc = TRUE )
  }
  all <- bind_rows(res)
  return(all)
}


#' Visualizes data frame with columns FPR, TPR, FDP
#' @export
#'
plot_FDR_summaries <- function(pStats, model_type = "mixed effects model"){
  summaryS <- pStats %>% dplyr::group_by(what) %>%
    dplyr::summarize(auc = ms_bench_auc(FPR, TPR), auc10 =  ms_bench_auc(FPR, TPR, 0.1), auc20 = ms_bench_auc(FPR, TPR, 0.2) )

  ftable <- flextable::flextable(summaryS) %>%
    flextable::set_caption(caption = paste0("AUC, and pAUC at 0.1 and 0.2 FPR for ", model_type))  %>%
    colformat_num(digits=2)

  sumd <- reshape2::melt(summaryS)
  barp <- ggplot(sumd, aes(x = what, y = value)) +
    geom_bar(stat = "identity",position = position_dodge()) +
    facet_wrap(~variable, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    coord_cartesian(ylim = c(floor(min(sumd$value)/10)*10,100))


  p1 <- ggplot(pStats , aes(x = FPR, y = TPR, color = what)) + geom_path()  + labs(tag = "A")
  p2 <- ggplot(pStats , aes(x = FPR, y = TPR, color = what)) + geom_path()  + labs(tag = "B") + xlim(0,0.2)
  p3 <- ggplot(pStats , aes(x = FDP, y = TPR, color = what)) + geom_path()  + labs(tag = "C") + xlim(0,0.2)
  rocp <- ggpubr::ggarrange(p1, p2, p3, nrow = 1,common.legend = TRUE, legend = "bottom")
  return(list(rocp = rocp, barp = barp, ftable = ftable, summaryS = summaryS))
}

#'summarise_missing_contrasts
#'@export
#'
summarise_missing_contrasts <- function(data,
                                        hierarchy = c("protein_Id"),
                                        what = "statistic"){
  xxA <- data %>%
    group_by_at(hierarchy) %>%
    summarize(n = n(), nr_na = sum(is.na(!!sym(what))))
  xx <- as.data.frame(table(xxA$nr_na))
  colnames(xx) <- c("nr_missing", "nr_Proteins")

  return(list(xx = xx, xxA = xxA ))
}

#' benchmark
#' @export
benchmark <- function(resXXmedpolishTSV,
                      relevantContrasts = NULL,
                      completeContrasts = TRUE,
                      toscale = c("p.value", "moderated.p.value"),
                      benchmark = c("pseudo_estimate", "estimate",  "statistic",  "scaled.p.value", "scaled.moderated.p.value" ),
                      model_type = "protein level measurments, linear model"
){
  res <- list()
  if (!is.null(relevantContrasts)) {
    resXXmedpolishTSV <- resXXmedpolishTSV %>% dplyr::filter(contrast %in% relevantContrasts)
  }

  prpr <- ms_bench_preprocess(resXXmedpolishTSV)
  smc <- summarise_missing_contrasts(prpr$data)
  res$smc <- smc
  prpr$data <- scale_probabilities(prpr$data, toscale = toscale)


  if ( completeContrasts ) {
    message("completing contrasts")
    print(dim(smc$xxA))
    print(dim(prpr$data))
    smc$xxA <- smc$xxA %>% filter(n == n - nr_na )
    prpr$data <- prpr$data %>% filter(protein_Id %in% smc$xxA$protein_Id)
    print(dim(smc$xxA))
    print(dim(prpr$data))

  }

  confusion <- do_confusion(prpr$data, arrangeby = benchmark)
  vissum <- plot_FDR_summaries(confusion,model_type = model_type)
  res$confusion <-  confusion
  res$vissum <- vissum
  return(res)
}
