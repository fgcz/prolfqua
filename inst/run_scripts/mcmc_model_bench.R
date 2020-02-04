rm(list = ls())
library(LFQServiceAnalysisTemplate)
library(tidyverse)
library(dplyr)

preprocess <- function(data){
  tmp <- data %>%
    ungroup() %>%
    mutate(ss  = case_when(
      grepl("HUMAN",protein_Id) ~ "HUMAN",
      grepl("ECOLI", protein_Id ) ~ "ECOLI",
      TRUE ~ "OTHER"))
  tmp <- tmp %>% filter(!ss == "OTHER")
  tmp <- tmp %>% mutate(TP = ss == "ECOLI")

  return(tmp)
}


source("../../R/tidyMS_benchmark.R")
if(TRUE){
  allresultsStan <- readRDS("rstandTruncMixed.RDS")
  allresultsStanf <- allresultsStan %>% dplyr::filter(!sapply(summary, is.null))
  allresultsStanf <- allresultsStan %>% select(protein_Id, summary) %>% unnest(cols="summary")
  allresultsStan <- left_join(dplyr::select(allresultsStan,protein_Id), allresultsStanf)
  tmp <- preprocess(allresultsStan)
  tmp <- tmp %>% mutate(stanR.mean.estimate = abs(mean))
  tmpStanRProb <- ms_bench_add_FPRTPR(tmp,arrangeby = "minp", type="probability")
  tmpStanEEstimate <- ms_bench_add_FPRTPR(tmp, arrangeby = "stanR.mean.estimate", type="foldchange", desc=TRUE)
}


allresults <- readRDS("allresults.Rds")

tmp <- preprocess(allresults$ropeca_P)
tmp <- tmp %>% mutate(median.estimate = abs(median.estimate))
tmpRopecaBBS <- ms_bench_add_FPRTPR(tmp,arrangeby = "beta.based.significance", type="probability")

tmpRopecaMedianEstimate <- ms_bench_add_FPRTPR(tmp, arrangeby = "median.estimate", type="foldchange", desc=TRUE)

tmp <- preprocess(allresults$resXXmedpolish)
tmp <- tmp %>% mutate(pseudo_estimate = abs(pseudo_estimate))

tmpMedpolishPValue <- ms_bench_add_FPRTPR(tmp,arrangeby = "moderated.p.value.adjusted", type="probability")
tmpMedpolishEstimate <- ms_bench_add_FPRTPR(tmp, arrangeby = "pseudo_estimate", type="foldchange", desc=TRUE)
tmpMedpolishEstimate <- tmpMedpolishEstimate %>% mutate(arrangeby = paste0("medpolish.",arrangeby))

tmp <- preprocess(allresults$resXXmixmodel)
tmp <- tmp %>% mutate(pseudo_estimate = abs(pseudo_estimate))

tmpMMpValueMod <- ms_bench_add_FPRTPR(tmp, arrangeby = "moderated.p.value.adjusted",type="probability")
tmpMMpValueMod <- tmpMMpValueMod %>% mutate(arrangeby = paste0("mm.",arrangeby))

tmpMMpValue <- ms_bench_add_FPRTPR(tmp, arrangeby = "p.value.adjusted",type="probability")
tmpMMMedianEstimate <- ms_bench_add_FPRTPR(tmp, arrangeby = "pseudo_estimate",type="foldchange", desc=TRUE)
tmpMMMedianEstimate <- tmpMMMedianEstimate %>% mutate(arrangeby = paste0("mm.",arrangeby))

if(FALSE){
  xx <- tmpMMMedianEstimate %>%  dplyr::filter(contrast == "dilution_(4.5/3)_1.5")
  plot(xx$FPR, xx$TPR, xlim=c(0.694 , 0.6944))
  xx <- xx %>%  filter(FPR>0.694 &  FPR < 0.6944 )
  View(inner_join(tmp, xx))
}

res <- bind_rows(list(tmpStanRProb,tmpStanEEstimate,
                      tmpRopecaBBS,tmpRopecaMedianEstimate,
                      tmpMedpolishPValue,tmpMedpolishEstimate,tmpMMpValueMod,
                      tmpMMpValue,tmpMMMedianEstimate))

res <- bind_rows(list(tmpStanEEstimate,
                      tmpRopecaMedianEstimate,
                      tmpMedpolishEstimate,
                      tmpMMMedianEstimate))

res <- bind_rows(list(tmpStanRProb,
                      tmpRopecaBBS,
                      tmpMedpolishPValue,
                      tmpMMpValue))



res %>% group_by(arrangeby, type) %>% summarize(n = n())
res <- res %>% arrange(FDP)
res %>% filter(arrangeby == "minp" & contrast == "dilution_(9/6)_1.5"	) %>% View
foo
ggplot(res, aes(x=FDP,y = TPR, color = arrangeby)) +
  geom_line(aes(linetype=type)) +
  facet_wrap(~contrast) + ylim(0,1)


ggplot(res, aes(x = FPR, y = TPR, color = arrangeby)) +
  geom_line(aes(linetype=type)) +
  facet_wrap(~contrast) + ylim(0.5,1)
ggplot(res, aes(x = FPR, y = TPR, color = arrangeby)) +
  geom_line(aes(linetype=type)) +
  facet_wrap(~contrast) + ylim(0.5,1)
ggplot(res, aes(x = FPR, y = TPR, color = arrangeby)) +
  geom_line(aes(linetype=type)) +
  facet_wrap(~contrast) + xlim(0,0.2)

summary <- res %>%
  group_by(arrangeby,contrast) %>%
  summarise(n=n(), auc = auc(FPR, TPR), auc01 = auc(FPR, TPR, 0.2), type=unique(type) )
head(summary)
ggplot(summary, aes(x=contrast, y = auc , color= type)) +
  geom_bar(stat="identity") +
  facet_wrap(~arrangeby) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(summary, aes(x=arrangeby, y = auc01, color=type )) +
  geom_bar(stat="identity") +
  facet_wrap(~contrast) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
