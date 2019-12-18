rm(list=ls())
library(LFQServiceAnalysisTemplate)
library(tidyverse)
allresults <- readRDS("allresults.Rds")

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

tmp <- preprocess(allresults$ropeca_P)
tmp <- tmp %>% mutate(median.estimate = abs(median.estimate))
tmpRopecaBBS <- add_FPRTPR(tmp,arrangby = "beta.based.significance", type="probability")
tmpRopecaMedianEstimate <- add_FPRTPR(tmp, arrangby = "median.estimate", type="foldchange", desc=TRUE)

tmp <- preprocess(allresults$resXXmedpolish)
tmp <- tmp %>% mutate(pseudo_estimate = abs(pseudo_estimate))

tmpMedpolishPValue <- add_FPRTPR(tmp,arrangby = "moderated.p.value.adjusted", type="probability")
tmpMedpolishEstimate <- add_FPRTPR(tmp, arrangby = "pseudo_estimate", type="foldchange", desc=TRUE)
tmpMedpolishEstimate <- tmpMedpolishEstimate %>% mutate(arrangeby = paste0("medpolish.",arrangeby))

tmp <- preprocess(allresults$resXXmixmodel)
tmp <- tmp %>% mutate(pseudo_estimate = abs(pseudo_estimate))

tmpMMpValue <- add_FPRTPR(tmp, arrangby = "p.value.adjusted",type="probability")
tmpMMMedianEstimate <- add_FPRTPR(tmp, arrangby = "pseudo_estimate",type="foldchange", desc=TRUE)
tmpMMMedianEstimate <- tmpMMMedianEstimate %>% mutate(arrangeby = paste0("mm.",arrangeby))

if(FALSE){
  xx <- tmpMMMedianEstimate %>%  dplyr::filter(contrast == "dilution_(4.5/3)_1.5")
  plot(xx$FPR, xx$TPR, xlim=c(0.694 , 0.6944))
  xx <- xx %>%  filter(FPR>0.694 &  FPR < 0.6944 )
  View(inner_join(tmp, xx))
}

res <- bind_rows(list(tmpRopecaBBS,tmpRopecaMedianEstimate,
                      tmpMedpolishPValue,tmpMedpolishEstimate,
                      tmpMMpValue,tmpMMMedianEstimate))
ggplot(res, aes(x = FPR, y = TPR, color = arrangeby)) +
  geom_line(aes(linetype=type)) +
  facet_wrap(~contrast) + ylim(0.5,1)
ggplot(res, aes(x = FPR, y = TPR, color = arrangeby)) +
  geom_line(aes(linetype=type)) +
  facet_wrap(~contrast) + xlim(0,0.2)


summary <- res %>%
  group_by(arrangeby,contrast) %>%
  summarise(n=n(), auc = auc(FPR, TPR), auc01 = auc(FPR, TPR, 0.2) )
head(summary)
ggplot(summary, aes(x=contrast, y = auc )) +
  geom_bar(stat="identity") +
  facet_wrap(~arrangeby) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(summary, aes(x=arrangeby, y = auc01 )) +
  geom_bar(stat="identity") +
  facet_wrap(~contrast) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
