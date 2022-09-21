library(tidyverse)
library(prolfqua)
istar <- prolfqua_data('data_ionstar')$normalized()
istar_data <- dplyr::filter(istar$data ,protein_Id %in% sample(protein_Id, 100))
pepIntensity <- istar_data
config <- istar$config$clone(deep = TRUE)

config$table$hierarchyDepth <- 2

ld <- LFQData$new(pepIntensity, config)
Contr <- c("dil.b_vs_a" = "dilution.a - dilution.b")

#debug(get_imputed_contrasts)

#debug(missigness_impute_factors_interactions)

get_imputed_contrasts_V2 <- function(pepIntensity,
                                     config,
                                     Contr,
                                     present = 0,
                                     global = TRUE){

  fun <- .missigness_impute_interactions(pepIntensity, config)
  long <- fun("long")
  x3 <- long |> filter(nrNAs == (max(long$nrNAs) - (present + 1))) |> pull(meanArea) |> mean(na.rm=TRUE)

  long <- long |> mutate(imputed_b = ifelse(is.na(meanArea), x3, meanArea))

  lt <- long
  imp <- lt |> pivot_wider(id_cols = config$table$hierarchyKeys(), names_from = interaction, values_from = imputed_b)
  lt <- lt |> mutate(nrNAs_b = ifelse( nrNAs == max(nrNAs) , nrNAs , 0) )
  nr <- lt |> pivot_wider(id_cols = config$table$hierarchyKeys(), names_from = interaction, values_from = nrNAs_b)

  imputed <- get_contrast(ungroup(imp), config$table$hierarchyKeys(), Contr)
  nrs <- get_contrast(ungroup(nr),  config$table$hierarchyKeys(), Contr)

  nrs <- nrs |> select(all_of(c(config$table$hierarchyKeys(),"contrast", "estimate" )))
  nrs <- nrs |> rename(indic = estimate)
  imputed <- inner_join(imputed, nrs)
  imputed2 <- imputed |> mutate(estimate = ifelse(indic < 0 & estimate < 0, 0, estimate))
  imputed2 <- imputed2 |> mutate(estimate = ifelse(indic > 0 & estimate > 0, 0, estimate))

  imputedProt <- aggregate_contrast(ungroup(imputed2),  subject_Id =  config$table$hkeysDepth())
  imputedProt$avgAbd <- (imputedProt$group_1 + imputedProt$group_2)/2
  imputedProt$group_1_name <- NULL
  imputedProt$group_2_name <- NULL
  imputedProt$group_1 <- NULL
  imputedProt$group_2 <- NULL
  return(imputedProt)
}




result = get_imputed_contrasts(
  ld$data,
  ld$config,
  Contr,
  global = TRUE)

result2 = get_imputed_contrasts_V2(
  ld$data,
  ld$config,
  Contr,
  global = TRUE)

dim(result)
dim(result2)
plot(result$estimate_median, result2$estimate_median, pch = "*")
abline(0,1, col=2)



result <- result2
result$isSingular <- TRUE
result <- select(result , -all_of(c("n","estimate_mad")))
var = summarize_stats(ld$data, ld$config)

pooled <- poolvar(var, ld$config, method = "V1")
pooled <- dplyr::select(pooled ,-all_of(c(ld$config$table$fkeysDepth()[1],"var")))
result <- dplyr::inner_join(result, pooled, by = ld$config$table$hkeysDepth())

resultNA <- result[result$n == 0, ]
resultnotNa <- result[result$n != 0,]
meandf <- resultnotNa |> summarize(n = 1, df = 1, sd = mean(sd),sdT = mean(sdT))
resultNA$n <- 0
resultNA$df <- 1
resultNA$sd <- meandf$sd
resultNA$sdT <- meandf$sdT
result <- bind_rows(resultNA, resultnotNa)

result <- dplyr::mutate(result, statistic = .data$estimate_median / .data$sdT,
                        p.value = 2*pt(abs(.data$statistic), df = .data$df, lower.tail = FALSE))

prqt <- -qt((1 - 0.95)/2, df = pmax(1,result$df))
result$conf.low <- result$estimate_median  - prqt * (result$sdT)
result$conf.high <- result$estimate_median + prqt * (result$sdT)




long2 <- tidyr::complete(long, tidyr::nesting(!!!syms(config$table$hierarchyKeys())), interaction)
dim(long2)
