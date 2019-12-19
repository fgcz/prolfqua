# Reminder of missigness analysis.
da <- summarised$results$pepIntensityNormalized
cf <- summarised$results$config_pepIntensityNormalized
p <- missigness_histogram(da, cf)

xx3 <- missingness_per_condition(data_cc, local_config)
xx4 <- missingness_per_condition_cumsum(data_cc, local_config)
missingPrec
ims <- interaction_missing_stats(da, cf)

tmp <- ims %>% filter(nrReplicates-nrNAs == 1)
xx <- inner_join( da , tmp )
hist(na.omit(xx$transformedIntensity))

tmp <- ims %>% filter(nrReplicates-nrNAs == 2)
xx <- inner_join( da , tmp )
xx %>% na.omit() -> xx
#xx %>% select(config$table$idVars(), config$table$fkeysLevel(), config$table$getWorkIntensity()) %>% tidyr::spread(c(config$table$idVars(), config$table$fkeysLevel()), config$table$getWorkIntensity())

xxx <- xx %>% group_by_at(c(config$table$hierarchyKeys(), config$table$fkeysLevel())) %>% mutate(index = row_number()) %>% ungroup() #%>%

xxx <- xxx %>% dplyr::select(c(config$table$hierarchyKeys(), config$table$fkeysLevel(), config$table$getWorkIntensity(),"index"))
xxx <- xxx %>% tidyr::spread("index","transformedIntensity")
plot(xxx$`1`,xxx$`2`)


tmp <- ims %>% filter(nrReplicates-nrNAs == 3)
xx <- inner_join( da , tmp )
hist(na.omit(xx$transformedIntensity))

