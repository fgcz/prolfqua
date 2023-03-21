library(prolfqua)
istar <- prolfqua_data('data_ionstar')
istar$config <- old2new(istar$config)
data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
lfqdata <- LFQData$new(data, istar$config)


# roll up to protein intensities
ag <- lfqdata$get_Aggregator()
lfqdata$hierarchy_counts()
bb <- ag$sum_topN(N = 10)

bb$rename_response("totalIntensity")

# compute protein level summaries
dall <- interaction_missing_stats(bb$data, bb$config, factors = NULL)
dfac <- interaction_missing_stats(bb$data, bb$config)
xd <- setdiff(colnames(dfac$data), colnames(dall$data))

for (i in xd) {
  dall$data[[i]] <- "ALL"
}

all <- bind_rows(dfac$data, dall$data)

nested <- all |> group_by(!!sym(lfqdata$config$table$factor_keys_depth())) |> tidyr::nest()

for (i in seq_len(nrow(nested))){
  nested$data[[i]] <- nested$data[[i]] |>
    arrange(.data$meanArea) |>
    mutate(id = row_number()) |>
    mutate(percentAbundance = meanArea/sum(meanArea)*100 ) |>
    mutate(cumulativePercAbundance = cumsum(percentAbundance)) |>
    mutate(percentProt = id / max(id) * 100)
}

tmp <- tidyr::unnest(nested, cols = "data")
tmp |> group_by(dilution.) |> summarize(sum = sum(percentAbundance))


topN <- tmp |>
  group_by(across(lfqdata$config$table$factor_keys_depth())) |>
  top_n(5, wt = meanArea)

if (FALSE) {
  ggplot(tmp , aes(x = percentProt, y = percentAbundance)) +
    geom_point() +
    facet_wrap(~dilution.) +
    geom_text(data = topN,  aes(label = !!sym(colnames(bb$hierarchy()))), size = 3) #+
  #scale_y_continuous(trans = log2_trans())

}


ggplot(tmp, aes(x = percentProt, y = cumulativePercAbundance, colour = dilution.)) +
  geom_point() +
  facet_wrap(~dilution.) +
  geom_text(data = topN,  aes(label = !!sym(colnames(bb$hierarchy()))), size = 3) #+

is_grouped_df(tmp)
tmp2 <- tmp |> filter(dilution. == "ALL")
approx(tmp2$percentProt, tmp2$cumulativePercAbundance, xout = 90)

tmp |> group_by(dilution.) |> summarize(
  yout99 = approx(percentProt, cumulativePercAbundance, xout = 99)$y,
  yout90 = approx(percentProt, cumulativePercAbundance, xout = 90)$y,
  yout50 = approx(percentProt, cumulativePercAbundance, xout = 50)$y,
  )

