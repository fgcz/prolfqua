library(ggplot2)
library(prolfqua)

istar <- prolfqua_data('data_ionstar')
istar$config <- old2new(istar$config)
data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
lfqdata <- LFQData$new(data, istar$config)

sr <- lfqdata$get_Summariser()
tmp <- sr$percentage_abundance()

tmp |> dplyr::group_by(dilution.) |> dplyr::summarize(sum = sum(percentAbundance))

head(tmp)
if (TRUE) {
  tmp <- tmp |> filter(dilution. == "ALL")
}
topN <- tmp |>
  dplyr::group_by(across(lfqdata$config$table$factor_keys_depth())) |>
  dplyr::top_n(5, wt = meanArea)

ggplot(tmp, aes(x = percent_prot, y = abundance_percent_cumulative)) +
  geom_point() +
  facet_wrap(~dilution.) +
  geom_text(data = topN,  aes(label = !!sym(colnames(bb$hierarchy()))), size = 3) #+

dplyr::is_grouped_df(tmp)
tmp2 <- tmp |> dplyr::filter(dilution. == "ALL")
approx(tmp2$percentProt, tmp2$cumulativePercAbundance, xout = 90)


tmp |> dplyr::group_by(dilution.) |> dplyr::summarize(
  percentSignal_for_least_Abundant_50 = approx(percentProt, cumulativePercAbundance, xout = 50)$y,
  percentSignal_for_mostAbundant_10 = 100 - approx(percentProt, cumulativePercAbundance, xout = 90)$y,
  percentSignal_for_mostAbundant_1 = 100 - approx(percentProt, cumulativePercAbundance, xout = 99)$y,
  )


