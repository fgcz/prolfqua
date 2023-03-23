library(ggplot2)
library(prolfqua)

istar <- prolfqua_data('data_ionstar')
istar$config <- old2new(istar$config)
data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
lfqdata <- LFQData$new(data, istar$config)

sr <- lfqdata$get_Summariser()
#tmp <- sr$percentage_abundance()




plot_abundance_vs_percent <- function(sr,
                                      top_percent = 5,
                                      factors = TRUE ,
                                      log = FALSE,
                                      colors = c("^REV_" =  "red", "^CON_" = "orange"),
                                      cumulative = TRUE) {
  tmp <- sr$percentage_abundance()
  protID <- sr$lfq$config$table$hierarchy_keys_depth()
  if (!factors) {
    tmp <- tmp |> dplyr::filter(!!sym(sr$lfq$config$table$factor_keys_depth()[1]) == "ALL")
  }
  colorV <- rep("black", nrow(tmp))

  for (i in seq_along(colors)) {
    colorV[grepl(names(colors)[i], tmp[[protID]])] <- colors[i]
  }
  tmp$color <- colorV

  if (cumulative) {
    topN <- tmp |>
      dplyr::group_by(across(sr$lfq$config$table$factor_keys_depth())) |>
      dplyr::filter(abundance_percent_cumulative > top_percent)
  } else {
    topN <- tmp |>
      dplyr::group_by(across(sr$lfq$config$table$factor_keys_depth())) |>
      dplyr::filter(abundance_percent > top_percent)
  }
  ggplot <- ggplot(tmp, aes(x = !!sym("percent_prot"),
                            y = !!sym(if (cumulative) {"abundance_percent_cumulative"} else {"abundance_percent"}))) +
    geom_point(color = tmp$color) +
    facet_wrap(~dilution.) +
    ggrepel::geom_label_repel(data = topN,  aes(label = !!sym(protID)), size = 3) +
    if (log) {scale_y_continuous(trans = 'log10')} else {NULL}
  return(ggplot)
}

debug(plot_abundance_vs_percent)
plot_abundance_vs_percent(sr, top_percent = 20, factors = FALSE, cumulative = TRUE)
plot_abundance_vs_percent(sr, top_percent = 3, factors = FALSE, cumulative = FALSE)
plot_abundance_vs_percent(sr, top_percent = 20, factors = TRUE, cumulative = TRUE)
plot_abundance_vs_percent(sr, top_percent = 20, factors = TRUE, cumulative = FALSE)



