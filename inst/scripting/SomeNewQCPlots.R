tmp <- tmp$data |> arrange(meanArea)
tmp$id <- 1:nrow(tmp)
tmp$txt <- paste0(tmp$protein_Id, tmp$peptide_Id)
top25 <- tmp |>
  group_by(across(lfqdata$config$table$factor_keys_depth())) |>
  top_n(5, wt = meanArea)
dim(top25)

top100 <- tmp |>
  group_by(across(lfqdata$config$table$factor_keys_depth())) |>
  top_n(100, wt = meanArea)
ggplot(tmp , aes(x = meanArea, y = meanArea)) +
  geom_point() +
  facet_wrap(~dilution.) +
  geom_text(data = top25,  aes(label =  txt), size = 3) #+
#  scale_y_continuous(trans = log2_trans())


