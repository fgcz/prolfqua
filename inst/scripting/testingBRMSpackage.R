rm(list = ls())
library(tidyverse)
options(mc.cores = parallel::detectCores())


summarised <- readRDS("aaa_summarized.RDA")
model <- "~ condition * batch + (1|peptide_Id)"

data <- summarised$results$pepIntensityNormalized
config <- summarised$results$config_pepIntensityNormalized
datacomp <- prolfqua::complete_cases(data, config)
data <- datacomp |> rename(condition = condition_, batch = batch_)


data_split <- data |>
  dplyr::group_by_at(config$table$hierarchy_keys()[1]) |>
  tidyr::nest()


data_split$data[[1]]

library("brms")

tmp <- data_split$data[[3]]
res <- vector(mode = "list", length = 5)
for (i in 1:5) {
  tmp_i <- tmp
  tmp_i$transformedIntensity[is.na(tmp_i$transformedIntensity)] <- rnorm(sum(is.na(tmp_i$transformedIntensity)), 0, 3)
  res[[i]] <- tmp_i
}

tmp |>
  group_by_at(c("condition", "batch")) |>
  summarize(nna = sum(is.na(transformedIntensity)))

fit1 <- brm_multiple(transformedIntensity ~ condition * batch + (1 | peptide_Id),
  data = res, iter = 4000
)
class(fit1)
summary(fit1)
marginal_effects(fit1)

postsample <- (posterior_samples(fit1, "^b"))
as.matrix(postsample) %*% as.matrix(c(1, 1, 0, 1))

colnames(postsample)


head(postsample)
plot(postsample$b_batchp2691, pch = ".")
dim(fitted(fit1))

fit1 <- brm(transformedIntensity | mi(sdy = 3) ~ condition * batch + (1 | peptide_Id),
  data = tmp, cores = 6, iter = 4000
)

tmp$transformedIntensity[is.na(tmp$transformedIntensity)] <- min(tmp$transformedIntensity, na.rm = TRUE)

fit1 <- brm(transformedIntensity | trunc() ~ condition * batch + (1 | peptide_Id),
  data = tmp, iter = 4000
)

summary(fit1)
marginal_effects(fit1)


fit1 <- brm(transformedIntensity | trunc(lb = -4) ~ condition * batch + (1 | peptide_Id),
  data = tmp, iter = 4000
)
summary(fit1)

marginal_effects(fit1)
xx <- posterior_samples(fit1, "^b")
head(xx)
plot(xx$b_Intercept)
plot(xx$b_batchp2691)

fit2 <- update(fit1, newdata = data_split$data[[4]])
marginal_effects(fit2)
fixef(fit2)
