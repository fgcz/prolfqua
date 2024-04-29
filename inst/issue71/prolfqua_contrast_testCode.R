library(prolfqua)
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")

#### data ####
mp_df <- read.csv("protein_intensities.csv")
long_df <- data.frame(mp_df) %>%
  pivot_longer(cols = -protein_Id, names_to = "raw.file",
               values_to = "medpolish", values_drop_na = TRUE)

annot <- read.csv("prolfqua-inputAnnotation.csv")

#### configuration ####
startdata <- dplyr::inner_join(long_df, annot, by = "raw.file")
atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("protein_Id")
atable$hierarchyDepth <- 1
atable$set_response("medpolish")
atable$factors[["condition_"]] = "condition"
atable$factors[["replicate"]] = "replicate"
atable$factorDepth <- 1

config <- prolfqua::AnalysisConfiguration$new(atable)
adata <- prolfqua::setup_analysis(startdata, config)
lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$hierarchy_counts()
lfqdata$remove_small_intensities()
lfqdata$factors()

#### contrasts ####
formula_Condition <-  strategy_lm("medpolish ~ condition_")

Contrasts <- c("Psfoam - Psfilm" = "condition_Psfoam - condition_Psfilm",
               "Psfoam - glucose" = "condition_Psfoam - condition_glucose",
               "Psfilm - glucose" = "condition_Psfilm - condition_glucose")
lfqdata <- lfqdata$complete_cases()

mod <- prolfqua::build_model(
  lfqdata$data,
  formula_Condition,
  subject_Id = lfqdata$config$table$hierarchy_keys())


contr <- prolfqua::Contrasts$new(mod, Contrasts, global = FALSE)

contr$get_linfct(global = FALSE)

contrdf <- contr$get_contrasts()

plotter <- contr$get_Plotter()
v1 <- plotter$volcano()

mod$modelDF
mod$modelDF$linear_model[1]

modg99273 <- mod$modelDF |> filter(protein_Id == "g99273")
modg99273$linear_model[1]


# check protein g99273
lfqdata$to_wide()$data %>% filter(protein_Id == "g99273")
v1$FDR$data %>% filter(protein_Id == "g99273")


