resXXmedpolish()

mod <- LFQService:::build_model(
  protdata$data,
  modelFunction,
  modelName = modelFunction$modelName,
  subject_Id = config$table$hkeysDepth())
mod$write_anova_figures(path = ps$modelling_path())
mod$write_coef_figures(path = ps$modelling_path())

contrast <- LFQService::Contrasts$new(
  mod$modelDF,
  Contrasts,
  modelFunction$contrast_fun,
  subject_Id = config$table$hkeysDepth(),
  modelName = modelFunction$model_name)

contrast$get_contrasts_sides()
contrast$get_linfct()
contrast$get_contrasts()
head(contrast$moderate())

bb <- contrast$ropeca()
tmp <- contrast$get_contrasts()
cp <- Contrasts_Plotter$new(tmp , contrast$subject_Id)
x <- cp$histogram()
x <- cp$histogram_estimate()
x <- cp$volcano()
x <- cp$volcano_plotly()
x <- cp$ma_plot()
x <- cp$ma_plotly()

cp$write_pdf(path = ps$modelling_path())
cp$write_plotly(path = ps$modelling_path())
head(cp$to_wide())
#source("c:/Users/wolski/prog/LFQService/R/tidyMS_application.R")
