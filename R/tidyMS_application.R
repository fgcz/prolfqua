
#' visualize output of `contrasts_linfct``
#' @export
#' @keywords internal
#' @family deprecated
#'
contrasts_linfct_vis <- function(contrasts,
                                 modelName = "Model",
                                 prefix = "Contrasts",
                                 subject_Id = "protein_Id",
                                 columns = c("p.value","p.value.adjusted"),
                                 estimate = "estimate",
                                 contrast = "lhs",
                                 fc = 1){
  res <- list()
  contrasts %>% tidyr::unite("label", subject_Id, sep = "~", remove = FALSE) -> contrasts
  #return(contrasts)
  # add histogram of p-values
  for (column in columns) {
    fig <- list()
    name <- paste0(prefix,"_Histogram_",column)
    fig$fname <- paste0(name, "_", modelName )
    fig$fig <- ggplot(data = contrasts, aes(x = !!sym(column))) +
      geom_histogram(breaks = seq(0, 1, by = 0.05)) +
      facet_wrap(vars(!!sym(contrast)))
    res[[name]] <- fig
  }
  message("histograms created")
  # add volcano plots
  for (column in columns) {
    message(column)
    fig <- list()
    name <- paste0(prefix,"_Volcano_",column)
    fig$fname <- paste0(name, "_", modelName )
    fig$fig <- LFQService:::.multigroupVolcano(contrasts,
                                               effect = estimate,
                                               p.value = column,
                                               condition = contrast,
                                               text = "label",
                                               xintercept = c(-fc, fc),
                                               colour = "isSingular",
                                               scales = "free_y")

    message("volcano1")
    fig$plotly <- contrasts %>% plotly::highlight_key(~label) %>%
      LFQService:::.multigroupVolcano(.,
                                      effect = estimate,
                                      p.value = column,
                                      condition = contrast,
                                      text = "label",
                                      xintercept = c(-fc, fc),
                                      colour = "isSingular",
                                      scales = "free_y") %>%
      plotly::ggplotly(tooltip = "label")
    message("volcano plotly")
    res[[name]] <- fig
  }
  message("volcanos_build")
  # add histogram of fold changes
  {
    fig <- list()
    name <- paste0(prefix,"_Histogram_FC_estimate")
    fig$fname <- paste0(name, "_", modelName )
    fig$fig <- ggplot(data = contrasts, aes(x = !!sym(estimate))) +
      geom_histogram(breaks = seq(floor(min(contrasts[[estimate]], na.rm = TRUE)),
                                  ceiling(max(contrasts[[estimate]], na.rm = TRUE)), by = 0.1)) +
      facet_wrap(vars(!!sym(contrast)))
    res[[name]] <- fig
  }
  # MA plot
  {
    ma_plot <- function(x, fc = 1){
      x <- ggplot(x , aes(x = (c1 + c2)/2,
                          y = !!sym(estimate),
                          text = !!sym("label"),
                          colour = !!sym("isSingular"))) +
        geom_point(alpha = 0.5) + scale_colour_manual(values = c("black", "red")) +
        facet_wrap(vars(!!sym(contrast))) + theme_light() +
        geom_hline(yintercept = c(-fc, fc), linetype = "dashed",colour = "red")
      return(x)
    }

    if (!is.null(contrasts$c1) && !is.null(contrasts$c2)) {
      fig <- list()
      name <- paste0(prefix,"_MA_FC_estimate")
      fig$fname <- paste0(name, "_", modelName )
      # pdf version
      fig$fig <- contrasts %>% ma_plot(fc = fc)

      # html version
      fig$plotly  <- contrasts %>%
        plotly::highlight_key(~label) %>%
        ma_plot() %>%
        plotly::ggplotly(tooltip = "label")

      res[[name]] <- fig
    }
  }
  return(res)
}


#' helper function to write the result of `contrasts_linfct_vis`
#'
#' used in p2901
#'
#' @export
#' @keywords internal
#' @family deprecated
contrasts_linfct_vis_write <- function(fig_list,
                                       path,
                                       fig.width = 10,
                                       fig.height = 10,
                                       format = c("pdf","html")){
  format <- match.arg(format)
  if (!is.null(path)) {
    for (fig in fig_list) {

      fpath <- file.path(path,paste0(fig$fname,".", format))


      if (format == "pdf") {
        message("Writing: ",fpath,"\n")
        pdf(fpath, width = fig.width, height = fig.height)
        print(fig$fig)
        dev.off()
      }else if (format == "html") {
        if (!is.null(fig$plotly)) {
          message("Writing: ",fpath,"\n")
          htmlwidgets::saveWidget(widget = fig$plotly, fig$fname, selfcontained = TRUE)
          file.rename(fig$fname, fpath)
        }
      }
    }
  }
}


#' Do contrast
#' @export
#' @keywords internal
#' @family deprecated
#' @examples
#'
workflow_contrasts_linfct_V2 <- function(models,
                                         contrasts,
                                         config,
                                         modelName = "Model",
                                         prefix = "Contrasts",
                                         contrastfun = LFQService::my_contest )
{

  # extract contrast sides
  tt <- contrasts[grep("-",contrasts)]
  tt <- tibble(lhs = names(tt) , contrast = tt)
  tt <- tt %>% mutate(contrast = gsub("[` ]","",contrast)) %>%
    tidyr::separate(contrast, c("c1", "c2"), sep = "-")

  models <- models %>% dplyr::filter(exists_lmer == TRUE)
  m <- get_complete_model_fit(models)
  linfct <- linfct_from_model(m$linear_model[[1]], as_list = FALSE)
  linfct <- unique(linfct) # needed for single factor models
  linfct_A <- linfct_matrix_contrasts(linfct, contrasts)

  subject_Id <- config$table$hkeysDepth()
  contrast_result <- contrasts_linfct_deprec(models,
                                      rbind(linfct, linfct_A),
                                      subject_Id = subject_Id,
                                      contrastfun = contrastfun )

  xx <- contrast_result %>% dplyr::select(subject_Id, "lhs", "estimate")
  xx <- xx %>% pivot_wider(names_from = "lhs", values_from = "estimate")

  contrast_result <- contrast_result %>% dplyr::filter(lhs %in% names(contrasts))

  get_contrast_cols <- function(i, contrast_results , contrast_table , subject_ID ){
    data.frame(lhs = contrast_table[i, "lhs"],
               dplyr::select_at(contrast_results, c( subject_ID, unlist(contrast_table[i,c("c1","c2")]))),
               c1_name = contrast_table[i,"c1", drop = T],
               c2_name = contrast_table[i,"c2", drop = T], stringsAsFactors = FALSE)
  }

  contrast_sides <- purrr::map_df(1:nrow(tt), get_contrast_cols, xx, tt, subject_Id)
  contrast_result <- inner_join(contrast_sides,contrast_result)
  contrast_result <- moderated_p_limma_long(contrast_result)

  prefix <- prefix
  modelName <- modelName


  res_fun <- function(path = NULL, columns = c("p.value",
                                               "p.value.adjusted",
                                               "moderated.p.value",
                                               "moderated.p.value.adjusted"),
                      DEBUG = FALSE){
    if (DEBUG) {
      return(list(contrast_result = contrast_result,
                  linfct_A = linfct_A,
                  modelName = modelName,
                  config = config,
                  prefix = prefix,
                  subject_Id = subject_Id,
                  columns = columns
      ))
    }

    visualization <- contrasts_linfct_vis(contrast_result,
                                          modelName ,
                                          prefix = prefix,
                                          subject_Id = subject_Id,
                                          columns = columns
    )

    relevant_columns <- c("lhs",
                          "c1_name",
                          "c1",
                          "c2_name",
                          "c2",
                          "sigma",
                          "df",
                          "isSingular",
                          "estimate",
                          "conf.low",
                          "conf.high") # other relevant columns.

    contrast_minimal <- contrast_result %>%
      dplyr::select(subject_Id, relevant_columns, columns )
    contrasts_wide <- pivot_model_contrasts_2_Wide(contrast_minimal,
                                                   subject_Id = subject_Id,
                                                   columns = c("estimate", columns))

    if (!is.null(path)) {
      if (FALSE) {
        contrasts_linfct_write(contrast_minimal,
                               config,
                               path = path,
                               modelName = modelName,
                               prefix = prefix,
                               columns = c("estimate", columns))
      }

      contrasts_linfct_vis_write(visualization, path = path, format = "pdf")
      contrasts_linfct_vis_write(visualization, path = path, format = "html")
    }

    res <- list(contrast_result = contrast_result,
                contrast_minimal = contrast_minimal,
                contrasts_wide = contrasts_wide,
                visualization = visualization,
                modelName = modelName,
                linfct_A = linfct_A,
                prefix = prefix)

    invisible(res)
  }
  return( res_fun )
}


#' Writes figures generated by `model_analyse_summarize_vis`
#'
#' @export
#' @keywords internal
#'
#' @family deprecated
model_analyse_summarize_vis_write <- function(modelling_result,
                                              path,
                                              fig.width = 10 ,
                                              fig.height = 10,
                                              all = FALSE){
  if (all) {
    fpath <- file.path(path, modelling_result$fname_histogram_coeff_p.values)
    message("Writing figure into : ", fpath, "\n")
    pdf(fpath, width = fig.width, height = fig.height )
    print(modelling_result$histogram_coeff_p.values)
    dev.off()

    fpath <- file.path(path, modelling_result$fname_VolcanoPlot)
    message("Writing figure into : ", fpath, "\n")
    pdf(fpath,
        width = fig.width , height = fig.height)
    print(modelling_result$VolcanoPlot)
    dev.off()

    fpath <- file.path(path, modelling_result$fname_Pairsplot_Coef)
    message("Writing figure into : ", fpath, "\n")
    pdf(fpath, width = fig.width , height = fig.height)
    print(modelling_result$Pairsplot_Coef)
    dev.off()
  }
  fpath <- file.path(path, modelling_result$fname_histogram_anova_p.values)
  message("Writing figure into : ", fpath, "\n")
  pdf(fpath, width = fig.width , height = fig.height)
  print(modelling_result$histogram_anova_p.values)
  dev.off()
}


#' writes results of `model_analyse`, anova table and all the coefficients with parameters.
#' @keywords internal
#' @family deprecated
#' @export
#'
model_analyse_summarize_write  <- function(
  modellingResult,
  path,
  all = FALSE){
  message("writing tables into :", path)
  if (all) {
    lfq_write_table(modellingResult$Model_Coeff,
                    path = path,
                    name  = modellingResult$fname_Model_Coeff )
  }
  lfq_write_table(modellingResult$Model_Anova,
                  path = path ,
                  name = modellingResult$fname_Model_Anova )
}



#' visualize workflow model analyse results
#'
#' used in p2901
#'
#' @export
#' @keywords internal
#' @family deprecated
#' @examples
#'
#' D <- LFQServiceData::ionstar$normalized()
#' D$data <- dplyr::filter(D$data, protein_Id %in% sample(protein_Id, 10))
#' D$config$table$getWorkIntensity()
#' modelName <- "f_condtion_r_peptide"
#' formula_randomPeptide <-
#'   make_custom_model_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id)",
#'    model_name = modelName)
#' modellingResult <-  LFQService:::model_analyse(
#'  D$data,
#'  formula_randomPeptide,
#'  D$config$table$hkeysDepth(),modelName)
#' tmp <- model_analyse_summarize(modellingResult$modelProtein)
#' res <- model_analyse_summarize_vis(tmp,
#'  D$config$table$hkeysDepth())
#'
#' res$histogram_coeff_p.values
#' res$VolcanoPlot
#' res$Pairsplot_Coef
#' res$histogram_anova_p.values
#'
model_analyse_summarize_vis <- function(modellingResult,
                                        subject_Id ="protein_Id") {
  Model_Coeff <- tidyr::unite(modellingResult$Model_Coeff, "subject_Id", subject_Id)
  Model_Anova <- tidyr::unite(modellingResult$Model_Anova, "subject_Id", subject_Id)
  modelName <- modellingResult$modelName
  fig <- list()

  ## Coef_Histogram
  fig$fname_histogram_coeff_p.values <- paste0("Coef_Histogram_",modelName,".pdf")
  fig$histogram_coeff_p.values <- ggplot(data = Model_Coeff, aes(x = Pr...t.., group = row.names.x.)) +
    geom_histogram(breaks = seq(0,1,by = 0.05)) +
    facet_wrap(~row.names.x.)

  ## Coef_VolcanoPlot
  fig$fname_VolcanoPlot <- paste0("Coef_VolcanoPlot_",modelName,".pdf")
  fig$VolcanoPlot <- Model_Coeff %>%
    dplyr::filter(row.names.x. != "(Intercept)") %>%
    LFQService::multigroupVolcano(
      effect = "Estimate",
      p.value = "Pr...t..",
      condition = "row.names.x.",
      label = "subject_Id" ,
      xintercept = c(-1, 1) ,
      colour = "isSingular" )

  ## Coef_Pairsplot
  forPairs <- Model_Coeff %>%
    dplyr::select(!!sym("subject_Id") , row.names.x. ,  Estimate ) %>%
    tidyr::spread(row.names.x.,Estimate )
  fig$fname_Pairsplot_Coef <- paste0("Coef_Pairsplot_",modelName,".pdf")
  fig$Pairsplot_Coef <-  GGally::ggpairs(forPairs, columns = 2:ncol(forPairs))

  ## Anova_p.values
  fig$fname_histogram_anova_p.values <- paste0("Anova_p.values_", modelName, ".pdf")
  fig$histogram_anova_p.values <-  modellingResult$Model_Anova %>% dplyr::filter(rownames.x. != "Residuals") %>%
    ggplot( aes(x = Pr..F., group = rownames.x.)) +
    geom_histogram(breaks = seq(0,1,by = 0.05)) +
    facet_wrap(~rownames.x.)

  return(fig)
}

#' Compute fold changes given Contrasts
#' @keywords internal
#' @family deprecated
#' @export
#'
#' @examples
#'
#' library(LFQService)
#' library(tidyverse)
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#' configur$parameter$qVal_individual_threshold <- 0.01
#' data <- LFQService::removeLarge_Q_Values(data, configur)
#' data <- complete_cases(data, configur)
#'
#' Contrasts <- c("TimeT168vsT2" = "TimeT168 - TimeT2","TimeT168vsT24" = "TimeT168 - TimeT24" )
#' message("missigness_impute_factors_interactions : imputed")
#' xx <- missigness_impute_factors_interactions(data, configur, value = "nrMeasured" )
#' imputed <- missigness_impute_contrasts(xx, configur, Contrasts)
#'
#' xx <- missigness_impute_factors_interactions(data, configur, value = "imputed" )
#' imputed <- missigness_impute_contrasts(xx, configur, Contrasts)
#' xx <- missigness_impute_factors_interactions(data, configur, value = "meanArea" )
#' message("missigness_impute_factors_interactions : meanArea")
#' mean <- missigness_impute_contrasts(xx, configur, Contrasts)
#' head(mean)
missigness_impute_contrasts <- function(data,
                                        config,
                                        contrasts,
                                        agg_fun = function(x){median(x, na.rm = TRUE)})
{
  for (i in 1:length(contrasts)) {
    message(names(contrasts)[i], "=", contrasts[i],"\n")
    data <- dplyr::mutate(data, !!names(contrasts)[i] := !!rlang::parse_expr(contrasts[i]))
  }

  if (!is.null(agg_fun)) {
    data <- data %>% group_by_at(c("value" , config$table$hkeysDepth())) %>%
      summarise_if(is.numeric, agg_fun)
  }
  return(data)
}


#' Compute fold changes given Contrasts 2
#' @keywords internal
#' @family deprecated
#' @export
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#' configur$parameter$qVal_individual_threshold <- 0.01
#' data <- LFQService::removeLarge_Q_Values(data, configur)
#' data <- complete_cases(data, configur)
#' unique(data$Time)
#' Contrasts <- c("timeT24-T8" = "TimeT24-TimeT8", "timeT72-T8" = "TimeT72 - TimeT8")
#' res <- workflow_missigness_impute_contrasts(data, configur, Contrasts)
#' res("long", "factors")
#' bb <- res("long", "all")
#' bb
#' bb$contrast
#' xx <- res("wide", "all")
#' plot((xx$`imputed~TimeT24` - xx$`imputed~TimeT8`) -  xx$`imputed~timeT24-T8`)
#' range((xx$`meanArea~TimeT24` - xx$`meanArea~TimeT8`) -  xx$`meanArea~timeT24-T8`)
workflow_missigness_impute_contrasts <- function(data,
                                                 config,
                                                 contrasts, agg_fun = function(x){median(x, na.rm = TRUE)},
                                                 global = TRUE){

  xx <- missigness_impute_factors_interactions(data, config, value = "imputed", global = global )
  message("missigness_impute_factors_interactions : imputed")
  imputed <- missigness_impute_contrasts(xx, config, contrasts, agg_fun = agg_fun)
  xx <- missigness_impute_factors_interactions(data, config, value = "meanArea" , global = global)
  message("missigness_impute_factors_interactions : meanArea")
  mean <- missigness_impute_contrasts(xx, config, contrasts, agg_fun = agg_fun)

  dd <- dplyr::bind_rows(imputed, mean)
  dd_long <- dd %>% tidyr::gather("contrast","int_val",
                                  colnames(dd)[sapply(dd, is.numeric)])

  res_fun <- function(value = c("long", "wide","raw"),
                      what = c("contrasts", "factors", "all"),
                      DEBUG = FALSE){
    value <- match.arg( value )
    what  <- match.arg( what  )
    if (DEBUG) {
      return(list(value = value,
                  what = what,
                  dd_long = dd_long,
                  contrasts = contrasts,
                  config = config ))
    }

    if (what == "contrasts") {
      dd_long <- dplyr::filter(dd_long, contrast %in% names(contrasts))
    }else if (what == "factors") {
      dd_long <- dplyr::filter(dd_long, !contrast %in% names(contrasts))
    }else if (what == "all") {
    }

    if (value == "long") {
      long_xxxx <- dd_long %>% spread(value, int_val)
      return(long_xxxx)
    }else if (value == "wide") {
      dd <- dd_long %>% unite(contrast.v , value, contrast, sep = "~") %>% spread(contrast.v, int_val)
      xxx_imputed <- inner_join(LFQService::summarize_hierarchy(data,config),dd)
      return(xxx_imputed)
    }else if (value == "raw") {
      return(dd_long)
    }
  }
  return(res_fun)
}


.columnsImputed <- function(all, contrasts) {
  getAST <- function(ee) purrr::map_if(as.list(ee), is.call, getAST)

  get_sides <- function(contrast, all_variables) {
    ast_list <- getAST(rlang::parse_expr(contrast))
    ast_array <- array(as.character(unlist(ast_list)))
    bb <- intersect(gsub("`","",ast_array),all_variables)
    return(bb)
  }

  all_variables <- c(names(contrasts), unique(all$contrast))
  res <- NULL

  for (i in 1:length(contrasts)) {
    cname <- names(contrasts)[i]
    cc <- get_sides(contrasts[i], all_variables)
    if (length(cc) != 2) {
      message("there are ", length(cc) , "> 2 elements")
      next;
    }

    tmp <- all %>% dplyr::filter( .data$contrast %in% c(cname,cc) )
    tmp <- tmp %>% dplyr::select(-.data$meanArea) %>%
      tidyr::spread(.data$contrast , .data$imputed)

    tmp <- tmp %>% add_column(lhs = cname,.after = 1)
    tmp <- tmp %>% add_column(c1_name = cc[1],.after = 2)
    tmp <- tmp %>% add_column(c2_name = cc[2],.after = 3)
    tmp <- tmp %>% dplyr::rename(c1 = !!sym(cc[1]), c2 = !!sym(cc[2]), estimate = !!sym(cname))
    res <- dplyr::bind_rows(res,tmp)
  }
  return(res)
}


#' summarize - compute anova and extract model coefficients  generated by `model_analyse`
#'
#' @export
#' @keywords internal
#' @family deprecated
#' @examples
#' rm(list = ls())
#' library(LFQService)
#' D <- LFQServiceData::ionstar$normalized()
#' D$data <- dplyr::filter(D$data ,protein_Id %in% sample(protein_Id, 10))
#' modelName <- "f_condtion_r_peptide"
#' formula_randomPeptide <-
#'   make_custom_model_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id)",
#'    model_name = modelName)
#' pepIntensity <- D$data
#' config <- D$config
#' config$table$hkeysDepth()
#' modellingResult <- LFQService:::model_analyse(
#'  pepIntensity,
#'  formula_randomPeptide,
#'  modelName = modelName,
#'  subject_Id = config$table$hkeysDepth())
#' names(modellingResult)
#' tmp <- model_analyse_summarize(modellingResult$modelProtein,
#' subject_Id = config$table$hkeysDepth(),
#'  modelName = modelName)
#' head(tmp$Model_Coeff)
#' head(tmp$Model_Anova)
model_analyse_summarize <- function(modelProteinF,
                                    subject_Id = "protein_Id",
                                    modelName = "Model"
){
  lmermodel <- "linear_model"

  modelProteinF <- get_complete_model_fit(modelProteinF)
  # modelProteinF <- modelProteinF %>% dplyr::filter(nrcoef == max(nrcoef))

  # Extract coefficients
  .coef_df <-  function(x){
    x <- coef(summary(x));
    x <- data.frame(row.names(x), x);
    return(x)
  }

  Model_Coeff <- modelProteinF %>%
    dplyr::mutate(!!"Coeffs_model" := purrr::map( !!sym(lmermodel),  .coef_df ))

  Model_Coeff <- Model_Coeff %>%
    dplyr::select(!!!syms(subject_Id), !!sym("Coeffs_model"), isSingular, nrcoef)

  # Problem with new version of tidyr.
  if (FALSE) {
    Model_Coeff <-  tidyr::unnest(Model_Coeff, cols = "Coeffs_model")
  } else {
    Model_Coeff <- tidyr::unnest_legacy(Model_Coeff)
  }



  # ANOVA
  .anova_df <- function(x){
    x <- anova(x)
    colnames(x) <- make.names(colnames(x))
    x <- data.frame(rownames(x), x)
    return(x)
  }

  Model_Anova <- modelProteinF %>% dplyr::mutate(!!"Anova_model" := purrr::map( !!sym(lmermodel),  .anova_df ))

  Model_Anova <- Model_Anova %>%
    dplyr::select(!!!syms(subject_Id), !!sym("Anova_model"), isSingular, nrcoef)
  Model_Anova <- tidyr::unnest_legacy(Model_Anova)
  #tidyr::unnest(cols = "Anova_model")


  return(list(
    modelName = modelName,
    Model_Coeff = Model_Coeff,
    fname_Model_Coeff =  paste0("Coef_",modelName),
    Model_Anova = Model_Anova,
    fname_Model_Anova =  paste0("ANOVA_",modelName)
  ))
}


### Do missing value imputation
#' decorate the data
#' @export
#' @keywords internal
#' @family deprecated
#' @examples
#'
#' library(LFQService)
#' library(tidyverse)
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#'
#' configur$parameter$qVal_individual_threshold <- 0.01
#' data <- LFQService::removeLarge_Q_Values(data, configur)
#' data <- complete_cases(data, configur)
#' xx <- LFQService::transform_work_intensity(data, configur, log2)
#' head(xx)
#' configur$table$getWorkIntensity()
#' unique(data$Time)
#' Contrasts <- c("timeT24-T8" = "TimeT24-TimeT8", "timeT72-T8" = "TimeT72 - TimeT8")
#' res <- workflow_missigness_impute_contrasts_V2(xx, configur, Contrasts)
#' head(res)
#' plot((res$c1+res$c2)/2 ,(res$c1 - res$c2))
#' abline(h=0)
#' plot((res$c1 - res$c2) - res$estimate)
#'
#' plot(res$c1, res$c2, log="xy")
#' abline(0,1, col=2)
#'
workflow_missigness_impute_contrasts_V2 <- function(
  data,
  config,
  contrasts,
  do_not_report = "",
  agg_fun = function(x){median(x, na.rm = TRUE)},
  global = TRUE){
  res_contrasts_imputed <- workflow_missigness_impute_contrasts(data,
                                                                config,
                                                                contrasts, agg_fun = agg_fun,
                                                                global = global)
  contrasts_xx_imputed <- res_contrasts_imputed("long",what = "all")
  contrasts_xx_imputed <- .columnsImputed(contrasts_xx_imputed,
                                          contrasts = contrasts[setdiff(names(contrasts) ,
                                                                        do_not_report)])

}


#' workflow_model_analyse
#' @export
#' @references function with paramter path
#' @keywords internal
#' @family deprecated
#'
#' @examples
#' D <- LFQServiceData::ionstar$normalized()
#' D$data <- dplyr::filter(D$data, protein_Id %in% sample(protein_Id,12))
#' modelName <- "f_condtion_r_peptide"
#' formula_randomPeptide <-
#'   make_custom_model_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id)",
#'   model_name = modelName)
#' modellingResult <-  workflow_model_analyse(D$data,
#'  formula_randomPeptide,
#'   modelName,
#'  subject_Id = D$config$table$hkeysDepth())
#' reslist <- modellingResult()
workflow_model_analyse <- function(data,
                                   modelFunction,
                                   modelName = modelFunction$model_name,
                                   subject_Id = "protein_Id"){

  modellingResult <- model_analyse(data,
                                   modelFunction,
                                   modelName = modelName,
                                   subject_Id)

  # delay write
  res_fun <- function(path = NULL, all = FALSE, DEBUG = FALSE){
    if (DEBUG) {
      return(list(
        modellingResult = modellingResult,
        modelName = modelName,
        subject_Id = subject_Id
      ))
    }

    summaryResult <- model_analyse_summarize(modellingResult$modelProtein,
                                             modelName = modelName,
                                             subject_Id = subject_Id)
    visualization <- model_analyse_summarize_vis(summaryResult, subject_Id)

    if (!is.null(path)) {
      model_analyse_summarize_write(summaryResult, path, all = all)
      model_analyse_summarize_vis_write(visualization, path, all = all)
    }
    return(list(modellingResult = modellingResult,
                summaryResult = summaryResult,
                visualization = visualization))
  }
  return(res_fun)
}



# merges contrasts and imputed contrasts
.makeResult_contrasts <- function(contrast_minimal,
                                  contrasts_xx_imputed,
                                  subject_Id,
                                  config,
                                  remove_imputed = TRUE ) {

  contrast_results <- dplyr::right_join( contrast_minimal,
                                         contrasts_xx_imputed,
                                         by = c(subject_Id,
                                                "lhs", "c1_name", "c2_name"), suffix = c("","_imputed"))

  contrast_results <- dplyr::rename(contrast_results, contrast = lhs) #
  contrast_results <- contrast_results %>%
    dplyr::mutate(pseudo_estimate = dplyr::case_when(is.na(estimate) ~ estimate_imputed, TRUE ~ estimate))
  contrast_results <- contrast_results %>%
    dplyr::mutate(is_pseudo_estimate = dplyr::case_when(is.na(estimate) ~ TRUE, TRUE ~ FALSE))

  if (remove_imputed) {
    contrast_results <- contrast_results %>%
      dplyr::mutate(c1 = dplyr::case_when(is.na(estimate) ~ c1_imputed, TRUE ~ c1))
    contrast_results <- contrast_results %>%
      dplyr::mutate(c2 = dplyr::case_when(is.na(estimate) ~ c2_imputed, TRUE ~ c2))
    contrast_results <- contrast_results %>%
      dplyr::select(-dplyr::contains("_imputed"))
  }

  separate_hierarchy(contrast_results, config) -> filtered_dd
  return(filtered_dd)
}

#' run the modelling using lmer or lm models
#' @keywords internal
#' @family deprecated
#'
#' @param data data
#' @param config AnalysisConfiguration
#' @param modelFunction modelling function
#' @param contrasts contrasts
#' @param modelling_dir directory to store modelling results
#' @param remove_imputed todo
#' @param do_not_report contrasts not to report
#' @param DEBUG default FALSE
#' @export
#' @examples
#'
application_run_modelling_V2 <- function(data,
                                         config,
                                         modelFunction,
                                         contrasts,
                                         modelling_dir = "modelling_results_protein" ,
                                         remove_imputed = TRUE,
                                         do_not_report = "",
                                         DEBUG = FALSE)
{
  # create result structure
  modelling_path <- modelling_dir

  ### make modeling  -----
  modellingResult_fun <- workflow_model_analyse(data,
                                                modelFunction,
                                                subject_Id = config$table$hkeysDepth())


  modellingResult <-  modellingResult_fun()
  modelProteinF <- modellingResult$modellingResult$modelProtein
  message("finished modelling, starting contrasts estimation")
  res_contrasts <- workflow_contrasts_linfct_V2(modelProteinF,
                                                contrasts,
                                                config,
                                                modelName = modelFunction$model_name,
                                                prefix =  "contrasts",
                                                contrastfun = modelFunction$contrast_fun)
  message("finished contrast estimateion, starting missing value imputation and contrast estimation")
  ### Do missing value imputation
  contrasts_xx_imputed <- workflow_missigness_impute_contrasts_V2(data,
                                                                  config,
                                                                  contrasts, do_not_report = do_not_report)

  #### Compute contrasts from model ####
  message("returning results")

  # RESULT FUNCTION
  res_fun <- function(do = c("result",
                             "write_modelling",
                             "write_contrasts"),
                      DEBUG = FALSE,
                      remove_imputed = TRUE) {

    do <- match.arg(do)
    if (DEBUG) {
      res <- list(modelFunction = modelFunction,
                  imputed = contrasts_xx_imputed,
                  remove_imputed = remove_imputed,
                  subject_Id = config$table$hkeysDepth(),
                  modelling_path = modelling_path,
                  modellingResult_fun = modellingResult_fun,
                  res_contrasts = res_contrasts
      )
      return(res)
    }

    if (do == "result") {
      contrast_minimal <- res_contrasts(columns = modelFunction$report_columns)$contrast_minimal
      result_table <- .makeResult_contrasts(contrast_minimal
                                            ,contrasts_xx_imputed,
                                            config$table$hkeysDepth(),
                                            config,
                                            remove_imputed = remove_imputed)

      return(result_table)
    } else if (do == "write_modelling") {
      modellingResult_fun(modelling_path)
    } else if (do == "write_contrasts") {
      filtered_dd <- res_contrasts(modelling_path, columns = modelFunction$report_columns)
      result_table <- .makeResult_contrasts(filtered_dd$contrast_minimal
                                            ,contrasts_xx_imputed,
                                            config$table$hkeysDepth(),
                                            config)
      lfq_write_table(result_table,
                      path = modelling_path,
                      name = "foldchange_estimates")

      return(result_table)
    }
  }
  return(res_fun)
}




#' Add Annotation to a data.frame in long format
#' @family setup
#'
#' for an usage example see run_script lfq_mixed_model_inference
#' @param intensityData data imported using ``
#' @param inputAnnotation annotation
#' @param fileName column name to join on.
#' @export
#'
add_annotation <- function(intensityData,
                                       inputAnnotation,
                                       fileName = "raw.file") {
  ## read the data
  {# add annotation
    if ( is.character(inputAnnotation) ) {
      annotation <- readxl::read_xlsx(inputAnnotation)
    } else {
      annotation <- inputAnnotation
    }
    noData <- annotation[!annotation[[fileName]] %in% intensityData[[fileName]],]
    if (nrow(noData)) {
      message("some files in annotation have no measurements")
      message(paste(noData, collapse = " "))
    }
    measSamples <- unique(intensityData[[fileName]])
    noAnnot <- measSamples[!measSamples %in% annotation[[fileName]] ]
    if (length(noAnnot) > 0 ) {
      message("some measured samples have no annotation!")
      message(paste(noAnnot,collapse = " "))
    }
    resPepProtAnnot <- inner_join(annotation, intensityData, by = fileName)
    ###  Setup analysis ####
  }
  return(resPepProtAnnot)
}



#' Used for metabolomics data analysis.
#' @keywords internal
#'
#' @export
#' @examples
#' #todo
#'
application_summarize_compound <- function(data,
                                           config,
                                           qc_path,
                                           prefix = c("ms")) {
  qc_apth <- qc_path
  prefix <- match.arg(prefix)

  results <- LFQService:::normalize_log2_robscale(data, config)


  res_fun <- function(do = c("plot", "write", "render", "print_compounds", "data"),
                      DEBUG = FALSE){
    do <- match.arg(do)
    if (DEBUG) {
      return(list(qc_path = qc_path, prefix = prefix, results = results ))
    }

    if ( do == "plot") {
      quants_write(results$data, results$config, qc_path)
    }else if (do == "render") {
      LFQService::render_MQSummary_rmd(results$data,
                                       results$config$clone(deep = TRUE),
                                       pep = TRUE,
                                       workdir = ".",
                                       dest_path = qc_path,
                                       dest_file_name = paste0(prefix, "_intensities_qc"),
                                       format = "html")
    }else if (do == "print_compounds") {
      figs <- plot_hierarchies_line_df(results$data, results$config )
      pdf(file.path(qc_path, paste0(prefix, "_intensities_inference_figures.pdf")))
      lapply(figs, print)
      dev.off()
    }else if (do == "data") {
      return(list(data  = results$data, config = results$config))
    }

  }
  return(res_fun)

}

