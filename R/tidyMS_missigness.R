# Functions - Missigness ----

#' compute missingness statistics per hierarchy and factor level
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param factors factor to include (default up to factor depth)
#' @param hierarchy hierarchy to include (default up to hierarchy depth)
#' @param workIntensity work intensity column
#' @export
#' @keywords internal
#' @examples
#' library(tidyverse)
#' library(LFQService)
#'
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#'
#' configur$parameter$qVal_individual_threshold <- 0.01
#' xx <- LFQService::removeLarge_Q_Values(data,
#'    configur)
#' xx <- complete_cases(xx, configur)
#' nrow(xx)
#' x <- interaction_missing_stats(xx, configur)$data %>% arrange(desc(nrNAs))
#' stopifnot(nrow(x) == 7416)
#' stopifnot(sum(is.na(x$meanArea)) == 249)
#' stopifnot(length(unique(x$protein_Id)) == 37)
#'
#' tmp <- interaction_missing_stats(xx, configur,
#'  factors= character(),
#'   hierarchy = configur$table$hierarchyKeys()[1])$data
#' stopifnot(nrow(tmp) == 37)
#'
#' tmp <- interaction_missing_stats(xx, configur,
#'   hierarchy = configur$table$hierarchyKeys()[1])$data
#' stopifnot(sum(is.na(tmp$nrMeasured))==0)
#'
interaction_missing_stats <- function(pdata,
                                      config,
                                      factors = config$table$fkeysDepth(),
                                      hierarchy = config$table$hierarchyKeys(),
                                      workIntensity = config$table$getWorkIntensity())
{
  pdata <- complete_cases(pdata, config)
  table <- config$table
  missingPrec <- pdata %>% group_by_at(c(factors,
                                         hierarchy,
                                         table$isotopeLabel
  ))
  missingPrec <- missingPrec %>%
    dplyr::summarize(nrReplicates = n(),
                     nrNAs = sum(is.na(!!sym(workIntensity))),
                     meanArea = mean(!!sym(workIntensity), na.rm = TRUE)) %>%
    mutate(nrMeasured = nrReplicates - nrNAs) %>% dplyr::ungroup()
  return(list(data = missingPrec,
              summaries = c("nrReplicates","nrNAs","nrMeasured","meanArea")))
}

#' Compute interaction averages and
#' impute data using mean of lowest 0.1 (default)
#'
#' used in Acetylation project p2916
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param factors factor to include (default up to factor depth)
#' @param probs quantile to take average from (default 0.1)
#' @param global global min value
#' @return function with parameter `value`
#' `c("long", "nrReplicates", "nrMeasured", "meanArea", "imputed", "allWide", "all")`
#' @export
#' @keywords internal
#' @return function
#' @examples
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#' configur$parameter$qVal_individual_threshold <- 0.01
#'
#' xx <- LFQService::removeLarge_Q_Values(data, configur)
#' xx <- complete_cases(xx, configur)
#'
#' tmp <- interaction_missing_stats(xx, configur)
#' fun <- .missigness_impute_interactions(xx, configur)
#'
#' dd <- fun("long")
#' head(dd)
#' xx <- fun(DEBUG=TRUE)
#' names(xx)
#' xx$long
#' sum(is.na(xx$long$nrReplicates))
#' xxx <- (fun("nrReplicates"))
#' head(xxx)
#' xxx <- fun("all")
#' xi <- fun("imputed")
#' xi
#' head(xxx)
#'
.missigness_impute_interactions <- function(pdata,
                                            config,
                                            factors = config$table$fkeysDepth(),
                                            probs = 0.1,
                                            global = TRUE){
  mstats <- interaction_missing_stats(pdata, config, factors = factors)
  x_summaries <- mstats$summaries
  xx <- mstats$data
  xx <- make_interaction_column(xx, factors, sep = ":")


  lowerMean <- function(meanArea, probs = probs){
    meanAreaNotNA <- na.omit(meanArea)
    small10 <- meanAreaNotNA[meanAreaNotNA < quantile(meanAreaNotNA, probs = probs)]
    meanArea[is.na(meanArea)] <- mean(small10)
    return(meanArea)
  }

  if (!global) {
    xx <- xx %>%
      group_by(interaction) %>%
      mutate(imputed = lowerMean(.data$meanArea,probs = 0.2))
  }else{
    xx <- xx %>%
      mutate(imputed = lowerMean(.data$meanArea,probs = 0.2))

  }
  res_fun <- function(value = c("long",
                                "nrReplicates",
                                "nrMeasured",
                                "meanArea",
                                "imputed",
                                "allWide",
                                "all" ),
                      add.prefix = TRUE,
                      DEBUG = FALSE){
    value <- match.arg(value)
    if (DEBUG) {
      return(list(value = value, long = xx , config = config ))
    }

    if (value == "long") {
      return(xx)
    }else{
      xx <- xx %>% dplyr::select(-one_of(factors))

      pid <- config$table$hkeysDepth()
      nrReplicates <- xx %>%
        dplyr::select( -one_of(c(setdiff(x_summaries,"nrReplicates"),"imputed") )) %>%
        tidyr::spread(interaction, nrReplicates, sep = ".nrReplicates.") %>%
        arrange(!!!syms(pid)) %>%
        dplyr::ungroup()
      nrMeasured <- xx %>% dplyr::select(-one_of(c(setdiff(x_summaries,"nrMeasured"),"imputed" ) )) %>%
        tidyr::spread(interaction, nrMeasured, sep = ".nrMeasured.") %>%
        arrange(!!!syms(pid)) %>% dplyr::ungroup()

      meanArea <- xx %>% dplyr::select(-one_of(c(setdiff(x_summaries,"meanArea"),"imputed" ) )) %>%
        tidyr::spread(interaction, meanArea, sep = ".meanArea.") %>%
        arrange(!!!syms(pid)) %>% dplyr::ungroup()

      meanAreaImputed <- xx %>% dplyr::select(-one_of(setdiff(x_summaries,"imputed" ) )) %>%
        tidyr::spread(interaction, imputed, sep = ".imputed.") %>%
        arrange(!!!syms(pid)) %>% dplyr::ungroup()

      allTables <- list(meanArea = meanArea,
                        nrMeasured = nrMeasured,
                        nrReplicates = nrReplicates,
                        meanAreaImputed = meanAreaImputed)

      if (value == "all") {
        allTables[["long"]] <- xx
        return(allTables)
      }else if (value == "allWide") {
        return(purrr::reduce(allTables, inner_join))
      }else if (value == "nrReplicates") {
        srepl <- if (add.prefix) {"nrRep."}else{""}
        colnames(nrReplicates) <- gsub("interaction.nrReplicates.", srepl ,colnames(nrReplicates))
        nrReplicates <- tibble::add_column( nrReplicates, "value" = value, .before = 1)
        return(nrReplicates)
      }else if (value == "nrMeasured") {
        srepl <- if (add.prefix) {"nrMeas."}else{""}
        colnames(nrMeasured) <- gsub("interaction.nrMeasured.", srepl ,colnames(nrMeasured))
        nrMeasured <- tibble::add_column( nrMeasured, "value" = value, .before = 1)
        return(nrMeasured)
      }else if (value == "meanArea") {
        srepl <- if (add.prefix) {"mean."}else{""}
        colnames(meanArea) <- gsub("interaction.meanArea.", srepl ,colnames(meanArea))
        meanArea <- tibble::add_column( meanArea, "value" = value, .before = 1)
        return(meanArea)
      }else if (value == "imputed") {
        srepl <- if (add.prefix) {"mean.imp."}else{""}
        colnames(meanAreaImputed) <- gsub("interaction.imputed.", srepl ,colnames(meanAreaImputed))
        meanAreaImputed <- tibble::add_column( meanAreaImputed, "value" = value, .before = 1)
        return(meanAreaImputed)
      }
    }
  }

  #  nrMeasured %>% dplyr::select(starts_with("interaction")) -> nrMeasuredM
  #  nrReplicates %>% dplyr::select(starts_with("interaction")) -> nrReplicatesM
  return(res_fun)
}


#' compute per group averages and impute values
#' should generalize at some stage
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param probs quantile to take average from (default 0.1)
#' @param value use default
#' @param add.prefix use default
#'
#' @export
#' @keywords internal
#' @examples
#'
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#' configur$parameter$qVal_individual_threshold <- 0.01
#' xx <- LFQService::removeLarge_Q_Values(data, configur)
#' xx <- complete_cases(xx, configur)
#' res <- missigness_impute_factors_interactions(xx, configur)
#' head(res)
#' res <- missigness_impute_factors_interactions(xx, configur, value = "imputed")
#' head(fun)
#' dim(fun)
#' dim(dplyr::distinct(fun[,1:6]))
#' fun <- missigness_impute_factors_interactions(xx, configur, value = "nrMeasured")
missigness_impute_factors_interactions <-
  function(pdata,
           config,
           probs = 0.1,
           value = c("nrReplicates", "nrMeasured", "meanArea", "imputed"),
           add.prefix = FALSE,
           global = TRUE){
    value <- match.arg(value)
    fac_fun <- list()
    fac_fun[["interaction"]] <- .missigness_impute_interactions(
      pdata,
      config,
      probs = probs,
      global = global)
    if (config$table$factorDepth > 1 ) { # if 1 only then done
      for (factor in config$table$fkeysDepth()) {
        fac_fun[[factor]] <- .missigness_impute_interactions(
          pdata,
          config,
          factors = factor,
          probs = probs,
          global = global)
      }
    }

    fac_res <- list()
    for (fun_name in names(fac_fun)) {
      fac_res[[fun_name]] <- fac_fun[[fun_name]](value, add.prefix = add.prefix)
    }

    intfact <- purrr::reduce(fac_res,
                             dplyr::inner_join,
                             by = c(config$table$hierarchyKeys(),
                                    config$table$isotopeLabel, "value"))
    return(intfact)
  }



#' Compute fold changes given Contrasts
#' @keywords internal
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

  xx <- missigness_impute_factors_interactions(data, config, value = "imputed", global=global )
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



### Do missing value imputation
#' decorate the data
#' @export
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


#' Histogram summarizing missigness
#' @export
#' @keywords internal
#' @family plotting
#' @examples
#' library(tidyverse)
#' library(LFQService)
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#'
#' xx <- complete_cases(data, configur)
#' configur$parameter$qVal_individual_threshold <- 0.01
#' xx <- LFQService::removeLarge_Q_Values(data, configur)
#' xx <- complete_cases(xx, configur)
#' missigness_histogram(xx, configur)
#'
#' missingPrec <- interaction_missing_stats(xx, configur)
#'
#' setNa <- function(x){ifelse(x < 100, NA, x)}
#' data %>% dplyr::mutate(Area = setNa(Area)) -> data
#' missigness_histogram(data, configur)
#'
missigness_histogram <- function(x, config, showempty = TRUE, factors = config$table$fkeysDepth()){
  table <- config$table
  missingPrec <- interaction_missing_stats(x, config , factors)$data
  missingPrec <- missingPrec %>%  dplyr::ungroup() %>% dplyr::mutate(nrNAs = as.factor(nrNAs))

  if (showempty) {
    if (config$table$is_intensity_transformed) {
      missingPrec <- missingPrec %>% dplyr::mutate(meanArea = ifelse(is.na(meanArea),1,meanArea))
    }else{
      missingPrec <- missingPrec %>% dplyr::mutate(meanArea = ifelse(is.na(meanArea),-20,meanArea))
    }

  }

  factors <- table$fkeysDepth()
  formula <- paste(table$isotopeLabel, "~", paste(factors, collapse = "+"))
  message(formula)

  p <- ggplot(missingPrec, aes(x = meanArea, fill = nrNAs, colour = nrNAs)) +
    geom_histogram(alpha = 0.2, position = "identity") +
    facet_grid(as.formula(formula)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  if (!config$table$is_intensity_transformed) {
    p <- p + scale_x_log10()
  }
  p
}

#' cumulative sums of missing
#' @export
#' @keywords internal
#' @family plotting
#' @examples
#'
#' setNa <- function(x){ifelse(x < 100, NA, x)}
#'
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#'
#' data %>% dplyr::mutate(Area = setNa(Area)) -> data
#' res <- missingness_per_condition_cumsum(data,configur)
#' stopifnot("ggplot" %in% class(res$figure))
#' print(res$figure)
#' res$data
missingness_per_condition_cumsum <- function(x,
                                             config,
                                             factors = config$table$fkeysDepth()){
  table <- config$table
  missingPrec <- interaction_missing_stats(x, config,factors)$data

  xx <- missingPrec %>% group_by_at(c(table$isotopeLabel, factors,"nrNAs","nrReplicates")) %>%
    dplyr::summarize(nrTransitions = n())

  xxcs <- xx %>% group_by_at( c(table$isotopeLabel,factors)) %>% arrange(nrNAs) %>%
    dplyr::mutate(cumulative_sum = cumsum(nrTransitions))
  res <- xxcs  %>% dplyr::select(-nrTransitions)

  formula <- paste(table$isotopeLabel, "~", paste(factors, collapse = "+"))
  message(formula)

  nudgeval = mean(res$cumulative_sum) * 0.05
  p <- ggplot(res, aes(x = nrNAs, y = cumulative_sum)) +
    geom_bar(stat = "identity", color = "black", fill = "white") +
    geom_text(aes(label = cumulative_sum), nudge_y = nudgeval, angle = -45) +
    facet_grid(as.formula(formula))

  res <- res %>% tidyr::spread("nrNAs","cumulative_sum")
  return(list(data = res, figure = p))
}

#' Summarize missing in condtion as barplot
#' @export
#' @keywords internal
#' @family plotting
#' @examples
#' setNa <- function(x){ifelse(x < 100, NA, x)}
#' bb <- LFQServiceData::skylinePRMSampleData_A
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#' data %>% dplyr::mutate(Area = setNa(Area)) -> data
#' res <- missingness_per_condition(data, configur)
#' names(res)
#' stopifnot(c(8,7) == dim(res$data))
#' stopifnot("ggplot" %in% class(res$figure))
#' print(res$figure)
#'
missingness_per_condition <- function(x, config, factors = config$table$fkeysDepth()){
  table <- config$table
  missingPrec <- interaction_missing_stats(x, config, factors)$data
  hierarchyKey <- tail(config$table$hierarchyKeys(),1)
  hierarchyKey <- paste0("nr_",hierarchyKey)
  xx <- missingPrec %>% group_by_at(c(table$isotopeLabel,
                                      factors,"nrNAs","nrReplicates")) %>%
    dplyr::summarize( !!sym(hierarchyKey) := n())

  formula <- paste(table$isotopeLabel, "~", paste(factors, collapse = "+"))
  #message(formula)

  nudgeval = max(xx[[hierarchyKey]]) * 0.05

  p <- ggplot(xx, aes_string(x = "nrNAs", y = hierarchyKey)) +
    geom_bar(stat = "identity", color = "black", fill = "white") +
    geom_text(aes(label = !!sym(hierarchyKey)), nudge_y = nudgeval, angle = 45) +
    facet_grid(as.formula(formula))
  xx <- xx %>% tidyr::spread("nrNAs",hierarchyKey)

  return(list(data = xx ,figure = p))
}




