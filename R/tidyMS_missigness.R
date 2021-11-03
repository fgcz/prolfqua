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
#' library(prolfqua)
#'
#' bb <- prolfqua::data_ionstar$filtered()
#' stopifnot(nrow(bb$data) == 25780)
#' stopifnot(nrow(bb$data) == 25780)
#' #saveRDS(bb, file="c:/users/wolski/__debugR/aaaaa.rds")
#' configur <- bb$config
#' data <- bb$data
#'
#' configur$parameter$qVal_individual_threshold <- 0.01
#' xx <- prolfqua::removeLarge_Q_Values(data,
#'    configur)
#' xx <- complete_cases(xx, configur)
#' x <- interaction_missing_stats(xx, configur)$data %>% arrange(desc(nrNAs))
#'
#' #readr::write_tsv(x, file="c:/users/wolski/__debugR/aaaaaa.tsv")
#' stopifnot(nrow(x) == 5540)
#' stopifnot(sum(is.na(x$meanArea)) == 206)
#' stopifnot(length(unique(x$protein_Id)) == 162)
#'
#' tmp <- interaction_missing_stats(xx, configur,
#'  factors= character(),
#'   hierarchy = configur$table$hierarchyKeys()[1])$data
#' stopifnot(nrow(tmp) == 162)
#'
#' tmp <- interaction_missing_stats(xx, configur,
#'   hierarchy = configur$table$hierarchyKeys()[1])$data
#' stopifnot(sum(is.na(tmp$nrMeasured))==0)
#' tmp
#'
#' interaction_missing_stats(xx, configur, factors = NULL)
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
                     meanArea = mean(!!sym(workIntensity), na.rm = TRUE),
                     medianArea = median(!!sym(workIntensity), na.rm = TRUE)) %>%
    mutate(nrMeasured = .data$nrReplicates - .data$nrNAs) %>% dplyr::ungroup()
  return(list(data = missingPrec,
              summaries = c("nrReplicates","nrNAs","nrMeasured","meanArea", "medianArea")))
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
#'
#' library(prolfqua)
#' bb <- prolfqua::data_ionstar$filtered()
#' stopifnot(nrow(bb$data) == 25780)
#' configur <- bb$config
#' data <- bb$data
#' configur$parameter$qVal_individual_threshold <- 0.01
#'
#' xx <- prolfqua::removeLarge_Q_Values(data, configur)
#' xx <- complete_cases(xx, configur)
#'
#' tmp <- interaction_missing_stats(xx, configur)
#' fun <- .missigness_impute_interactions(xx, configur)
#'
#' long <- fun("long")
#' head(long)
#' debugData <- fun(DEBUG=TRUE)
#' names(debugData)
#' head(debugData$long)
#' sum(is.na(xx$long$nrReplicates))
#' xxx <- (fun("nrReplicates"))
#' alldata <- fun("all")
#' head(xxx)
#'
#' imputed <- fun("imputed")
#' missing <- fun("nrMeasured")
#'
#'  meanArea <- fun("mean")
#'  print(sum(is.na(meanArea$mean.dilution.a)))
#'  #stopifnot(sum(is.na(meanArea$mean.dilution.a)) == 59)
#'  stopifnot(sum(is.na(imputed$mean.imp.dilution.a))==0)
#'
.missigness_impute_interactions <- function(pdata,
                                            config,
                                            factors = config$table$fkeysDepth(),
                                            probs = 0.1,
                                            global = TRUE){
  mstats <- interaction_missing_stats(pdata, config, factors = factors)
  x_summaries <- mstats$summaries
  mstats <- mstats$data
  mstats <- make_interaction_column(mstats, factors, sep = ":")


  lowerMean <- function(meanArea, probs = probs){
    meanAreaNotNA <- na.omit(meanArea)
    small10 <- meanAreaNotNA[meanAreaNotNA < quantile(meanAreaNotNA, probs = probs)]
    meanArea[is.na(meanArea)] <- mean(small10)
    return(meanArea)
  }

  if (!global) {
    mstats <- mstats %>%
      group_by(interaction) %>%
      mutate(imputed = lowerMean(.data$meanArea,probs = probs))
  }else{
    mstats <- mstats %>%
      mutate(imputed = lowerMean(.data$meanArea,probs = probs))

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
      return(list(value = value, long = mstats , config = config ))
    }

    if (value == "long") {
      return(mstats)
    }else{
      mstats <- mstats %>% dplyr::select(-one_of(factors))

      pid <- config$table$hkeysDepth()
      nrReplicates <- mstats %>%
        dplyr::select( -one_of(c(setdiff(x_summaries,"nrReplicates"),"imputed") )) %>%
        tidyr::spread(interaction, nrReplicates, sep = ".nrReplicates.") %>%
        arrange(!!!syms(pid)) %>%
        dplyr::ungroup()
      nrMeasured <- mstats %>% dplyr::select(-one_of(c(setdiff(x_summaries,"nrMeasured"),"imputed" ) )) %>%
        tidyr::spread(interaction, nrMeasured, sep = ".nrMeasured.") %>%
        arrange(!!!syms(pid)) %>% dplyr::ungroup()

      meanArea <- mstats %>% dplyr::select(-one_of(c(setdiff(x_summaries,"meanArea"),"imputed" ) )) %>%
        tidyr::spread(interaction, meanArea, sep = ".meanArea.") %>%
        arrange(!!!syms(pid)) %>% dplyr::ungroup()

      meanAreaImputed <- mstats %>% dplyr::select(-one_of(setdiff(x_summaries,"imputed" ) )) %>%
        tidyr::spread(interaction, .data$imputed, sep = ".imputed.") %>%
        arrange(!!!syms(pid)) %>% dplyr::ungroup()

      allTables <- list(meanArea = meanArea,
                        nrMeasured = nrMeasured,
                        nrReplicates = nrReplicates,
                        meanAreaImputed = meanAreaImputed)

      if (value == "all") {
        allTables[["long"]] <- mstats
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
#' @family imputation
#' @examples
#'
#' bb <- prolfqua::data_ionstar$filtered()
#' stopifnot(nrow(bb$data) == 25780)
#' configur <- bb$config
#' data <- bb$data
#' xx <- complete_cases(data, configur)
#'
#' res <- missigness_impute_factors_interactions(xx, configur)
#' head(res)
#' res <- missigness_impute_factors_interactions(xx, configur, value = "imputed")
#' head(res)
#' fun <- missigness_impute_factors_interactions(xx, configur, value = "nrMeasured")
#' fun
missigness_impute_factors_interactions <-
  function(pdata,
           config,
           probs = 0.03,
           value = c("nrReplicates", "nrMeasured", "meanArea", "imputed"),
           add.prefix = FALSE,
           global = TRUE)
  {
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

    fac_res <- vector(mode = "list", length = length(fac_fun))
    names(fac_res) <- names(fac_fun)
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
#' @family imputation
#' @export
#'
#' @examples
#'
#' library(prolfqua)
#' library(tidyverse)
#' bb <- prolfqua::data_ionstar$normalized()
#' configur <- bb$config
#' data <- bb$data
#'
#' Contrasts <- c("dilution.b-a" = "dilution.b - dilution.a", "dilution.c-e" = "dilution.c - dilution.e")
#' mean <- missigness_impute_factors_interactions(data, configur, value = "meanArea" )
#' mean <- get_contrast(mean, configur$table$hierarchyKeys(), Contrasts)
#' meanProt <- aggregate_contrast(mean,  subject_Id =  configur$table$hkeysDepth())
#'
#' imputed <- missigness_impute_factors_interactions(data, configur, value = "imputed" )
#' imputed <- get_contrast(imputed, configur$table$hierarchyKeys(), Contrasts)
#' head(imputed)
#'
#' imputedProt <- aggregate_contrast(imputed,  subject_Id =  configur$table$hkeysDepth())
#' plot(imputedProt$c1 - imputedProt$c2, imputedProt$estimate_median)
#' abline(c(0,1), col=2, pch = "*")
#' dim(meanProt)
#' sum(is.na(meanProt$estimate_median)) == 0
#' sum(is.na(imputedProt$estimate_median)) == 0
#' plot(meanProt$estimate_median - imputedProt$estimate_median )
#'
aggregate_contrast <- function(
  data,
  subject_Id ,
  agg_func = list(median = function(x){ stats::median(x, na.rm = TRUE) },
                   mad = function(x){ stats::mad(x, na.rm = TRUE)} ),
  contrast = "contrast")
{
  grouping_columns <- c(contrast, subject_Id, "c1_name","c2_name")
  dataG <- data %>%
    group_by(!!!syms(grouping_columns))

  resN <- dataG %>% dplyr::summarise(n = n(), .groups="drop")
  resE <- dataG %>% dplyr::summarise(across(.data$estimate,
                                     agg_func
                                    ), .groups = "drop")
  agg_func_c <- agg_func[1]
  resC <- dataG %>% dplyr::summarise(across(all_of(c("c1", "c2")),
                                     agg_func_c,
                                     .names = "{col}"
                                     ),
                              .groups = "drop")
  #res <- dplyr::full_join(ungroup(resC),ungroup(resE), by = grouping_columns)
  res <- Reduce(function(x,y){ dplyr::full_join(x , y , by = grouping_columns)},
                list(resN,resE,resC) )
  return(res)

}
#' Compute fold changes given Contrasts
#' @keywords internal
#' @family imputation
#' @param data data.frame
#' @param data hierarchyKeys of Analysis Configuration
#' @param contrasts list of contrasts to compute
#' @export
#' @examples
#'
#'
#' library(prolfqua)
#' library(tidyverse)
#' bb <- prolfqua::data_ionstar$filtered()
#' stopifnot(nrow(bb$data) == 25780)
#' configur <- bb$config
#' data <- bb$data
#' data <- complete_cases(data, configur)
#'
#' Contrasts <- c("aVSe" = "dilution.a - dilution.e","aVSb" = "dilution.a - dilution.b" )
#' message("missigness_impute_factors_interactions : imputed")
#' xx <- missigness_impute_factors_interactions(data, configur, value = "nrMeasured" )
#' imputed <- get_contrast(xx, configur$table$hierarchyKeys(), Contrasts)
#'
#' xx <- missigness_impute_factors_interactions(data, configur, value = "imputed" )
#'
#' imputed <- get_contrast(xx, configur$table$hierarchyKeys(), Contrasts)
#' head(imputed)
#'
get_contrast <- function(data,
                         hierarchyKeys,
                         contrasts)
{
  getAST <- function(ee) purrr::map_if(as.list(ee), is.call, getAST)
  get_sides <- function(contrast, all_variables) {
    ast_list <- getAST(rlang::parse_expr(contrast))
    ast_array <- array(as.character(unlist(ast_list)))
    bb <- intersect(gsub("`","",ast_array),all_variables)
    return(bb)
  }


  for (i in 1:length(contrasts)) {
    message(names(contrasts)[i], "=", contrasts[i],"\n")
    data <- dplyr::mutate(data, !!names(contrasts)[i] := !!rlang::parse_expr(contrasts[i]))
  }
  res <- vector(mode = "list", length(contrasts))
  names(res) <- names(contrasts)
  for (i in 1:length(contrasts)) {
    sides <- get_sides(contrasts[i], colnames(data))
    df  <- select(data , c( hierarchyKeys, c1 = sides[1], c2 = sides[2], estimate = names(contrasts)[i]))
    df$c1_name <- sides[1]
    df$c2_name <- sides[2]
    df$contrast <-  names(contrasts)[i]

    res[[names(contrasts)[i]]] <- df
  }
  res <- dplyr::bind_rows(res)

  return(ungroup(res))
}

#' compute contrasts based on peptide fold changes
#'
#' Computes median of fold peptide fold change to obtain protein fold change.
#' peptide fold change is computed based on the average of intensities in the conditions.
#'
#'
#' @param data data.frame
#' @param config AnalysisConfiguration
#' @param contrast list of contrasts
#' @keywords internal
#' @family missingness
#' @family modelling
#' @export
#' @examples
#'
#'
#' library(prolfqua)
#' library(tidyverse)
#' bb <- prolfqua::data_ionstar$normalized()
#' configur <- bb$config
#' data <- bb$data
#' configur$parameter$qVal_individual_threshold <- 0.01
#' data <- prolfqua::removeLarge_Q_Values(data, configur)
#' data <- complete_cases(data, configur)
#'
#' Contrasts <- c("dilution.b-a" = "dilution.b - dilution.a", "dilution.c-e" = "dilution.c - dilution.e")
#' res <- get_imputed_contrasts(data, configur, Contrasts)
#' head(res)
#'
#' if(FALSE){
#' debug(get_imputed_contrasts)
#' config <- configur
#' contrasts <- Contrasts
#' imputed <- missigness_impute_factors_interactions(data, config, value = "imputed" )
#' head(imputed)
#' imputed <- get_contrast(ungroup(imputed), config$table$hierarchyKeys(), contrasts)
#' head(imputed)
#' imputedProt <- aggregate_contrast(imputed,  subject_Id =  config$table$hkeysDepth())
#' head(imputedProt)
#' }
get_imputed_contrasts <- function(data, config, contrasts, probs = 0.03, global = TRUE){
  imputed <- missigness_impute_factors_interactions(data, config, value = "imputed" ,probs = probs, global = global)
  imputed <- get_contrast(ungroup(imputed), config$table$hierarchyKeys(), contrasts)
  imputedProt <- aggregate_contrast(ungroup(imputed),  subject_Id =  config$table$hkeysDepth())
  return(imputedProt)
}
#' Histogram summarizing missigness
#' @export
#' @keywords internal
#' @family plotting
#' @family imputation
#' @examples
#' library(tidyverse)
#' library(prolfqua)
#' bb <- prolfqua::data_ionstar
#' configur <- bb$config
#' data <- bb$data
#' xx <- complete_cases(data, configur)
#' missigness_histogram(xx, configur)
#'
#' missingPrec <- interaction_missing_stats(xx, configur)
#'
#' bx <- prolfqua::data_ionstar$normalized()
#' configur <- bx$config
#' data <- bx$data
#' data <- complete_cases(data, configur)
#' missingPrecNorm <- interaction_missing_stats(data, configur)
#'
#' missigness_histogram(data, configur, showempty=FALSE)
#' missigness_histogram(data, configur, showempty=TRUE)
missigness_histogram <- function(x,
                                 config,
                                 showempty = FALSE,
                                 factors = config$table$fkeysDepth()){
  table <- config$table
  missingPrec <- interaction_missing_stats(x, config , factors)$data
  missingPrec <- missingPrec %>%
    dplyr::ungroup() %>% dplyr::mutate(nrNAs = as.factor(.data$nrNAs))

  if (showempty) {
    if (config$table$is_intensity_transformed) {
      missingPrec <- missingPrec %>%
        dplyr::mutate(meanArea = ifelse(is.na(.data$meanArea), min(.data$meanArea, na.rm = TRUE) - 1,
                                        .data$meanArea))
    }else{
      missingPrec <- missingPrec %>%
        dplyr::mutate(meanArea = ifelse(is.na(.data$meanArea),min(.data$meanArea, na.rm = TRUE) - 20,.data$meanArea))
    }

  }

  factors <- table$fkeysDepth()
  formula <- paste(table$isotopeLabel, "~", paste(factors, collapse = "+"))
  message(formula)

  p <- ggplot(missingPrec, aes(x = .data$meanArea, fill = .data$nrNAs, colour = .data$nrNAs)) +
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
#' @family imputation
#' @examples
#'
#'
#' bb <- prolfqua::data_ionstar$filtered()
#' stopifnot(nrow(bb$data) == 25780)
#' configur <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
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
    dplyr::mutate(cumulative_sum = cumsum(.data$nrTransitions))
  res <- xxcs  %>% dplyr::select(-.data$nrTransitions)

  formula <- paste(table$isotopeLabel, "~", paste(factors, collapse = "+"))
  message(formula)

  nudgeval = mean(res$cumulative_sum) * 0.05
  p <- ggplot(res, aes(x = .data$nrNAs, y = .data$cumulative_sum)) +
    geom_bar(stat = "identity", color = "black", fill = "white") +
    geom_text(aes(label = .data$cumulative_sum), nudge_y = nudgeval, angle = -45) +
    facet_grid(as.formula(formula))

  res <- res %>% tidyr::spread("nrNAs","cumulative_sum")
  return(list(data = res, figure = p))
}

#' Summarize missing in condtion as barplot
#' @export
#' @keywords internal
#' @family plotting
#' @family imputation
#' @examples
#' bb <- prolfqua::data_ionstar$filtered()
#' stopifnot(nrow(bb$data) == 25780)
#' configur <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' res <- missingness_per_condition(data, configur)
#' stopifnot(c(5,8) == dim(res$data))
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




