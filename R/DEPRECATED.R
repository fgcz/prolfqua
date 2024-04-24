get_impute_contrasts_V1 <- function(
    lfqdata, contrasts,
    method = "V1",
    minds = 1 ,
    present = 1,
    confint = 0.95,
    .p.adjust = prolfqua::adjust_p_values,
    all = FALSE) {

  result = get_imputed_contrasts(
    lfqdata$data,
    lfqdata$config,
    contrasts,
    present = present)

  # compute statistics using pooled variance
  result$isSingular <- TRUE
  result <- select(result , -all_of(c("nr_children","estimate_mad")))

  var = summarize_stats(lfqdata$data, lfqdata$config)

  pooled <- poolvar(var, lfqdata$config, method = method)
  pooled <- dplyr::select(pooled ,-all_of(c(lfqdata$config$table$factor_keys_depth()[1],"var")))

  result <- dplyr::inner_join(result, pooled, by = lfqdata$config$table$hierarchy_keys_depth())

  result_sd_zero <- result[result$df == 0, ]
  resultnot_zero <- result[result$df > 0,]
  meandf <- resultnot_zero |> summarize(
    n = 1, df = 1,
    sd = quantile(sd, prob = 0.75, na.rm = TRUE),
    sdT = quantile(sdT, prob = 0.75, na.rm = TRUE))

  meandf$sd <-  ifelse(meandf$sd > 0, meandf$sd, minsd)
  meandf$sdT <-  ifelse(meandf$sdT > 0, meandf$sdT, minsd)

  result_sd_zero$df <- 1
  result_sd_zero$sd <- meandf$sd
  result_sd_zero$sdT <- meandf$sdT
  result <- bind_rows(result_sd_zero, resultnot_zero)

  result <- result |> mutate(sd = ifelse(sd > 0 , sd, meandf$sd))
  result <- result |> mutate(sdT = ifelse(sdT > 0 , sdT, meandf$sdT))

  result <- dplyr::mutate(result, statistic = .data$estimate_median / .data$sdT,
                          p.value = 2*pt(abs(.data$statistic), df = .data$df, lower.tail = FALSE))

  prqt <- -qt((1 - confint)/2, df = result$df)
  result$conf.low <- result$estimate_median  - prqt * (result$sdT)
  result$conf.high <- result$estimate_median + prqt * (result$sdT)
  result <- .p.adjust(result, column = "p.value", group_by_col = "contrast", newname = "FDR")

  if (!all) {
    result <- select(result, -all_of( c("isSingular", "nrMeasured" , "mean" ,"n.groups", "n", "meanAll") ) )
  }
  return(result)
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
#' @param present estimate lod from group with default 1 present values.
#' @keywords internal
#' @family missingness
#' @family modelling
#' @return data.frame with name of contrast, ID, number of observations (if with aggregation), median, mad, and avgAbd
#' @export
#' @examples
#'
#'
#' istar <- sim_lfq_data_peptide_config()
#' config <- istar$config
#' analysis <- istar$data
#' data <- complete_cases(analysis, config)
#'
#' Contrasts <- c("dilution.b-a" = "group_A - group_B", "dilution.c-e" = "group_A - group_Ctrl")
#' #debug(get_imputed_contrasts)
#' res <- get_imputed_contrasts(data, config, Contrasts)
#' config <- config
#' contrasts <- Contrasts
#' imputed <- missigness_impute_factors_interactions(data, config, value = "imputed" )
#' imputed <- get_contrast(imputed, config$table$hierarchy_keys(), contrasts)
#' imputedProt <- aggregate_contrast(imputed,  subject_Id =  config$table$hierarchy_keys_depth())
#'
get_imputed_contrasts <- function(pepIntensity,
                                  config,
                                  Contr,
                                  present = 1){
  message("!!!!deprecated!!!!")
  if (!present > 0) {
    stop("At least 1 observation in interaction to infer LOD.")
  }
  long <- missigness_impute_factors_interactions(pepIntensity, config, value = "long" )
  # determine limit of detection
  LOD <- long |> filter(nrNAs == nrReplicates - present) |> pull(meanAbundance) |> median(na.rm=TRUE)

  long <- tidyr::complete(long, tidyr::nesting(!!!syms(config$table$hierarchy_keys())), interaction)
  long <- long |> mutate(imputed_b = ifelse(is.na(meanAbundance), LOD, meanAbundance))

  lt <- long
  imp <- lt |> pivot_wider(id_cols = config$table$hierarchy_keys(), names_from = interaction, values_from = imputed_b)
  lt <- lt |> mutate(is_missing = ifelse( nrNAs == nrReplicates , 1 , 0) )
  nr <- lt |> pivot_wider(id_cols = config$table$hierarchy_keys(), names_from = interaction, values_from = is_missing)

  imputed <- get_contrast(ungroup(imp), config$table$hierarchy_keys(), Contr)
  nrs <- get_contrast(ungroup(nr),  config$table$hierarchy_keys(), Contr)

  nrs <- nrs |> select(all_of(c(config$table$hierarchy_keys(),"contrast", "estimate" )))
  nrs <- nrs |> rename(indic = estimate)
  imputed <- inner_join(imputed, nrs)
  imputed2 <- imputed |> mutate(estimate = ifelse(indic < 0 & estimate < 0, 0, estimate))
  imputed2 <- imputed2 |> mutate(estimate = ifelse(indic > 0 & estimate > 0, 0, estimate))

  imputedProt <- aggregate_contrast(ungroup(imputed2),  subject_Id =  config$table$hierarchy_keys_depth())
  imputedProt$avgAbd <- (imputedProt$group_1 + imputedProt$group_2)/2
  imputedProt$group_1_name <- NULL
  imputedProt$group_2_name <- NULL
  imputedProt$group_1 <- NULL
  imputedProt$group_2 <- NULL
  return(imputedProt)
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
#' istar <- sim_lfq_data_peptide_config(weight_missing = 2)
#' config <- istar$config
#' analysis <- istar$data
#'
#' xx <- complete_cases(analysis, config)
#' res <- missigness_impute_factors_interactions(xx, config)
#' res <- missigness_impute_factors_interactions(xx, config, value = "imputed")
#' res <- missigness_impute_factors_interactions(xx, config, value = "nrMeasured")
#' #debug(missigness_impute_factors_interactions)
#' long <- missigness_impute_factors_interactions(xx, config, value = "long")
#' head(long)
#' plot(long$meanAbundance, long$imputed)
missigness_impute_factors_interactions <- function(
    pdata,
    config,
    probs = 0.03,
    value = c("long", "nrReplicates", "nrMeasured", "meanAbundance", "imputed"),
    add.prefix = FALSE,
    global = TRUE)
{
  message("!!!! deprecated !!!!")
  value <- match.arg(value)
  fac_fun <- list()
  fac_fun[["interaction"]] <- .missigness_impute_interactions(
    pdata,
    config,
    probs = probs,
    global = global)
  if (config$table$factorDepth > 1 ) { # if 1 only then done
    for (factor in config$table$factor_keys_depth()) {
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
  if (value == "long") {
    intfact <- dplyr::bind_rows(fac_res)
  } else {
    intfact <- purrr::reduce(fac_res,
                             dplyr::inner_join,
                             by = c(config$table$hierarchy_keys(),
                                    config$table$isotopeLabel, "value"))

  }
  return(dplyr::ungroup(intfact))
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
#' `c("long", "nrReplicates", "nrMeasured", "meanAbundance", "imputed", "allWide", "all")`
#' @export
#' @keywords internal
#' @return function
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config(Nprot = 20,weight_missing = 2)
#' config <- istar$config
#' analysis <- istar$data
#' xx <- complete_cases(analysis, config)
#' nrPepTimesDilution <- length(unique(paste0(xx$protein_Id, xx$peptide_Id))) *
#'     length(unique(xx$group_))
#' funx <- .missigness_impute_interactions(xx, config)
#' long <- funx("long")
#' alldata <- funx("all")
#' stopifnot(length(names(alldata)) == 5)
#'
#' imputed <- funx("imputed")
#' stopifnot(nrow(imputed) == length(unique(paste0(xx$protein_Id, xx$peptide_Id))))
#' missing <- funx("nrMeasured")
#' stopifnot(nrow(missing) == length(unique(paste0(xx$protein_Id, xx$peptide_Id))))
#'
#'  meanAbundance <- funx("mean")
#' stopifnot(nrow(meanAbundance) == length(unique(paste0(xx$protein_Id, xx$peptide_Id))))
#'  stopifnot(sum(is.na(imputed$mean.imp.group_A))==0)
#'
.missigness_impute_interactions <- function(pdata,
                                            config,
                                            factors = config$table$factor_keys_depth(),
                                            probs = 0.1,
                                            global = TRUE){
  message("deprecated")
  mstats <- interaction_missing_stats(pdata, config, factors = factors)
  x_summaries <- mstats$summaries
  mstats <- mstats$data
  mstats <- make_interaction_column(mstats, factors, sep = ":")

  lowerMean <- function(meanAbundance, probs = probs){
    meanAbundanceNotNA <- na.omit(meanAbundance)
    small10 <- meanAbundanceNotNA[meanAbundanceNotNA < quantile(meanAbundanceNotNA, probs = probs)]
    meanAbundance[is.na(meanAbundance)] <- mean(small10)
    return(meanAbundance)
  }

  if (!global) {
    mstats <- mstats |>
      group_by(interaction) |>
      dplyr::mutate(imputed = lowerMean(.data$meanAbundance,probs = probs))
  }else{
    mstats <- mstats |>
      dplyr::mutate(imputed = lowerMean(.data$meanAbundance,probs = probs))
  }

  res_fun <- function(value = c("long",
                                "nrReplicates",
                                "nrMeasured",
                                "meanAbundance",
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
      mstats <- mstats |> dplyr::select(-one_of(factors))

      pid <- config$table$hierarchy_keys_depth()
      nrReplicates <- mstats |>
        dplyr::select( -one_of(c(base::setdiff(x_summaries,"nrReplicates"),"imputed") )) |>
        tidyr::spread(interaction, nrReplicates, sep = ".nrReplicates.") |>
        dplyr::arrange(!!!syms(pid)) |>
        dplyr::ungroup()
      nrMeasured <- mstats |> dplyr::select(-one_of(c(base::setdiff(x_summaries,"nrMeasured"),"imputed" ) )) |>
        tidyr::spread(interaction, nrMeasured, sep = ".nrMeasured.") |>
        dplyr::arrange(!!!syms(pid)) |> dplyr::ungroup()

      meanAbundance <- mstats |> dplyr::select(-one_of(c(base::setdiff(x_summaries,"meanAbundance"),"imputed" ) )) |>
        tidyr::spread(interaction, meanAbundance, sep = ".meanAbundance.") |>
        dplyr::arrange(!!!syms(pid)) |> dplyr::ungroup()

      meanAbundanceImputed <- mstats |> dplyr::select(-one_of(base::setdiff(x_summaries,"imputed" ) )) |>
        tidyr::spread(interaction, .data$imputed, sep = ".imputed.") |>
        dplyr::arrange(!!!syms(pid)) |> dplyr::ungroup()

      allTables <- list(meanAbundance = meanAbundance,
                        nrMeasured = nrMeasured,
                        nrReplicates = nrReplicates,
                        meanAbundanceImputed = meanAbundanceImputed)

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
      }else if (value == "meanAbundance") {
        srepl <- if (add.prefix) {"mean."}else{""}
        colnames(meanAbundance) <- gsub("interaction.meanAbundance.", srepl ,colnames(meanAbundance))
        meanAbundance <- tibble::add_column( meanAbundance, "value" = value, .before = 1)
        return(meanAbundance)
      }else if (value == "imputed") {
        srepl <- if (add.prefix) {"mean.imp."}else{""}
        colnames(meanAbundanceImputed) <- gsub("interaction.imputed.", srepl ,colnames(meanAbundanceImputed))
        meanAbundanceImputed <- tibble::add_column( meanAbundanceImputed, "value" = value, .before = 1)
        return(meanAbundanceImputed)
      }
    }
  }

  #  nrMeasured |> dplyr::select(starts_with("interaction")) -> nrMeasuredM
  #  nrReplicates |> dplyr::select(starts_with("interaction")) -> nrReplicatesM
  return(res_fun)
}

