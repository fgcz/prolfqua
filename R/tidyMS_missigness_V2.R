#' compute group mean by LOD
#'
#' weight lod by nr of NA's $(LOD * nrNas + meanAbundance *nrObs)/(nrMeasured)$
#'
#' @export
#' @examples
#' Contrasts <- c("group.b-a" = "group_A - group_B", "group.a-ctrl" = "group_A - group_Ctrl")
#' dd <- prolfqua::sim_lfq_data_protein_config(Nprot = 100,weight_missing = 2)
#' mh <- prolfqua::MissingHelpers$new(dd$data, dd$config, prob = 0.8,weighted = TRUE)
#' mh$get_stats()
#' mh$get_LOD()
#' mh$impute_weighted_lod()
#' mh$impute_lod()
#' mh$get_poolvar()
#' bb <- mh$get_contrast_estimates(Contrasts)
#' mh$get_contrasts(Contrasts)
#'
#' dd <- prolfqua::sim_lfq_data_protein_2Factor_config(Nprot = 100,weight_missing = 0.1)
#'
#' Contrasts <- c("c1" = "TreatmentA - TreatmentB",
#'                "C2" = "BackgroundX- BackgroundZ",
#'                "c3" = "`TreatmentA:BackgroundX` - `TreatmentA:BackgroundZ`",
#'                "c4" = "`TreatmentB:BackgroundX` - `TreatmentB:BackgroundZ`"
#'                )
#' mh <- prolfqua::MissingHelpers$new(dd$data, dd$config, prob = 0.8,weighted = TRUE)
#' mh$get_stats()$interaction |> table()
#' mh$get_contrast_estimates(Contrasts)
#'
MissingHelpers <- R6::R6Class(
  "MissingHelpers",

  public = list(
    #' @field data data
    data = NULL,
    #' @field config config
    config = NULL,
    #' @field prob quantile of groups with one observed value to estimate LOD
    prob = 0.5,
    #' @field keep stats which might be time consuming to compute
    stats = NULL,
    #' @field weighted should we weight the LOD
    weighted = TRUE,
    # constructor/initializer
    #' @description
    #' initialize
    #' @param data data
    #' @param config config
    #' @param prob default 0.5, median of groups with one observed value
    #' @param weighted should group average be computed used weighting, default TRUE.
    initialize = function(data, config, prob = 0.5, weighted = TRUE)
    {
      self$data = data
      self$config = config
      self$prob = prob
      self$weighted = weighted
    },
    get_stats = function(){
      if (is.null(self$stats)) {
        self$stats = prolfqua::summarize_stats_factors(self$data, self$config)
      }
      return(self$stats)
    },
    get_LOD = function(){
      LOD <- self$get_stats() |> dplyr::filter(nrMeasured == 1) |>
        dplyr::summarize(LOD = quantile(meanAbundance, probs = self$prob ,na.rm = TRUE)) |>
        pull()
      return(LOD)
    },
    impute_weighted_lod = function(){
      toimp <- self$get_stats()
      toimp$meanAbundanceZero <- ifelse(is.na(toimp$meanAbundance), 0, toimp$meanAbundance)
      impDat <- toimp |> mutate(meanAbundanceImp = (.data$nrMeasured * .data$meanAbundanceZero + .data$nrNAs * self$get_LOD()) / .data$nrReplicates  )
      return(impDat)
    },
    impute_lod = function(){
      toimp <- self$get_stats()
      toimp$meanAbundanceImp <- ifelse(is.na(toimp$meanAbundance), self$get_LOD(), toimp$meanAbundance)
      return(toimp)
    },
    get_poolvar = function(prob = 0.75) {
      if (self$weighted) {
        impDat <- self$impute_weighted_lod()
      } else {
        impDat <- self$impute_lod()
      }
      pooled <- prolfqua::poolvar(impDat, self$config, method = "V1")
      pooled <- dplyr::select(pooled ,-all_of(c(self$config$table$factor_keys_depth()[1],"var")))

      pooled_zero <- pooled[pooled$df > 0,]
      meandf <- pooled_zero |> summarize(
        n = 1, df = 1,
        sd = quantile(sd, prob = prob, na.rm = TRUE),
        sdT = quantile(sdT, prob = prob, na.rm = TRUE))

      minsd <- 1
      meandf$sd <-  ifelse(meandf$sd > 0, meandf$sd, minsd)
      meandf$sdT <-  ifelse(meandf$sdT > 0, meandf$sdT, minsd)

      pooled <- pooled |> mutate(sd = ifelse(is.na(sd) ,meandf$sd, sd))
      pooled <- pooled |> mutate(sdT = ifelse(is.na(sdT) , meandf$sdT, sdT ))
      pooled <- pooled |> mutate(df = ifelse(df == 0, 1, df))
      return(pooled)
    },
    #' @description
    #' get contrast estimates
    #' @param Contrasts named array with contrasts
    get_contrast_estimates = function(
      Contrasts
    ){
      if (self$weighted) {
        lt <- self$impute_weighted_lod()
      } else {
        lt <- self$impute_lod()
      }
      abundance_column = "meanAbundanceImp"
      hierarchy_keys <- self$config$table$hierarchy_keys()
      imp <- lt |> pivot_wider(id_cols = hierarchy_keys,
                               names_from = interaction,
                               values_from = !!sym(abundance_column))


      imputed <- prolfqua::get_contrast(ungroup(imp), hierarchy_keys, Contrasts)
      imputed$avgAbd <- (imputed$group_1 + imputed$group_2)/2
      imputed <- imputed |> dplyr::rename(
        !!paste0(abundance_column,"_group_1") := "group_1",
        !!paste0(abundance_column, "_group_2") := "group_2")

      nr <- lt |> mutate(is_missing = ifelse( .data$nrNAs == .data$nrReplicates , 1 , 0) )
      nr <- nr |> pivot_wider(id_cols = hierarchy_keys, names_from = interaction, values_from = .data$is_missing)
      nrs <- prolfqua::get_contrast(ungroup(nr), hierarchy_keys, Contrasts)

      nrs <- nrs |> select(all_of(c(hierarchy_keys,"contrast", "estimate" )))
      nrs <- nrs |> rename(indic = estimate)
      imputed <- inner_join(imputed, nrs, by = c(hierarchy_keys, "contrast"))


      nrMeasured <- lt |> pivot_wider(id_cols = hierarchy_keys,
                                      names_from = interaction,
                                      values_from = nrMeasured)
      nrMeasured <- prolfqua::get_contrast(ungroup(nrMeasured),  hierarchy_keys, Contrasts)
      nrMeasured <- nrMeasured |> select(all_of(c(hierarchy_keys, "contrast", nrMeasured_group_1 = "group_1", nrMeasured_group_2  = "group_2")))
      imputed <- inner_join(imputed, nrMeasured, by = c(hierarchy_keys, "contrast"))

      imputed2 <- imputed |> mutate(estimate = ifelse(.data$indic < 0 & .data$estimate < 0, 0, .data$estimate))
      imputed2 <- imputed2 |> mutate(estimate = ifelse(.data$indic > 0 & .data$estimate > 0, 0, .data$estimate))

      if (FALSE) {
        imputed2 <- prolfqua::aggregate_contrast(
          ungroup(imputed2),
          subject_Id =  self$config$table$hierarchy_keys_depth())
      }
      return(imputed2)
    },

    get_contrasts = function(Contrasts, confint = 0.95, all = FALSE) {
      imputed <- self$get_contrast_estimates(Contrasts = Contrasts)
      pooled <- self$get_poolvar()
      result <- inner_join(imputed, pooled, by = self$config$table$hierarchy_keys())
      result <- private$add_p_values(result, confint = confint)
      if (!all) {
        result <- select(result, -all_of( c("nrMeasured" , "mean" ,"n.groups", "n", "meanAll") ) )
      }
      return(result)
    }

  ),
  private = list(
    add_p_values = function(result, confint = 0.95, all = TRUE){

      result <- dplyr::mutate(result, statistic = .data$estimate / .data$sdT,
                              p.value = 2*pt(abs(.data$statistic), df = .data$df, lower.tail = FALSE))
      prqt <- -qt((1 - confint)/2, df = result$df)
      result$conf.low <- result$estimate  - prqt * (result$sdT)
      result$conf.high <- result$estimate + prqt * (result$sdT)

      return(result)
    }
  )
)





