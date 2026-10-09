#' Group statistics and limit of detection
#'
#' Computes per-group statistics of the abundances and, from them, the limit of
#' detection (LOD) used by \code{\link{build_model_impute}} to impute groups
#' without observations.
#'
#' @return An R6 class generator.
#' @export
#' @examples
#' dd <- prolfqua::sim_lfq_data_protein_config(Nprot = 100,weight_missing = 2)
#' mh <- prolfqua::MissingHelpers$new(dd$data, dd$config, prob = 0.8)
#' xx <- mh$get_stats()
#' xx <- mh$get_lod()
#'
#' dd <- prolfqua::sim_lfq_data_2factor_config(Nprot = 100,weight_missing = 0.1)
#' mh <- prolfqua::MissingHelpers$new(dd$data, dd$config, prob = 0.8)
#' mh$get_stats()$interaction |> table()
#' mh$get_lod()
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
    #' @field stats data.frame with group statistics
    stats = NULL,
    #' @description
    #' initialize
    #' @param data data
    #' @param config config
    #' @param prob default 0.5, median of groups with one observed value
    initialize = function(data, config, prob = 0.5) {
      self$data <- data
      self$config <- config
      self$prob <- prob
    },
    #' @description
    #' get data.frame with statistics
    #' @return data.frame
    get_stats = function() {
      if (is.null(self$stats)) {
        self$stats <- prolfqua::summarize_stats_factors(
          LFQData$new(self$data, self$config)
        )
      }
      return(self$stats)
    },
    #' @description
    #' determine limit of detection
    #' computes quantile of abundances in groups with a single observation;
    #' without such groups, of the partly observed groups with the fewest
    #' observations, and when every group is complete or empty, the lowest
    #' group mean
    #' @return numeric LOD
    get_lod = function() {
      if (is.null(private$.lod_cache)) {
        private$.lod_cache <- .lod_from_group_stats(self$get_stats(), self$prob)
      }
      return(private$.lod_cache)
    }
  ),
  private = list(
    .lod_cache = NULL
  )
)

# Groups seen once are the nearest the data comes to its detection limit. A
# nearly complete dataset can have none while one feature still misses a whole
# group and needs a value to refit, so the partly observed groups with the
# fewest observations stand in, and when every group is complete or empty, the
# lowest group mean.
.lod_from_group_stats <- function(stats, prob) {
  measured <- dplyr::filter(stats, .data$nrMeasured > 0)
  partial <- dplyr::filter(measured, .data$nrMeasured == 1 | .data$nrMeasured < .data$nrReplicates)
  if (nrow(partial) == 0) {
    return(if (nrow(measured) == 0) NA_real_ else min(measured$meanAbundance))
  }
  fewest <- dplyr::filter(partial, .data$nrMeasured == min(.data$nrMeasured))
  unname(stats::quantile(fewest$meanAbundance, probs = prob, na.rm = TRUE))
}
