# ContrastsMissing----
#' Compute contrasts with group mean imputation (DEPRECATED)
#'
#' `r lifecycle::badge("deprecated")` (placeholder — prolfqua does not
#' yet depend on the lifecycle package; consider this a soft
#' deprecation notice. See \emph{Deprecated} below.)
#'
#' If there are no observations in one of the groups for some of the
#' proteins, the group mean cannot be estimated. Therefore, assuming
#' that the observation is missing because the protein abundance is
#' below the detection limit, we substitute the unobserved group with
#' the median of protein abundances observed only in one sample of the
#' group. The variance of a protein is estimated using the pooled
#' variance of all observations of all groups.
#'
#' @section Deprecated:
#' `ContrastsMissing` is a pre-fitting group-mean substitution: it does
#' not fit a per-protein model, just stamps the group-mean delta with a
#' pooled-variance t-test. It is superseded by
#' \code{\link{build_model_impute}} which refits failed/singular
#' proteins with LOD-imputed values and borrowed variance per protein.
#' New code should prefer the \code{lm_impute} facade
#' (\code{\link{ContrastsLMImputeFacade}}) or any of its
#' \code{limma_impute} / \code{limma_voom_impute} cousins. Constructing
#' \code{ContrastsMissing} emits a \code{.Deprecated} warning at
#' \code{initialize} time.
#'
#' Still reachable via \code{model = "lm_missing"} for users who want
#' to explicitly run the group-mean fallback; the construction emits a
#' \code{.Deprecated} warning each time. Will be removed once a
#' merge-style replacement (e.g. an \code{lm_impute_missing} facade
#' composing \code{build_model} with \code{build_model_impute}) lands.
#'
#' @seealso \code{\link{build_model_impute}},
#'   \code{\link{ContrastsLMImputeFacade}}
#'
#' @family modelling
#' @return An R6 class generator.
#' @keywords internal
#' @export
#' @examples
#' Nprot <- 120
#' istar <- prolfqua::sim_lfq_data_protein_config(Nprot = Nprot,weight_missing = .4)
#' istar$data$abundance |> is.na() |> sum()
#' protIntensity <- istar$data
#' config <- istar$config
#'
#'
#' lProt <- LFQData$new(protIntensity, config)
#' lProt$rename_response("transformedIntensity")
#'
#' Contr <- c("dil.b_vs_a" = "group_A - group_Ctrl")
#' csi <- ContrastsMissing$new(lProt, contrasts = Contr)
#' csi$get_contrast_sides()
#'
#' res <- csi$get_contrasts()
#'
#' stopifnot(nrow(res) ==  (protIntensity$protein_Id |> unique() |> length()))
#' res$contrast |> table()
#' stopifnot((res$p.value |> is.na() |> sum()) == 0)
#' plot(res$diff, -log10(res$p.value), pch = ".")
#' csi$column_description()
#' x<- csi$get_Plotter()
#' p <- x$volcano()
#' pdf(file = NULL)
#' print(p)
#' dev.off()
#'
#' dd <- prolfqua::sim_lfq_data_2factor_config(Nprot = 100,weight_missing = 0.1)
#'
#' Contrasts <- c("c1" = "TreatmentA - TreatmentB",
#'                "C2" = "BackgroundX- BackgroundZ",
#'                "c3" = "`TreatmentA:BackgroundX` - `TreatmentA:BackgroundZ`",
#'                "c4" = "`TreatmentB:BackgroundX` - `TreatmentB:BackgroundZ`"
#'                )
#' lProt <- LFQData$new(dd$data, dd$config)
#' lProt$rename_response("transformedIntensity")
#'
#' csi <- ContrastsMissing$new(lProt, contrasts = Contrasts)
#' res <- csi$get_contrasts()
#' pl <- csi$get_Plotter()
#' pdf(file = NULL)
#' pl$volcano()
#' dev.off()
ContrastsMissing <- R6::R6Class(
  "ContrastsMissing",
  inherit = ContrastsInterface,
  private = list(
    method = "V1"
  ),
  public = list(
    #' @field subject_id subject_id e.g. protein_ID column
    subject_id = character(),
    #' @field contrasts array with contrasts (see example)
    contrasts = character(),
    #' @field model_name model name
    model_name = character(),
    #' @field contrast_result data frame with results of contrast computation
    contrast_result = NULL,
    #' @field lfqdata data frame
    lfqdata = NULL,
    #' @field confint confidence interval
    confint = 0.95,
    #' @field p.adjust function to adjust p-values
    p.adjust = NULL,
    #' @field global Take global or local values for imputation
    global = logical(),
    #' @field present default 1, presence in interaction to infer limit of detection.
    present = 1,
    #' @field minsd default 1, if standard deviation can not be estimated, what is the prior minimum sd, default = 1s
    minsd = 1,
    #' @description
    #' initialize
    #' @param lfqdata LFQData
    #' @param contrasts array of contrasts (see example)
    #' @param confint confidence interval
    #' @param p.adjust method for p-value adjustment - default Benjamini Hochberg
    #' @param model_name default "groupAverage"
    initialize = function(
      lfqdata,
      contrasts,
      confint = 0.95,
      p.adjust = prolfqua::adjust_p_values,
      model_name = "groupAverage"
    ) {
      .Deprecated(
        new = "build_model_impute() + ContrastsLMImputeFacade",
        msg = paste(
          "ContrastsMissing is deprecated: it substitutes group means",
          "rather than fitting a model. Prefer build_model_impute",
          "(LOD-imputed per-protein refit with borrowed variance) via",
          "the lm_impute / limma_impute facades. See ?ContrastsMissing",
          "for details."
        )
      )
      self$subject_id <- lfqdata$relevant_hierarchy_keys()
      self$contrasts <- contrasts
      self$model_name <- model_name
      self$lfqdata <- lfqdata
      self$confint <- confint
      self$p.adjust <- p.adjust
    },
    #' @description
    #' get contrasts sides
    #'
    get_contrast_sides = function() {
      parse_contrast_sides(self$contrasts)
    },
    #' @description
    #' table with results of contrast computation
    #' @param all FALSE, do not show all columns (default)
    get_contrasts = function(all = FALSE) {
      if (is.null(self$contrast_result)) {
        if (self$lfqdata$get_config()$hierarchy_depth < length(self$lfqdata$get_config()$hierarchy_keys())) {
          abort_bad_argument(
            "lfqdata",
            "be aggregated first (hierarchy_depth must equal the number of hierarchy keys)",
            not = paste0(
              "hierarchy_depth = ",
              self$lfqdata$get_config()$hierarchy_depth,
              "; hierarchy_keys = ",
              length(self$lfqdata$get_config()$hierarchy_keys())
            )
          )
        } else {
          mh1 <- prolfqua::MissingHelpers$new(
            self$lfqdata$data_long(),
            self$lfqdata$get_config(),
            prob = 0.5,
            weighted = TRUE
          )
          result <- mh1$get_contrasts(Contrasts = self$contrasts, confint = self$confint, all = all)
          result <- self$p.adjust(result, column = "p.value", group_by_col = "contrast", newname = "FDR")
        }
        result <- result |> rename(diff = estimate, sigma = sd, std.error = sdT)
        # group-mean substitution, not a model fit: every row is a fallback
        result$estimate_type <- "missing_fallback"
        result <- mutate(result, modelName = self$model_name, .before = 1)
        result <- dplyr::relocate(result, "estimate_type", .after = "modelName")
        self$contrast_result <- ungroup(result)
      }
      res <- self$contrast_result
      stopifnot(all(super$column_description()$column_name %in% colnames(res)))
      invisible(res)
    },
    #' @description
    #' get ContrastsPlotter
    #' @return Contrast_Plotter
    get_Plotter = function() {
      res <- ContrastsPlotter$new(
        self$get_contrasts(),
        subject_id = self$subject_id,
        volcano = list(list(score = "p.value", thresh = 0.1), list(score = "FDR", thresh = 0.1)),
        histogram = list(list(score = "p.value", xlim = c(0, 1, 0.05)), list(score = "FDR", xlim = c(0, 1, 0.05))),
        modelName = "modelName",
        diff = "diff",
        contrast = "contrast"
      )
      return(res)
    },
    #' @description
    #' convert contrast results to wide format
    #' @param columns value column default p.value
    #' @return data.frame
    to_wide = function(columns = c("p.value", "FDR", "statistic")) {
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_to_wide(
        contrast_minimal,
        subject_id = self$subject_id,
        columns = c("diff", columns),
        contrast = "contrast"
      )
      return(contrasts_wide)
    }
  )
)
