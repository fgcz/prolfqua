# ContrastsFacadeBase ----

#' Base class for contrast analysis facades
#'
#' The extension point for new modelling backends. Every method of \code{\link{build_contrast_analysis}} is a
#' facade class that inherits from \code{ContrastsFacadeBase}; a package adds a method by writing such a class and
#' registering it with \code{\link{register_facade}} (usually from its \code{.onLoad()}). The base class provides
#' the uniform API: \code{get_contrasts()} (stamped with the method name), \code{get_missing()},
#' \code{get_Plotter()}, \code{to_wide()} and \code{get_config()}.
#'
#' A subclass implements \code{initialize(lfqdata, modelstr, contrasts, ...)}:
#' \itemize{
#'   \item call the private \code{.setup(lfqdata, contrasts, key, modelstr)} first; it checks the input shape,
#'     the contrasts and the design, stores a copy of the input, and returns the full model formula string
#'   \item set \code{self$model} to the fitted model and \code{self$contrast} to the contrast object the facade
#'     delegates to (e.g. \code{\link{ContrastsModerated}})
#'   \item set the private field \code{.needs} to \code{"nested"} if the facade takes child-level input
#'     (e.g. peptides) and reports parent-level contrasts; the default is \code{"aggregated"}
#' }
#' The limma family additionally sets \code{self$.drop_na_diff <- TRUE} to drop rows whose fold change could not
#' be estimated.
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @seealso \code{\link{register_facade}}, \code{\link{list_facades}}
#' @examples
#' # a facade: unweighted linear model per protein, moderated contrasts
#' ContrastsUnweightedLMFacade <- R6::R6Class(
#'   "ContrastsUnweightedLMFacade",
#'   inherit = ContrastsFacadeBase,
#'   public = list(
#'     initialize = function(lfqdata, modelstr, contrasts) {
#'       full_formula <- private$.setup(lfqdata, contrasts, "lm_unweighted", modelstr)
#'       self$model <- build_model(lfqdata, strategy_lm(full_formula))
#'       self$contrast <- ContrastsModerated$new(Contrasts$new(self$model, contrasts))
#'     }
#'   )
#' )
#' istar <- sim_lfq_data_protein_config(n_proteins = 10)
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' facade <- ContrastsUnweightedLMFacade$new(lfqdata, "~ group_", c(AvsCtrl = "group_A - group_Ctrl"))
#' head(facade$get_contrasts())
#' # registered from a package's .onLoad(), build_contrast_analysis() finds the method by name:
#' # register_facade("lm_unweighted", class = "ContrastsUnweightedLMFacade", needs = "aggregated",
#' #                 package = "mypackage")
#' @srrstats {RE1.3} Facade results retain subject identifiers, and get_missing reports absent subject/contrast pairs.
#' @srrstats {RE4.8} Facades copy LFQData; matrix fits retain E, weights, observed masks and sample annotation.
ContrastsFacadeBase <- R6::R6Class(
  "ContrastsFacadeBase",
  inherit = ContrastsInterface,
  public = list(
    #' @field model fitted model object
    model = NULL,
    #' @field contrast inner contrast object the facade delegates to
    contrast = NULL,
    #' @field .lfqdata snapshot of the input LFQData (used by get_missing)
    .lfqdata = NULL,
    #' @field .contrast_names names of the requested contrasts
    .contrast_names = NULL,
    #' @field facade_name registry key stamped into the modelName column
    facade_name = NULL,
    #' @field .drop_na_diff drop rows with NA diff (limma family); default FALSE
    .drop_na_diff = FALSE,
    #' @description get the \code{\link{ContrastConfiguration}}; defaults to
    #'   the inner contrast object's configuration
    get_config = function() {
      if (is.null(self$config)) {
        self$config <- self$contrast$get_config()
      }
      self$config
    },
    #' @description get contrast results, stamped with the facade key
    #' @param ... passed to the inner contrast object's get_contrasts
    get_contrasts = function(...) {
      res <- .stamp_facade_identity(self$contrast$get_contrasts(...), self$facade_name)
      if (isTRUE(self$.drop_na_diff) && "diff" %in% colnames(res)) {
        res <- res[!is.na(res$diff), , drop = FALSE]
      }
      res
    },
    #' @description get subject x contrast pairs that could not be estimated
    get_missing = function() {
      .compute_missing(
        self$.lfqdata,
        self$.contrast_names,
        self$get_contrasts(),
        contrast_col = self$get_config()$contrast_col
      )
    },
    #' @description get \code{\link{ContrastsPlotter}}
    #' @param ... passed to the inner contrast object's get_Plotter
    get_Plotter = function(...) self$contrast$get_Plotter(...),
    #' @description convert results to wide format
    #' @param ... passed to the inner contrast object's to_wide
    to_wide = function(...) self$contrast$to_wide(...),
    #' @description problems and dropped items of this analysis, read from the input data, the aggregated data
    #'   (\code{hurdle_nested}, \code{limpa_nested}) and the model(s), plus the subjects absent from the contrasts
    #'   with the reason
    #' @return tibble with columns source, step, event, n, ids (list), reason, failure, response
    provenance = function() {
      absent <- .absent_provenance(self$get_missing(), self$model, self$.lfqdata$subject_id())
      dplyr::bind_rows(self$.recorded_provenance(), dplyr::mutate(absent, source = "contrasts", .before = 1))
    },
    #' @description internal: the entries recorded in the input data, the aggregated data and the model(s),
    #'   without evaluating contrasts (used by \code{print()})
    #' @return tibble as \code{provenance()}, without the absent subjects
    .recorded_provenance = function() {
      sources <- list(`input data` = private$.input$provenance())
      if (!identical(self$.lfqdata, private$.input)) {
        sources[["aggregated data"]] <- self$.lfqdata$provenance()
      }
      models <- if (is.null(self$model)) {
        list(`model (intensity)` = self$contrast$intensity$model, `model (detection)` = self$contrast$detection$model)
      } else {
        list(model = self$model)
      }
      for (name in names(models)) {
        sources[[name]] <- models[[name]]$provenance()
      }
      dplyr::bind_rows(sources, .id = "source")
    }
  ),
  private = list(
    # the input LFQData copy made by .setup(); hurdle_nested and limpa_nested replace .lfqdata afterwards
    .input = NULL,
    # input shape this class takes: "aggregated" (subject_id equals the hierarchy) or "nested" (child-level input);
    # nested facades override it. The registry's `needs` mirrors it for callers that look up a method by name.
    .needs = "aggregated",
    # Shared initialize() prologue: check the input shape required by `.needs`, store the input and the contrast
    # names, record the facade key, and return the full model formula string.
    .setup = function(lfqdata, contrasts, key, modelstr) {
      if (!inherits(lfqdata, "LFQData")) {
        abort_bad_argument("lfqdata", "be an LFQData object")
      }
      .assert_model_input(lfqdata)
      .assert_named_contrasts(contrasts)
      formula <- .as_model_formula(modelstr)
      if (!is.character(modelstr) || length(formula) != 2) {
        abort_bad_argument("modelstr", "be one right-hand-side formula string, e.g. `~ group_`")
      }
      if (identical(private$.needs, "nested")) {
        .assert_nested_facade_input(lfqdata, class(self)[1])
      } else {
        .assert_aggregated_facade_input(lfqdata, class(self)[1])
      }
      .assert_contrast_groups(lfqdata, contrasts)
      .assert_full_rank_design(lfqdata, modelstr)
      self$.lfqdata <- lfqdata$get_copy()
      private$.input <- self$.lfqdata
      self$.contrast_names <- names(contrasts)
      self$facade_name <- key
      paste(lfqdata$response(), modelstr)
    }
  )
)


#' Provenance entries for the subjects absent from the contrasts
#'
#' The reason is the failure text when the subject's fit failed, otherwise "contrast not estimable".
#' @param missing get_missing() result: subject columns and contrast
#' @param model fitted model (ModelInterface) or NULL
#' @param subject_id subject key columns
#' @noRd
.absent_provenance <- function(missing, model, subject_id) {
  log <- Provenance$new()
  if (nrow(missing) == 0L) {
    return(log$entries())
  }
  missing$id <- do.call(paste, c(as.data.frame(missing)[, subject_id, drop = FALSE], sep = "~"))
  missing$reason <- "contrast not estimable"
  if (!is.null(model)) {
    entries <- model$provenance()
    failed <- entries[entries$event == "fits failed", ]
    for (i in seq_len(nrow(failed))) {
      missing$reason[missing$id %in% failed$ids[[i]]] <- paste("fit failed:", failed$reason[i])
    }
  }
  contrast_col <- setdiff(colnames(missing), c(subject_id, "id", "reason"))[1]
  groups <- unique(missing[, c(contrast_col, "reason")])
  for (i in seq_len(nrow(groups))) {
    rows <- missing[[contrast_col]] == groups[[contrast_col]][i] & missing$reason == groups$reason[i]
    log$record(
      "get_missing",
      paste0("subjects absent from contrast ", groups[[contrast_col]][i]),
      missing$id[rows],
      groups$reason[i]
    )
  }
  log$entries()
}
