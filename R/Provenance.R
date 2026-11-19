# Provenance ----

#' Record of problems and dropped items in one object
#'
#' Each LFQData and model object keeps its own record; nothing is copied between objects. A result that is the
#' same data (a transformation, removed rows, a clone) continues from a copy of the record. `report()` turns the
#' entries recorded since the last report into one summary: a `prolfqua_fit_failed` warning for failures, one
#' `prolfqua_message` for everything else.
#' @noRd
Provenance <- R6::R6Class(
  "Provenance",
  private = list(
    .entries = NULL,
    .reported = logical(0)
  ),
  public = list(
    initialize = function() {
      private$.entries <- tibble::tibble(
        step = character(0),
        event = character(0),
        n = integer(0),
        ids = list(),
        reason = character(0),
        failure = logical(0),
        response = character(0)
      )
    },
    # @param step function or method that recorded the entry
    # @param event what happened, e.g. "rows removed"
    # @param ids affected IDs
    # @param reason why
    # @param failure TRUE for failures (fits, convergence), reported as a warning
    # @param n count; defaults to the number of distinct IDs
    # @param response response column of an LFQData, NA otherwise
    record = function(step, event, ids, reason, failure = FALSE, n = length(unique(ids)), response = NA_character_) {
      if (n == 0L) {
        return(invisible(self))
      }
      private$.entries <- dplyr::bind_rows(
        private$.entries,
        tibble::tibble(
          step = step,
          event = event,
          n = as.integer(n),
          ids = list(unique(as.character(ids))),
          reason = reason,
          failure = failure,
          response = response
        )
      )
      private$.reported <- c(private$.reported, FALSE)
      invisible(self)
    },
    entries = function() {
      private$.entries
    },
    report = function() {
      new <- private$.entries[!private$.reported, ]
      private$.reported[] <- TRUE
      failures <- new[new$failure, ]
      if (nrow(failures) > 0L) {
        rlang::warn(.format_provenance(failures), class = "prolfqua_fit_failed")
      }
      others <- new[!new$failure, ]
      if (nrow(others) > 0L) {
        inform_prolfqua(.format_provenance(others))
      }
      invisible(self)
    }
  )
)

#' One line per entry: "step: n event (e.g. ids): reason"
#' @param entries provenance entries
#' @noRd
.format_provenance <- function(entries) {
  lines <- vapply(
    seq_len(nrow(entries)),
    function(i) {
      ids <- entries$ids[[i]]
      example <- toString(utils::head(ids, 5))
      if (length(ids) > 5L) {
        example <- paste0(example, ", ...")
      }
      paste0(entries$step[i], ": ", entries$n[i], " ", entries$event[i], " (e.g. ", example, "): ", entries$reason[i])
    },
    character(1)
  )
  paste(c(lines, "See $provenance()."), collapse = "\n")
}

#' Distinct IDs of key columns, joined with "~"
#' @param data data frame
#' @param keys key columns
#' @noRd
.key_ids <- function(data, keys) {
  unique(do.call(paste, c(as.data.frame(data)[, keys, drop = FALSE], sep = "~")))
}
