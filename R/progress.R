# Injectable progress reporting for the per-protein / per-contrast loops.
#
# `model_analyse()` and the contrast loops historically created a
# `progress::progress_bar` directly. That bar renders only to an interactive
# terminal (it is silently disabled on a non-tty stderr, e.g. a docker
# `Rscript` run) and it writes via carriage-return repaint rather than to a
# structured log. `.make_progress()` keeps that terminal bar as the default but
# lets a batch consumer (e.g. prolfquapp) opt in to a reporter that emits
# through `message()` or a user callback, so progress lands in the watched log.
#
# Every reporter exposes the same duck-typed contract the strategy
# `model_fun()`s and the contrast loops already call: `$tick(len = 1,
# tokens = list())`. This keeps prolfqua free of a hard `logger` dependency --
# it only ever invokes a user-supplied `function(i, total, label)`.

# A reporter that does nothing -- used for empty loops (total <= 0).
.null_reporter <- function() {
  list(tick = function(len = 1, tokens = list()) invisible(NULL))
}

# Format a number of seconds as a short human-readable duration ("45s", "2.3m",
# "1.4h").
.format_duration <- function(secs) {
  if (!is.finite(secs) || secs < 0) {
    return("?")
  }
  if (secs < 60) {
    return(sprintf("%ds", round(secs)))
  }
  if (secs < 3600) {
    return(sprintf("%.1fm", secs / 60))
  }
  sprintf("%.1fh", secs / 3600)
}

# A reporter that throttles output: it emits only on the final tick, or once at
# least `min_interval` seconds have elapsed since the last emit, or once at
# least `min_pct` percent of additional progress has accrued. This avoids tens
# of thousands of chatty lines while still proving liveness.
#
# `callback` mode calls `callback(i, total, label)` -- the documented public
# contract for the `prolfqua.progress` option. `message_mode` emits a
# `message()` heartbeat with elapsed time and ETA.
.throttled_reporter <- function(
  total,
  label = NULL,
  callback = NULL,
  message_mode = FALSE,
  min_interval = 10,
  min_pct = 2
) {
  label <- if (is.null(label) || identical(label, "")) "fit" else label
  state <- new.env(parent = emptyenv())
  state$i <- 0L
  state$start <- Sys.time()
  state$last_time <- state$start
  state$last_pct <- 0

  emit <- function(now) {
    if (!is.null(callback)) {
      callback(state$i, total, label)
      return(invisible(NULL))
    }
    elapsed <- as.numeric(difftime(now, state$start, units = "secs"))
    eta <- if (state$i > 0) elapsed / state$i * (total - state$i) else NA_real_
    message(sprintf(
      "%s: %d/%d (%d%%), elapsed %s, ETA ~%s",
      label,
      state$i,
      total,
      round(100 * state$i / total),
      .format_duration(elapsed),
      .format_duration(eta)
    ))
    invisible(NULL)
  }

  list(
    tick = function(len = 1, tokens = list()) {
      state$i <- state$i + len
      now <- Sys.time()
      pct <- 100 * state$i / total
      due_time <- as.numeric(difftime(now, state$last_time, units = "secs")) >= min_interval
      due_pct <- (pct - state$last_pct) >= min_pct
      if (state$i >= total || due_time || due_pct) {
        state$last_time <- now
        state$last_pct <- pct
        emit(now)
      }
      invisible(NULL)
    }
  )
}

#' Build a progress reporter for the per-protein / per-contrast loops
#'
#' Internal factory used by `model_analyse()` and the contrast loops. The
#' returned object exposes `$tick(len = 1, tokens = list())`, the same contract
#' as `progress::progress_bar`, so call sites are agnostic to which reporter
#' they received.
#'
#' The `reporter` argument (defaulting to the `prolfqua.progress` option)
#' selects the mode:
#' \itemize{
#'   \item `NULL` -- the legacy `progress::progress_bar` (terminal, unchanged).
#'   \item a `function(i, total, label)` -- a throttled user callback, e.g. one
#'     that calls `logger::log_info()`.
#'   \item `"message"` -- a throttled `message()` heartbeat with elapsed/ETA.
#' }
#'
#' @param total number of iterations; `<= 0` yields a no-op reporter.
#' @param label short label naming the pass (e.g. "firth multi-peptide").
#' @param reporter reporter selector; see Details. Defaults to
#'   `getOption("prolfqua.progress", NULL)`.
#' @return an object with a `$tick()` method.
#' @keywords internal
.make_progress <- function(total, label = NULL, reporter = getOption("prolfqua.progress", NULL)) {
  if (length(total) != 1 || is.na(total) || total <= 0) {
    return(.null_reporter())
  }
  if (is.null(reporter)) {
    return(progress::progress_bar$new(total = total))
  }
  if (is.function(reporter)) {
    return(.throttled_reporter(total, label, callback = reporter))
  }
  if (identical(reporter, "message")) {
    return(.throttled_reporter(total, label, message_mode = TRUE))
  }
  stop("unknown prolfqua.progress reporter: use NULL, a function(i, total, label), or 'message'")
}
