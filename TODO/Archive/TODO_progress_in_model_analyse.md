# Progress reporting in `model_analyse` (and the firth/contrast loops)

## TL;DR

A progress bar **already exists** in `model_analyse()` and is correctly ticked
by every strategy's `model_fun`. The problem is **it is silently suppressed in
non-interactive / docker runs**, and it writes to `stderr` rather than to the
`logger` log file that the app and users actually watch. So during a long
`firth_nested` fit (see WU347684 / WU347680: 5–7 h, 34 cores pegged, log frozen
at the "model formula" line) there is **zero signal** that anything is
progressing. The fix is not "add a progress bar" — it is "make the existing one
visible in batch runs and/or emit coarse progress to the log."

## Evidence / current state

- `R/tidyMS_build_model.R:637` — `model_analyse()` creates
  `pb <- progress::progress_bar$new(total = nrow(nested_proteins))` and passes
  `pb` into `purrr::map(data, model_strategy$model_fun, pb = pb)` (`:639`).
  The map is **sequential** over proteins, so tick-based progress is reliable.
- Every strategy ticks it:
  - `StrategyLM$model_fun` — `R/tidyMS_R6_Modelling.R:85-86`
  - `StrategyRLM$model_fun` — `R/tidyMS_R6_Modelling.R:189-190`
  - `StrategyLmer$model_fun` — `R/tidyMS_R6_Modelling.R:294-295`
  - `StrategyLogistf$model_fun` — `R/logistf.R:361-362`
- A second bar exists in the firth contrast step:
  `contrasts_linfct_firth()` — `R/logistf.R:35,48`.

### Why it is invisible in production
1. `progress::progress_bar$new()` is called **without `force = TRUE`** and with
   the default `stream = stderr()`. The `progress` package only renders when the
   stream is a terminal (`isatty`). In the docker `Rscript --vanilla` run
   (stderr piped to logs, not a tty) the bar is **silently disabled**.
2. Even if rendered, it goes to `stderr` via carriage-return repaint. The
   `prolfqua_*.log` the user tails is the **`logger`** stream
   (`prolfquapp/R/R6_DEAnalyse.R:240` `logger::log_info("model formula: ...")`).
   prolfqua itself does **not** depend on `logger`, so the bar and the log are
   two disconnected channels.
3. For `firth_nested`, `build_model_logistf()` calls `model_analyse()` **twice**
   (`R/logistf.R:199` models2 = multi-peptide proteins; `:213` models1 =
   single-peptide), so there are two independent, unlabelled bars plus the
   contrast bar — and the expensive one (models2, formula augmented with the
   peptide key `~ G_ + Subject_ + <peptide_Id>`) gives no hint of its cost.

## Goal

When a fit runs for minutes/hours in a batch/docker context, the log should show
that it is alive and roughly how far along it is (count, %, ETA) — without
requiring an interactive terminal, and without prolfqua taking a hard `logger`
dependency.

## Proposed approach

- **A. Injectable progress callback.** Give `model_analyse(..., progress = NULL)`
  an optional callback invoked every N proteins (or every X seconds). Default
  keeps current bar behaviour; prolfquapp passes a callback that does
  `logger::log_info()`. Keeps prolfqua `logger`-free and makes progress land in
  the real log.
- **B. Coarse milestone logging at the prolfquapp call site.** In
  `prolfquapp` (`R6_DEAnalyse.R`, around the facade build at `:240`) log
  `start fitting N proteins`, then completion with elapsed time. Guarantees a
  heartbeat in the watched log even if a single model fit takes a long time
  between ticks. Combine with A for intra-fit granularity.
- **C. Optional forced terminal progress.** For interactive/diagnostic use,
  support `force = TRUE` (and a sane `format`/`show_after`) on the legacy
  `progress::progress_bar` path. This is useful at a terminal, but it is not the
  production fix: carriage-return repaint on `stderr` can be noisy and still does
  not land in the structured `logger` stream users watch.

### Recommendation
Implement A+B first: a prolfqua callback/reporter API plus prolfquapp
`logger::log_info()` heartbeat messages. Keep the current terminal progress bar
as the default path, and add forced terminal progress only as an opt-in
diagnostic mode.

## Touch points

- `R/tidyMS_build_model.R:637` — `model_analyse()` bar/reporter (label with
  `model_name`, add `progress` callback param; optionally support forced
  terminal progress on the legacy bar path).
- `R/logistf.R:35` — `contrasts_linfct_firth()` bar.
- `R/ContrastFirth.R:97,118` — `ContrastsFirth$get_linfct()` has separate
  model1/model2 linfct loops used by the nested Firth facade; instrument these
  too so the post-fit contrast setup is not silent.
- `R/logistf.R:199,213` — two `model_analyse()` calls; pass distinct labels
  (e.g. "firth multi-peptide" vs "firth single-peptide") so the two bars are
  distinguishable.
- `prolfquapp/R/R6_DEAnalyse.R:~240` — add start/elapsed `logger::log_info()`
  around the facade build (milestone heartbeat), and/or supply the progress
  callback from option A.

## Acceptance criteria

- A `firth_nested` (or any) run in docker writes periodic progress to the
  watched log (count/%, and ideally ETA), proving liveness.
- The two firth `model_analyse` passes are individually identifiable in output.
- No new hard dependency forced on prolfqua (callback/log progress keeps prolfqua
  `logger`-free; the logging happens in prolfquapp).
- Interactive behaviour unchanged (bar still renders in a terminal).

## Detailed implementation plan — callback/log progress

### Why callback/log progress
Confirmed from the WU347684 workunit log: prolfqua `message()` output **does**
reach the captured docker/slurm log (`completing cases`, `starting aggregation`,
`Joining ...` all appear), but the `progress::progress_bar` does **not** (non-tty
suppression). So a callback that emits through `message()` / `logger` lands
exactly in the log the operator is tailing. This keeps prolfqua free of a hard
`logger` dependency while making progress visible in batch.

### Design: a duck-typed "reporter" with `$tick()`
The strategies already call `pb$tick()` (`tidyMS_R6_Modelling.R:86,190,295`,
`logistf.R:362`). Keep that contract; just pass an object that *also* exposes
`$tick(len = 1, tokens = list())` (same signature as `progress::progress_bar`),
so call sites can share the same interface.

Important semantic detail: the current strategy methods tick **before** the
expensive fit starts. For batch logs, progress should represent completion (or
explicitly say "starting"). Prefer moving ticks to after each successful/failed
fit attempt, or emit paired `starting protein i` / `finished protein i` events in
the reporter. Otherwise a run can report 100% and then spend a long time inside
the final model fit.

Add an internal constructor:

```r
# default reporter = current progress_bar; override via option or arg
.make_progress <- function(total, label = NULL,
                           reporter = getOption("prolfqua.progress", NULL)) {
  if (is.null(reporter)) {                       # legacy path, unchanged
    return(progress::progress_bar$new(total = total))
  }
  if (is.function(reporter)) {                   # user callback(i, total, label)
    return(.callback_reporter(total, label, reporter))  # throttled by wall time
  }
  if (identical(reporter, "message")) {          # simple batch-visible path
    return(.message_reporter(total, label))      # message() every X s / P %
  }
  stop("unknown prolfqua.progress reporter")
}
```

- `.callback_reporter` / `.message_reporter` track an internal counter + start
  time (`Sys.time()`), and only emit when **>= N seconds elapsed since last emit**
  (e.g. 10 s) or **every P%** (e.g. 2%). This avoids 21k chatty lines.
- Emit content: `"<label>: protein i/total (P%), elapsed Xm, ETA ~Ym"` where
  `ETA = elapsed/i * (total - i)`.
- `total == 0` → no-op reporter.

### Wiring
1. `model_analyse(..., progress = getOption("prolfqua.progress", NULL))`
   (`tidyMS_build_model.R:625-639`): build reporter via `.make_progress(total,
   label = model_name, reporter = progress)` and pass into the existing
   `purrr::map(..., pb = pb)`. Keep the NULL/default path behaviour unchanged.
2. Thread an optional `label` so callers can name the pass.
3. Ensure ticks are completion-oriented. For `model_fun()` implementations, wrap
   the fit in a local result and call `pb$tick()` after the `tryCatch()` returns,
   or add separate reporter methods for started/finished if start messages are
   wanted.
4. `build_model_logistf` (`logistf.R:199,213`): pass distinct labels
   `"firth multi-peptide"` / `"firth single-peptide"` so the two `model_analyse`
   passes are distinguishable in the log.
5. `contrasts_linfct_firth` (`logistf.R:35`): replace its standalone
   `progress_bar$new` with `.make_progress(..., label = "firth contrasts")`.
6. `ContrastsFirth$get_linfct()` (`ContrastFirth.R:97,118`): replace the model1
   and model2 standalone bars with labelled reporters, e.g. `"firth linfct
   single-peptide"` and `"firth linfct multi-peptide"`.

### prolfquapp side (keeps prolfqua logger-free)
- In `prolfquapp` startup / `CMD_DEA_V2.R` (or `R6_DEAnalyse$build_*`), set once:
  ```r
  options(prolfqua.progress = function(i, total, label) {
    label <- if (is.null(label) || identical(label, "")) "fit" else label
    logger::log_info("{label}: {i}/{total} ({round(100*i/total)}%) ...")
  })
  ```
  prolfqua only ever calls a user-supplied function — no `logger` import.
- Add a start/elapsed heartbeat around the facade build in
  `R6_DEAnalyse.R:~240` (already logs the model formula there) so even a
  zero-tick stall is bracketed by `start fitting N proteins` / `done in Xs`.

### Tests (prolfqua `tests/testthat/`)
- Callback reporter: collect `(i,total)` into an env; assert ticks are monotonic
  and final `i == total` after a small `build_contrast_analysis(..., "firth")`.
- `"message"` reporter emits in a non-interactive session: capture via
  `capture.output(type = "message", ...)` / `withCallingHandlers`.
- Default (`NULL`) path: assert results identical to pre-change (regression).

### Risk / scope
- `model_analyse` is core to **all** Wald facades — keep the `NULL` default path
  untouched and gate every new behaviour behind a non-NULL reporter/option.
- Reporter must be concurrency-safe **if** the protein loop is ever parallelised
  (it is currently a sequential `purrr::map`, so fine today — note it).
- Reporter `$tick()` must accept `len`/`tokens` args for drop-in compatibility
  with `progress_bar`.

### Phasing
1. `.make_progress` + reporters + unit tests (no call-site change).
2. Wire into `model_analyse` (default unchanged).
3. Apply to the firth loops + labels.
4. prolfquapp option + heartbeat.
5. Docs: document `prolfqua.progress` option and the `callback(i, total, label)`
   contract.

## Notes / related

- The 34-core / 129-thread usage comes from multithreaded BLAS inside each fit,
  **not** from parallelising the protein loop — the loop is a sequential
  `purrr::map`, which is exactly why per-protein ticks are trustworthy.
- Related diagnosis: `prolfquapp/TODO` discussion of `firth_nested` being
  pathologically slow on large peptide-level DIA-NN input with a paired design.
  Progress reporting does not fix the cost, but it turns a "frozen, is it
  hung?" experience into an informed "it's protein 4,000/21,000, ETA Xh"
  decision.
