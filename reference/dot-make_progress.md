# Build a progress reporter for the per-protein / per-contrast loops

Internal factory used by \`model_analyse()\` and the contrast loops. The
returned object exposes \`\$tick(len = 1, tokens = list())\`, the same
contract as \`progress::progress_bar\`, so call sites are agnostic to
which reporter they received.

## Usage

``` r
.make_progress(
  total,
  label = NULL,
  reporter = getOption("prolfqua.progress", NULL)
)
```

## Arguments

- total:

  number of iterations; \`\<= 0\` yields a no-op reporter.

- label:

  short label naming the pass (e.g. "firth multi-peptide").

- reporter:

  reporter selector; see Details. Defaults to
  \`getOption("prolfqua.progress", NULL)\`.

## Value

an object with a \`\$tick()\` method.

## Details

The \`reporter\` argument (defaulting to the \`prolfqua.progress\`
option) selects the mode:

- \`NULL\` – the legacy \`progress::progress_bar\` (terminal,
  unchanged).

- a \`function(i, total, label)\` – a throttled user callback, e.g. one
  that calls \`logger::log_info()\`.

- \`"message"\` – a throttled \`message()\` heartbeat with elapsed/ETA.
