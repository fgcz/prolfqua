# Two facade bugs surfaced by the WU345302 harness

Both surfaced after the `needs` interface was collapsed to `"same"` / `"nested"`
and the prolfquapp peptide-to-protein dispatch was unified into a single
`run_dea`. Reproduce via:

```bash
cd integration_test && make wu345302-facades
cat test-outputs/wu345302_facades/status.tsv
```

15 of 17 facades pass; the two below fail.

## 1. `ContrastsLimpaFacade` over-enforces `config$opt_se`

### Symptom

```
ERROR ContrastsLimpaFacade requires LFQData with config$opt_se set.
      Use AggregateLimpa to produce the input.
```

Log: `integration_test/test-outputs/wu345302_facades/logs/limpa.stdout.log`.

### What the caller did

`prolfquapp::run_dea(software = "prolfquapp.DIANN", model = "limpa")` ran the
normal `"same"` path:

1. `ProteinDataPrep$aggregate()` with the configured aggregator
   (`medpolish` in this fixture).
2. `ProteinDataPrep$transform_data()` (`none`).
3. `ContrastsLimpaFacade$new(lfq, modelstr, contrasts)` via
   `build_contrast_analysis(..., method = "limpa")`.

Step 3 raises because the constructor checks `config$opt_se` and the
input wasn't produced by `AggregateLimpa`.

### Intended semantics (from W.W.)

`limpa` is a **`"same"` facade**. It should accept whatever LFQData prolfqua's
normal aggregation pipeline produced and operate at that hierarchy level —
protein-aggregated input → protein-level FC, or precursor-level input →
precursor-level FC (precursor imputation lives in the same code path; the
hierarchy of the input is what determines the output level).

`limpa_nested` is the facade that owns the `AggregateLimpa` step: it takes
nested peptide/precursor input, runs `AggregateLimpa` internally, and emits
protein-level contrasts.

### Fix

In `R/ContrastsFacades.R::ContrastsLimpaFacade`:

- Drop the `config$opt_se` requirement from the constructor.
- Make the facade work directly on the supplied (already-aggregated) LFQData,
  the same way every other `"same"` facade does.
- If the limpa fit needs a `SummarizedExperiment`-shaped input, build it from
  the supplied LFQData on the fly inside the facade — do not require the
  caller to pre-set `config$opt_se`.

In `R/ContrastsChildToParentFacades.R::ContrastsLimpaNestedFacade`:

- This is where the `AggregateLimpa` pre-step belongs. Confirm it already runs
  `AggregateLimpa` on the nested input and produces an `opt_se`-equipped
  LFQData for the limpa fit. If not, move that logic here.

### Verification

```bash
cd integration_test && bash scripts/run_wu345302_facades.sh
# limpa row should switch from fail → ok in status.tsv
```

`limpa_nested` already passes today; do not regress it.

## 2. `ContrastsLmerNestedFacade` errors with "subscript out of bounds"

### Symptom

```
INFO  model formula: normalized_abundance ~ G_
determine linear functions:
ERROR subscript out of bounds
ERROR Stack trace:
ERROR No traceback available
```

Log: `integration_test/test-outputs/wu345302_facades/logs/lmer_nested.stdout.log`.

### What the caller did

`prolfquapp::run_dea(software = "prolfquapp.DIANN_PEPTIDE", model = "lmer_nested")`.
Peptide-level LFQData reached the facade through the normal nested path
(`hierarchy_depth = 1`, lfq_data_peptide branch). All other nested facades
(`firth_nested`, `limpa_nested`, `ropeca_nested`) succeed on the same input.

### Likely root cause

The error is raised inside the lmer-nested contrast/linfct computation.
Generic message, no traceback in the harness log. Most likely a vector
indexing assumption that doesn't hold for the WU345302 fixture (e.g.
a contrast row indexing into a coefficient vector that's shorter than
expected when some random-effect levels collapse).

### Investigation

1. Re-run in an interactive R session with the same fixture:
   ```r
   options(error = function() traceback(3))
   pkgload::load_all("prolfquapp")
   prolfquapp::run_dea(
     indir   = "integration_test/fixtures/diann_wu345302/out-DIANN",
     dataset = "integration_test/fixtures/diann_wu345302/dataset_saint.csv",
     software = "prolfquapp.DIANN_PEPTIDE",
     config   = prolfquapp::get_config(
       "integration_test/test-outputs/wu345302_facades/configs/config_lmer_nested.yaml"
     )
   )
   ```
2. Capture the traceback. Most likely landing inside
   `contrasts_linfct()` or `compute_contrast()` operating on the lmer
   coefficient matrix.
3. Fix the indexing or pad the missing levels.

### Verification

`lmer_nested` row in the harness status flips from fail → ok.

## Out of scope

- The `needs` interface itself stays at two values. These bugs are
  per-facade implementation issues, not interface issues.
- prolfquapp-side `run_dea` dispatch is correct (proven by the 15 passes
  including `limpa_nested`, `firth_nested`, `ropeca_nested`). No change there.
