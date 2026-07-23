# Review: `plan_binomial.md` — native `binomial_nested` model

Reviewer notes on `TODO/plan_binomial.md`. Findings were checked against the current
`prolfqua` source (version 1.6.3) and cross-verified with an adversarial claim-by-claim pass.

## Overall assessment

The *structural* half of the plan is sound and well-reasoned — reuse-first data prep, the
choice of `quasibinomial()`, and avoiding bespoke `Model`/`Contrasts` classes mostly hold up
against the code. But there are **two concrete code-level defects**, a **misleading set of
facade templates**, and a cluster of **statistical decisions the plan presents as settled that
need to be made explicitly**. Nothing is fatal to the approach, but the plan as written would
lead an implementer into a wrong turn on the facade and a silently-wrong moderated statistic.

## ✅ Related pre-existing bug — found during this review and FIXED

Investigating defect #1 (below) surfaced that the **same bug class already shipped in the
`rlm` backend**, independent of the binomial plan. It has been **fixed** as part of this work:

- **Symptom:** `StrategyRLM$sigma()` returned `stats::sigma(model)` (the ordinary-residual
  scale) while `MASS::rlm` builds its `vcov()` / `summary()` standard errors from the robust
  scale `model$s`. Because `ContrastsModerated` rescales the Wald statistic by
  `sigma / sqrt(var.post)`, the mismatch left a spurious per-protein factor in every
  **moderated** `rlm` statistic, p-value, FDR and confidence bound.
- **Scope:** only the `rlm` facade's *moderated* output was affected. `lm` (deviance = RSS, so
  `stats::sigma` is correct) and `rfit` (custom coherent `sigma.rfit_prolfqua = tauhat`) were
  already fine; `rlm`'s *unmoderated* `diff` / `std.error` were always correct.
- **Fix (root cause, one line):** `StrategyRLM$sigma <- function(model) model$s` in
  `R/tidyMS_R6_Modelling.R`, with a reproducing coherence test in
  `tests/testthat/test-Model.R` (asserts `strategy_rlm(...)$sigma(fit) == fit$s` and
  `sigma * unscaled == sqrt(diag(vcov))`) and a `NEWS.md` bullet under 1.6.3. All model,
  contrast, facade, DEqMS and vectorize test suites pass.

**Implication for the binomial plan:** the fix establishes the `sigma ↔ vcov` coherence
invariant (and its regression test). For a quasibinomial glm that invariant is *already*
satisfied by `stats::sigma` (R returns the Pearson dispersion for glm objects), so
`StrategyBinomial` needs no change here — see defect #1, which is **retracted**.

## Post-implementation recheck (verified against the merged code)

The `binomial_nested` implementation (`R/BinomialNested.R`, facade in
`R/ContrastsChildToParentFacades.R`, moderation floor in `R/tidyMS_moderation.R`) was
re-verified against this review:

- **Defect #1 — RETRACTED.** `sigma = stats::sigma` is correct for the quasibinomial glm
  (empirically `stats::sigma == sqrt(summary$dispersion)`, coherent with `vcov()`).
- **Defect #5 — addressed.** The facade sets a family-aware `ContrastConfiguration`
  (`supports_dea_qc = FALSE`); roxygen documents `diff` = log odds ratio and `avgAbd` = average
  linear predictor on the log-odds scale, so log-odds output is not mislabeled as abundance.
- **Variance floor — implemented as analyzed, including the ordering subtlety.** `binomial_bound`
  maps to `variance_floor = 1`; floored rows get `df.prior <- Inf` *after* the pre-existing
  "all `df.prior` infinite" fallback, so it is not clobbered.
- **Tests green.** `test-BinomialNested.R` passes 65/65.

The remaining "statistical decisions" below were deliberate author choices (documented in
`NEWS.md`): the symmetric pseudo-count for separated fits and the posterior-dispersion bound.

## What the plan gets right

- **Reusing `.prepare_logistf_lfqdata()`** (plan §"Reuse-first design"). The helper at
  `R/logistf.R:87-94` is already backend-neutral (`get_copy` → `complete_cases` →
  `encode_bin_resp` → `set_config_value`). `complete_cases()` runs *before* encoding, so the
  undetected (`0`) cells exist. The only in-repo coupling is one
  `prolfqua:::.prepare_logistf_lfqdata` reference in `tests/testthat/test-ContrastsFacades.R:286`
  — trivial to update on rename.
- **The count summary as the only new data operation.** No existing function produces a
  per-parent×sample detected/undetected count — `HierarchyCountsSample` is sample-level,
  `summarize_stats` is condition-level. `available = n()` is a legitimate binomial denominator
  because the completed grid gives each parent×sample exactly `nr_children` rows.
- **Choosing `quasibinomial()` rather than `binomial()`** — the cleverest part of the plan.
  Because dispersion is estimated, `summary.glm` emits a **t-table (`Pr(>|t|)` → `Pr...t..`)**
  and exposes `df.residual`/`sigma`, exactly the lm-shaped surface that
  `Model`/`Contrasts`/`AnovaExtractor` (`get_anova_df(test="F")` → `Pr..F.`) expect. A
  `Pr...z..` mismatch worry applies only to plain `binomial()` and is avoided by this choice.
- **The contrast core is family-agnostic.** `linfct_from_model` / `.model_coeff_matrix`
  (`R/tidyMS_contrasts.R:23-49`) read `coef()`/`model.matrix()` by name, and `.compute_contrast`
  reads only `coefficients/vcov/df.residual/sigma` — all provided by a glm. `diff` is a correct
  log-odds ratio; `avgAbd` a correct mean-logit.
- **The moderation-floor mechanism is computationally safe.** `pt(x, df=Inf) = pnorm`,
  `qt(p, df=Inf) = qnorm`, and FDR reads only `p.value`, so `df → Inf` cannot break the
  downstream.

## Code-level defects (fix before implementing)

### 1. ~~`sigma = stats::sigma` is the wrong scale for a quasibinomial glm~~ — RETRACTED (my claim was WRONG)

**This claim was incorrect and is withdrawn.** It was rejected during implementation, and an
empirical check confirms the rejection was right.

The moderation rescale at `R/tidyMS_moderation.R`,
`moderated.statistic = statistic * sigma / sqrt(var.post)`, is only coherent when
`std.error = sigma × (unscaled)`. I asserted that `stats::sigma(glm)` returns the
deviance-based scale `sqrt(deviance/df.residual)` and would therefore break this. **It does
not.** For a glm, R's `sigma()` returns the *Pearson dispersion* — the same quantity `vcov()`
uses. Verified on a quasibinomial fit:

```
stats::sigma(fit)          = 1.149071
sqrt(summary$dispersion)   = 1.149071   <- equal: Pearson, what vcov() uses
sqrt(deviance/df.residual) = 1.177689   <- NOT what stats::sigma returns
stats::sigma * unscaled    == sqrt(diag(vcov(fit)))   <- reconstructs the SE exactly
```

So `sigma = function(model) stats::sigma(model)` (as implemented in `StrategyBinomial`) is
**correct and coherent** for the quasibinomial glm; no change is needed.

My error was over-generalizing a *class-specific* behavior: the identical reasoning is genuinely
correct for `rlm` (where `stats::sigma` returns the ordinary-residual scale, ≠ the robust scale
`$s` that `vcov()` uses — the real bug we fixed) but wrong for `glm` (where `stats::sigma` already
equals the Pearson dispersion). The `sigma ↔ vcov` coherence invariant still matters; `glm` simply
satisfies it via `stats::sigma` out of the box. The `rlm` fix and its regression test stand.

### 2. The cited facade templates are the wrong ones (major)

The plan says the facade "matches how firth_nested / ropeca_nested facades are built." It does
not:

- `firth_nested` (`R/ContrastsChildToParentFacades.R:267-276`) uses a **bespoke `ContrastsFirth`
  adapter** — no generic `Contrasts$new`, no `ContrastsModerated`. Following "reuse firth" as
  "copy firth's facade" leads straight to the bespoke `ContrastsBinomialNested` adapter the plan
  explicitly wants to avoid.
- `ropeca_nested` (`:138-151`) wraps in `ContrastsROPECA` (not `ContrastsModerated`), overrides
  `get_contrasts/get_Plotter/to_wide`, and overrides `config$subject_id`.

The actual structural analog is **`ContrastsLmerNestedFacade`** (`:82-94`): `assert_nested_input
→ build_model → ContrastsModerated$new(Contrasts$new(...)) → set self$config <-
self$contrast$get_config()`, registered `needs="nested"` at `R/ContrastsFacades.R:1058`. The plan
should name **lmer_nested** as the template.

### 3. Be explicit about where count aggregation happens (new wiring, undersold)

Unlike every nested facade (which feeds the *original* nested LFQData to `build_model`), the
binomial facade must first collapse `bin_resp` across children into a parent×sample
detected/undetected table, repackage it as an LFQData, and fit with `subject_Id = parent` so
`build_model` nests correctly and `df.residual = n_samples − n_coef`. Do this in the facade
**before** `build_model` (feeding raw child rows would force `model_fun` to aggregate
internally). The plan's steps imply this but never state the repackaging or the `subject_Id`
choice — and no cited analog demonstrates it.

### 4. The moderation floor lives in `moderated_p_limma`, not `ContrastsModerated`

The `squeezeVar` call is at `R/tidyMS_moderation.R:22`, one level below the R6 class. The floor
threads through three sites (`moderated_p_limma`, `moderated_p_limma_long`, and a new field on
`ContrastsModerated$initialize` + its hardcoded `get_contrasts` call at
`R/ContrastsModerated.R:96-100`). Two caveats:

- `df → Inf` is **not** automatic — `:35` computes `df.total = df + df.prior`, and the floor
  branch must override it *after* bypassing the pre-existing "all `df.prior` infinite" fallback
  at `:25-27` that would otherwise replace `Inf` with a finite value.
- The exported `df` column will then contain `Inf`, so any report that prints/uses it must
  tolerate that.
- The floor must be **supplied by the caller** — the decorator can't infer family.

### 5. Set a family-aware `ContrastConfiguration` in the facade (minor–major)

The generic `Contrasts` never sets `self$config`, so it inherits the default that labels
`avgAbd` as *abundance* and `diff` as a symmetric effect. Here `avgAbd` is a **mean logit
(routinely negative)** and `diff` a **log-odds ratio** — feeding the default config into
prolfquapp's DEA/volcano/abundance-filter path mislabels axes and breaks any `avgAbd > 0`
assumption. The adding-models skill is explicit that the facade should set `self$config` with
correct roles; the plan's "expose the wrapped contrast configuration" would expose the
mislabeled default.

### 6. Minor: `report_columns` and `is_mixed`

Every existing strategy declares them and a test asserts them, though `build_model` /
`model_analyse` don't functionally read them. Set them anyway for consistency (cheap, avoids the
test surprise) — the plan's strategy field list omits both.

## Statistical decisions to make explicit

The plan presents these as closed; each deserves an explicit, defended choice.

- **`prior_count = 0.1` is not a substitute for Firth.** A flat symmetric pseudocount bounds the
  coefficient but doesn't *solve* separation the way the Jeffreys penalty in `logistf` does — and
  it is added regardless of the binomial denominator, so effect-size shrinkage becomes
  `nr_children`-dependent (low-peptide proteins attenuated more than high-peptide ones). Under
  genuine separation the reported log-odds ratio is partly an artifact of `prior_count`.
  Aggregation to counts mitigates this when `nr_children` is large, but single-child proteins
  reduce to Bernoulli where it bites hardest. State whether a denominator-scaled prior or a
  `bayesglm`/Firth-style penalty was considered and rejected.
- **"nested" is arguably a misnomer.** Collapsing to per-sample counts treats a protein's
  precursors as independent binomial trials, discarding child identity — *less* nested than
  `firth_nested`, which carries `peptide_Id` as a covariate. It ignores within-protein detection
  correlation and will tend to overstate precision.
- **quasibinomial dispersion at small residual df.** With ~2–4 samples the dispersion is very
  noisy, and for single-child proteins it is non-identifiable from binary data. This is the input
  to the moderation.
- **`squeezeVar` on a binomial dispersion.** It assumes Gaussian residual variances with
  scaled-χ² sampling; a quasibinomial `σ²` does not follow that model, so the pooled prior and
  moderated p-values are heuristic. The floor-at-1 is a defensible guard against sub-binomial
  overconfidence, but it is grafted onto a Gaussian shrinkage engine. (This compounds defect #1:
  even after fixing the Pearson-scale bug, the statistical justification for shrinking this
  quantity deserves a sentence.)

## Suggested test additions

The plan's test list is strong (it already covers rank deficiency, insufficient residual df, and
"absence of msqrob2"). Add:

- **moderated statistic parity** — assert `moderated.statistic` matches a hand-computed rescale
  using `sqrt(summary(fit)$dispersion)`; this is the regression test for defect #1.
- **separation behavior** — a protein detected in all of one group / none of the other; assert
  the estimate is finite and pin how `prior_count` bounds it.
- **single-child (`nr_children == 1`) protein** — confirm the Bernoulli-per-sample case produces
  sane df/dispersion or is handled.
- **config semantics** — assert the facade's `ContrastConfiguration` labels `diff`/`avgAbd` as
  log-odds, not abundance.

## Bottom line

Proceed with the structural design (reuse prep, `quasibinomial`, generic contrast path,
lmer_nested-style facade), but:

1. fix `sigma` to the Pearson dispersion (`sqrt(summary(fit)$dispersion)`);
2. rewrite the facade section to name `lmer_nested` and describe the count→LFQData repackaging +
   `subject_Id = parent`;
3. set a family-aware `ContrastConfiguration`;
4. convert the four statistical assumptions from assertions into explicit, defended choices —
   especially the `prior_count`-vs-Firth separation question, which is the substantive risk to
   the method's validity.
