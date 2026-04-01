# Limpa Integration Analysis & Plan

> Based on: limpa source code (`/Users/wolski/projects/prolfqua_fml/limpa`), both limpa papers (Li & Smyth 2023, Li, Cobbold & Smyth 2025), prolfqua's limma integration pattern, and the SKILL.md for adding models.

## References

- https://github.com/SmythLab/limpa
- https://bioconductor.org/packages/release/bioc/html/limpa.html
- Source: `/Users/wolski/projects/prolfqua_fml/limpa`
- Papers: `/Users/wolski/projects/prolfqua_fml/addinfo/OxforSmythLimpa_ref1.pdf` (2023 DPC theory), `/Users/wolski/projects/prolfqua_fml/addinfo/SmythLimpaBioarxv.pdf` (2025 package paper)

---

## Answers to Questions from TODO_limpa_integration.md

### Q1: Does limpa do imputation on precursor/peptide level, or only protein level?

**Both, depending on the function used:**

- **`dpcQuant(y, protein.id, dpc)`** — Takes precursor/peptide-level input (matrix with NAs, rows=precursors, cols=samples) and produces **protein-level** output. This is the main pipeline: precursor-to-protein aggregation with probabilistic missing value handling.

- **`dpcQuantByRow(y, dpc)`** — Operates at **whatever level the input is** (precursor, peptide, protein). No aggregation; each row is quantified/imputed independently. Useful for PTM/isoform analyses or when you already have protein-level data but want uncertainty estimates.

So limpa fully supports both peptide-level and protein-level analysis. The DPC itself is fitted at precursor level and is a global parameter for the entire dataset.

### Q2: Is limpa an Aggregator or an Imputer in prolfqua terms?

**It is primarily an Aggregator (precursor → protein) that handles missing values as part of the aggregation.**

`dpcQuant()` combines two responsibilities that prolfqua currently separates:
1. Peptide-to-protein rollup (like `get_Aggregator()`)
2. Missing value handling (like `get_Imputer()`)

However, `dpcQuantByRow()` is purely an imputer (no aggregation).

**Recommendation for prolfqua integration:** Limpa should be integrated as a **modelling method** (facade), not as an Aggregator. The reason: limpa's aggregation is tightly coupled to its probabilistic model, and the output feeds directly into its DE pipeline. Separating the aggregation from the DE would lose the standard error propagation that is limpa's key advantage. The facade should accept precursor-level LFQData and handle the entire pipeline internally (DPC estimation → protein quantification → DE testing).

An alternative simpler integration: use `dpcQuant()` as a standalone aggregator and then feed the result into any prolfqua model. This loses the vooma precision weighting but is still a good aggregation method.

COMMENT: But all of the prolfqua models accept weights. Is this different from the precision weights?

**Response to COMMENT:** Yes, prolfqua models accept weights — but there's a distinction between what prolfqua currently passes and what limpa computes:

- **prolfqua's current weights (`nr_children`)**: Integer peptide counts per protein per sample. Simple proxy: "more peptides = more reliable estimate."

- **limpa's standard errors from `dpcQuant()`**: Continuous, per-protein, per-sample posterior SDs from the Bayesian quantification. These encode *which specific precursors* were missing, the DPC-based probability of those missing values, and the posterior uncertainty. Much more informative than simple counts.

- **limpa's vooma precision weights** (from `voomaLmFitWithImputation()`): Not the raw SEs themselves, but **derived from a bivariate lowess trend** of `sqrt(residual SD) ~ average_intensity + SE_predictor`. This is analogous to how `voom` computes weights from the mean-variance relationship rather than using raw counts directly. The weights end up as `w_ij = 1/f(fitted_value)^4`.

**So there are two feasible integration levels:**

1. **Minimal (Aggregator only):** Run `dpcQuant()` as aggregator, then pass the SE column as weights to any prolfqua model via `strategy_lm(..., weights = "limpa_se")` or `strategy_limma(..., weights = "limpa_se")`. This is easy to implement and already compatible with the existing weight machinery. We'd just need to convert limpa's `EList$other$standard.error` back to a tidy column in the LFQData. The weights would be `1/SE^2` (inverse variance weighting).

2. **Full (Facade):** Run the complete limpa pipeline (`dpc` → `dpcQuant` → `voomaLmFitWithImputation` → `eBayes`). This gets the bivariate vooma precision weights, the DF correction for fully imputed groups, and the full limpa methodology. More code, but statistically better.

**Recommendation:** Start with option 2 (full facade) since it's the principled approach and ~200-300 lines. But also expose `dpcQuant` as an aggregation strategy (option 1) for users who want to mix limpa aggregation with other models. The two are not mutually exclusive.

### Q3: Is imputation a separate step from protein abundance inference in limpa?

**No — they are unified.** Limpa does not impute missing values in the traditional sense. Instead, missing values are **probabilistically integrated out** during protein inference using the Detection Probability Curve (DPC). Each missing observation contributes its log-probability of being missing to the likelihood, rather than being filled with a fabricated number.

The workflow stages are:
1. **DPC estimation** (`dpc()`) — global, dataset-wide, fully independent
2. **Protein quantification** (`dpcQuant()`) — per-protein Bayesian maximum posterior estimation; missing values contribute via DPC, not imputation
3. **Differential expression** (`dpcDE()` / `voomaLmFitWithImputation()`) — limma-based with precision weights from step 2

Steps 1 and 2 are independent from step 3. Step 2 produces a complete matrix (no NAs) with standard errors.

### Q4: Is protein inference separated from the modelling?

**Yes, completely.** `dpcQuant()` produces a standard limma `EList` object with:
- `$E` — protein-level log-intensities (complete, no NAs)
- `$other$standard.error` — per-protein, per-sample standard errors
- `$other$n.observations` — observation counts
- `$genes$NPrec` — number of precursors per protein

This `EList` can be passed to **any** downstream model: limpa's own `dpcDE()`, standard `limma::lmFit()`, or any other method. The protein quantification and DE testing are fully decoupled.

### Q5: Can I use limpa imputation and then fit a different model (limma, prolfqua lm)?

**Yes.** After `dpcQuant()`, the output is a standard `EList`. You can:

```r
# Option 1: Standard limma (ignoring SEs)
fit <- limma::lmFit(y.protein, design)
fit <- limma::eBayes(fit)

# Option 2: limpa's vooma (using SEs as precision weights)
fit <- limpa::voomaLmFitWithImputation(y.protein, design)
fit <- limma::eBayes(fit)

# Option 3: Convert back to long format for prolfqua lm/rlm
# (would need to convert EList back to tidy data.frame)
```

The **key trade-off**: using standard limma or prolfqua's lm after `dpcQuant()` discards the standard error information. Limpa's advantage comes from propagating those SEs as precision weights via `voomaLmFitWithImputation()`. Without that, you still get a good aggregation (better than median polish for missing-heavy data) but lose the variance modelling benefit.

### Q6: Can I disable moderation and apply DEqMS moderation instead?

**Partially.** Limpa itself does not apply empirical Bayes moderation — that's done by limma's `eBayes()` downstream. Limpa only applies **precision weighting** (vooma-style). So:

- You can skip `eBayes()` and apply DEqMS instead
- `dpcQuant()` output includes `$genes$NPrec` (number of precursors per protein) which is exactly what DEqMS needs for count-dependent variance moderation
- In prolfqua terms: wrap `ContrastsLimpa` with `ContrastsModeratedDEqMS` using `NPrec` as the count column

However, the vooma precision weights and DEqMS serve somewhat different purposes:
- **Vooma weights**: Per-sample, per-protein precision from quantification uncertainty
- **DEqMS**: Per-protein variance moderation based on peptide count

They are complementary and could potentially be combined (vooma weights during fitting, DEqMS for variance moderation).

### Q7: Once inferred missingness is handled, is limpa any different from limma?

**The DE testing is essentially limma with two additions:**

1. **Precision weights from quantification uncertainty** — `voomaLmFitWithImputation()` extends limma's vooma by using a bivariate variance predictor: both average intensity AND standard error from `dpcQuant()`. Standard limma only uses average intensity for the mean-variance trend.

2. **DF correction for fully imputed proteins** — When all values in a group are missing (fully imputed), limpa adjusts the residual degrees of freedom using the hat matrix, preventing the model from being overconfident about imputed groups.

The output is a standard `MArrayLM` object. All limma downstream tools work: `eBayes()`, `topTable()`, `contrasts.fit()`, `duplicateCorrelation()`, `makeContrasts()`, etc. The innovation is in the weights and DF correction, not in the testing framework itself.

### Q8: Similarity of internal data structures to limma

**Very high.** Limpa uses limma's own data structures throughout:

| Stage | Input | Output |
|-------|-------|--------|
| DPC estimation | numeric matrix | list (beta0, beta1, hyperparams) |
| Protein quantification | EList/matrix + protein IDs | **EList** (limma class) |
| DE testing | EList + design matrix | **MArrayLM** (limma class) |

This means we can reuse significant code from the limma integration in prolfqua. The `ModelLimma` and `ContrastsLimma` classes already know how to work with `MArrayLM` objects.

---

## Architecture Comparison: limma vs limpa in prolfqua

| Aspect | limma integration | Proposed limpa integration |
|--------|-------------------|---------------------------|
| Input level | Protein (aggregated) | Precursor/peptide (pre-aggregation) |
| Strategy | `StrategyLimma` (formula, trend, robust, weights) | `StrategyLimpa` (formula, DPC params, priors) |
| Builder | `build_model_limma()` | `build_model_limpa()` |
| Model class | `ModelLimma` | `ModelLimpa` (or reuse `ModelLimma`) |
| Contrasts class | `ContrastsLimma` | `ContrastsLimpa` (or reuse `ContrastsLimma`) |
| Facade | `ContrastsLimmaFacade` | `ContrastsLimpaFacade` |
| Wide pivot | `lfqdata$to_wide()` | `lfqdata$to_wide()` (at precursor level) |
| Design matrix | `model.matrix()` from annotation | Same |
| Fit function | `limma::lmFit()` | `limpa::dpcQuant()` → `limpa::voomaLmFitWithImputation()` |
| Moderation | `limma::eBayes()` | `limma::eBayes()` (same) |
| Output | `MArrayLM` | `MArrayLM` (same) |

### Key Difference: Input Level

The limma facade expects **aggregated** (protein-level) LFQData where `subject_Id == hierarchy_keys`. Limpa needs **precursor-level** data. This is the main architectural challenge.

**Options:**
1. **Facade accepts precursor-level LFQData directly** — The facade handles the pivot at precursor level, passes the matrix to `dpcQuant()`, then to `voomaLmFitWithImputation()`. This is the cleanest approach but requires the facade to handle a different input shape than other facades.

2. **Two-step: Aggregator + Facade** — Create a `LimpaAggregator` that wraps `dpcQuant()` and returns aggregated LFQData with SE columns, then pass to a limma-like facade. Cleaner separation but loses the tight coupling between aggregation and DE.

**Recommendation: Option 1** — A single facade that accepts precursor-level input and handles the full pipeline. The facade documentation should clearly state this difference from other facades.

---

## Proposed Integration Design

### Layer 1: Strategy

```r
StrategyLimpa <- R6::R6Class("StrategyLimpa",
  public = list(
    formula = NULL,        # R formula object
    model_name = NULL,     # "limpa"
    trend = NULL,          # passed to eBayes
    robust = NULL,         # passed to eBayes
    # limpa-specific:
    dpc_slope = NULL,      # optional fixed DPC slope (NULL = estimate)
    prior_logFC = NULL,    # optional prior on fold changes
    b1_upper = NULL,       # upper bound for DPC slope
    initialize = function(modelstr, model_name = "limpa",
                          trend = FALSE, robust = FALSE,
                          dpc_slope = NULL, prior_logFC = NULL,
                          b1_upper = 1) { ... }
  )
)

strategy_limpa <- function(modelstr, ...) {
  StrategyLimpa$new(modelstr, ...)
}
```

### Layer 2: Builder + Model

```r
build_model_limpa <- function(lfqdata, strategy) {
  # 1. Pivot precursor-level data to wide matrix
  wide <- lfqdata$to_wide(as.matrix = TRUE)
  expr_matrix <- wide$data           # precursors x samples
  annotation <- wide$annotation      # sample metadata
  rowdata <- wide$rowdata            # precursor + protein IDs

  # 2. Extract protein IDs from rowdata
  protein_id <- rowdata[[protein_id_column]]

  # 3. Estimate DPC
  dpc_est <- limpa::dpc(expr_matrix, b1.upper = strategy$b1_upper)

  # 4. Quantify proteins
  y_protein <- limpa::dpcQuant(expr_matrix, protein.id = protein_id, dpc = dpc_est)

  # 5. Build design matrix
  rhs <- formula(delete.response(terms(strategy$formula)))
  design <- model.matrix(rhs, data = annotation)

  # 6. Fit DE model
  fit <- limpa::voomaLmFitWithImputation(y_protein, design)
  fit <- limma::eBayes(fit, trend = strategy$trend, robust = strategy$robust)

  # 7. Build dummy lm for contrast extraction
  dummy_model <- .build_dummy_lm(design, annotation, strategy$formula)

  # 8. Return ModelLimpa
  ModelLimpa$new(fit, design, strategy$formula, subject_Id,
                 strategy$model_name, rowdata_protein,
                 strategy$trend, strategy$robust, dummy_model)
}
```

**Key insight:** Since `voomaLmFitWithImputation()` returns a standard `MArrayLM`, we may be able to **reuse `ModelLimma`** directly. The only difference is in the builder function. This would significantly reduce code duplication.

### Layer 3: Contrasts

If `ModelLimpa` is essentially `ModelLimma` with a different builder, we can **reuse `ContrastsLimma`** entirely. The contrast computation is identical since both work on `MArrayLM` objects.

### Layer 4: Facade

```r
ContrastsLimpaFacade <- R6::R6Class("ContrastsLimpaFacade",
  inherit = ContrastsInterface,
  public = list(
    model = NULL,
    contrast = NULL,
    initialize = function(lfqdata, modelstr, contrasts,
                          trend = FALSE, robust = FALSE,
                          dpc_slope = NULL, ...) {
      # NOTE: lfqdata must be at precursor/peptide level (not aggregated)
      response <- lfqdata$config$get_response()
      full_formula <- paste(response, modelstr)
      strat <- strategy_limpa(full_formula, trend = trend, robust = robust,
                              dpc_slope = dpc_slope, ...)
      self$model <- build_model_limpa(lfqdata, strat)
      self$contrast <- ContrastsLimma$new(self$model, contrasts)
    },
    get_contrasts = function() { self$contrast$get_contrasts() },
    get_missing = function() { ... },
    get_Plotter = function() { ContrastsPlotter$new(self) },
    to_wide = function() { ... }
  )
)
```

### Registration in build_contrast_analysis()

```r
# In build_contrast_analysis.R, add to method choices:
method = c("lm", "rlm", "lmer", "lm_missing", "lm_impute",
           "limma", "limma_impute", "deqms", "ropeca", "firth",
           "limpa")   # NEW

# In switch:
limpa = ContrastsLimpaFacade$new(lfqdata, modelstr, contrasts, ...),
```

---

## Implementation Considerations

### 1. Input Level Mismatch

All current facades expect aggregated (protein-level) LFQData. The limpa facade needs precursor-level data. This creates a usability concern: users must know NOT to aggregate before calling the limpa method.

**Mitigation options:**
- Validate in the facade: check that `lfqdata$subject_Id()` is NOT equal to the full hierarchy (i.e., data is nested)
- Document clearly
- Consider a helper that checks input level

### 2. Hierarchy Mapping

Limpa needs a `protein.id` vector mapping each precursor row to its protein. In prolfqua, this is encoded in the `AnalysisConfiguration`:
- `config$table$hierarchy_keys()` gives all levels (protein, peptide, precursor)
- The protein-level key is typically `config$table$hierarchy_keys()[1]`
- After pivoting to wide, the protein ID column is in `rowdata`

### 3. Standard Error Propagation

Limpa's key advantage is propagating SEs from `dpcQuant()` into DE. If we split this across prolfqua steps (aggregation, then modelling), we need to ensure SEs flow through. The cleanest approach: keep everything in a single facade.

### 4. DESCRIPTION Dependencies

```
Suggests:
    limpa (>= 1.0.0)
```

Conditional availability check in the facade:
```r
if (!requireNamespace("limpa", quietly = TRUE)) {
  stop("Package 'limpa' is required for method='limpa'.")
}
```

Same pattern as limma dependency (check how limma is declared in DESCRIPTION).

### 5. Reuse Potential

Since `voomaLmFitWithImputation()` returns a standard `MArrayLM`:
- **Reuse `ModelLimma`** — same class, different builder
- **Reuse `ContrastsLimma`** — identical contrast computation
- Only need new: `StrategyLimpa`, `build_model_limpa()`, `ContrastsLimpaFacade`

This is ~200-300 lines of new code, not a major addition.

---

## Alternative Integration: Limpa as Aggregator Only

For users who want limpa's protein quantification but their own DE method:

```r
# Use limpa for aggregation, then any prolfqua method
lfqdata_protein <- limpa_aggregate(lfqdata_precursor)  # wraps dpcQuant()
result <- build_contrast_analysis(lfqdata_protein, modelstr, contrasts, method = "limma")
```

This would require:
- A new aggregation strategy class (like `AggregateMedpolish`, `AggregateTopN`)
- Converting limpa's `EList` output back to tidy LFQData format
- Mapping `NPrec` and `standard.error` to appropriate columns

**Downside:** Loses the vooma precision weighting. The SEs from `dpcQuant()` would be discarded unless we add a mechanism to carry them through as weights.

---

## Summary of Answers

| Question | Answer |
|----------|--------|
| Precursor or protein level? | **Both**: `dpcQuant()` for precursor→protein, `dpcQuantByRow()` for same-level |
| Aggregator or Imputer? | **Primarily Aggregator** with built-in probabilistic missing handling |
| Imputation separate from inference? | **No** — unified via DPC likelihood (not traditional imputation) |
| Protein inference separate from modelling? | **Yes** — `dpcQuant()` output is a standard EList |
| Can use different model after limpa imputation? | **Yes** — any limma-compatible method works on the EList |
| Can disable moderation and use DEqMS? | **Yes** — skip eBayes, use DEqMS with NPrec counts |
| Is limpa different from limma after imputation? | **Mostly limma** + vooma precision weights from SEs + DF correction |
| Data structures similar to limma? | **Identical** — uses limma's EList and MArrayLM throughout |

## Files to Create/Modify

| File | Action | Content |
|------|--------|---------|
| `R/ContrastsLimpa.R` | **Create** | `StrategyLimpa`, `strategy_limpa()`, `build_model_limpa()` |
| `R/ContrastsFacades.R` | **Modify** | Add `ContrastsLimpaFacade` |
| `R/build_contrast_analysis.R` | **Modify** | Add `"limpa"` to method choices and switch |
| `DESCRIPTION` | **Modify** | Add `limpa` to Suggests |
| `tests/testthat/test-ContrastsLimpa.R` | **Create** | Tests for limpa integration |

## Open Design Questions

1. **Should we also expose limpa as an Aggregator?** (separate from the facade) This allows `dpcQuant()` followed by any downstream method but loses SE propagation.

2. **Should `ModelLimpa` be a new class or reuse `ModelLimma`?** Since both wrap `MArrayLM`, reuse is likely cleaner. The only question is whether we need limpa-specific methods on the model.

3. **How to handle the input-level difference?** The facade needs precursor-level data while others need protein-level. Should we validate, warn, or document?

4. **Should we support `dpcQuantByRow()` separately?** This enables precursor-level imputation without aggregation — useful if users want to run precursor-level DE (unusual but possible).

5. **DPC parameter tuning:** Should users be able to supply a pre-estimated DPC, or should the facade always estimate it fresh? Pre-estimated DPC allows sharing across analyses.
