# ContrastsMissing / SimpleImpute — Algorithm Documentation & Redesign Plan

## Current Algorithm (MissingHelpers + ContrastsMissing)

### Overview

`ContrastsMissing` (in `R/ContrastsSimpleImpute.R`) is a group-average-based contrast method designed for data with many missing values. Unlike `Contrasts` (which fits per-protein linear models) or `ContrastsLimma` (which fits a limma matrix model), `ContrastsMissing` computes imputed group means directly, forms contrasts by subtraction, and derives test statistics from a pooled variance estimate.

The heavy lifting is done by the `MissingHelpers` R6 class (`R/tidyMS_missigness_V2.R`).

### Files

- `R/ContrastsSimpleImpute.R` — `ContrastsMissing` R6 class (thin wrapper)
- `R/tidyMS_missigness_V2.R` — `MissingHelpers` R6 class (engine)
- `R/tidyMS_stats.R` — `summarize_stats()`, `summarize_stats_factors()`, `poolvar()`, `compute_pooled()`, `pooled_V1()`
- `R/tidyMS_missigness.R` — `get_contrast()`, `complete_cases()`

---

### Step 1: Compute Group Statistics

**Function:** `summarize_stats(pdata, config)` (tidyMS_stats.R:168)

First calls `complete_cases(pdata, config)` to ensure all sample×feature combinations exist (adding NA rows for missing measurements).

Then groups by `hierarchy_keys × factor_keys` and computes per group:

| Column | Formula |
|--------|---------|
| `nrReplicates` | `n()` — total samples in group |
| `nrMeasured` | `sum(!is.na(intensity))` |
| `nrNAs` | `sum(is.na(intensity))` |
| `sd` | `sd(intensity, na.rm=TRUE)` |
| `var` | `var(intensity, na.rm=TRUE)` |
| `meanAbundance` | `mean(intensity, na.rm=TRUE)` — NaN if all NA |
| `medianAbundance` | `median(intensity, na.rm=TRUE)` |
| `interaction` | paste of factor columns (group label) |

For multi-factor designs, `summarize_stats_factors()` calls `summarize_stats()` once for the full interaction and once per individual factor, then `bind_rows()` all results. This allows contrasts at both the interaction level and marginal factor levels.

### Step 2: Estimate LOD (Limit of Detection)

**Method:** `MissingHelpers$get_LOD()` (tidyMS_missigness_V2.R:68)

```r
LOD <- stats |>
  filter(nrMeasured == 1) |>
  summarize(LOD = quantile(meanAbundance, probs = prob))
```

- Selects groups with **exactly one** observed value
- Takes the `prob` quantile (default 0.5 = median) of their mean abundances
- **Rationale:** A single observation in a group is likely just above the detection limit. The median of these values estimates the LOD.
- Returns a single scalar (global LOD, not per-protein or per-group)

### Step 3: Impute Group Means

Two strategies, controlled by `weighted` parameter (default TRUE):

#### 3a. Weighted LOD (default) — `impute_weighted_lod()` (line 78)

```
meanAbundanceZero = ifelse(is.na(meanAbundance), 0, meanAbundance)
meanAbundanceImp = (nrMeasured × meanAbundanceZero + nrNAs × LOD) / nrReplicates
```

This is a weighted average: observed measurements contribute their mean, missing measurements contribute LOD. When a group has all observations, `meanAbundanceImp = meanAbundance`. When all missing, `meanAbundanceImp = LOD`.

#### 3b. Simple LOD — `impute_lod()` (line 91)

```
meanAbundanceImp = ifelse(is.na(meanAbundance), LOD, meanAbundance)
```

Only replaces fully-missing groups with LOD; partially-observed groups keep their observed mean unchanged (NAs ignored in original mean).

### Step 4: Pooled Variance Estimation

**Method:** `MissingHelpers$get_poolvar()` → `poolvar()` → `compute_pooled()` → `pooled_V1()`

#### Per-protein pooled variance (pooled_V1, tidyMS_stats.R:46):

For each protein, given k groups with nrMeasured > 1:

```
SS = Σ (n_i - 1) × var_i     # sum of squares
pool.var = SS / (Σn_i - k)    # pooled variance
df = Σn_i - k                 # degrees of freedom
sd = sqrt(pool.var)
```

**Standard error for contrasts (two-sample):**

```
sdT = sqrt(pool.var × 2 / (pool.n / n.groups))
```

Where `pool.n = Σn_i` and `n.groups` = number of groups with >1 observation. This is an equal-n t-test SE formula using the average group size.

#### Fallback for proteins with insufficient data:

```
1. Filter to proteins with df > 0 AND sd > 0
2. Compute global fallback:
     median_sd  = quantile(all valid sds,  prob=0.75)
     median_sdT = quantile(all valid sdTs, prob=0.75)
3. Replace NA/zero sd  → median_sd
4. Replace NA/zero sdT → median_sdT
5. Replace df=0 → df=1
```

The 75th percentile (not median!) is used as fallback — a conservative choice that overestimates variance for proteins with no variance estimate.

### Step 5: Contrast Estimation

**Method:** `MissingHelpers$get_contrast_estimates()` (tidyMS_missigness_V2.R:129)

1. **Pivot** imputed group means to wide format: rows = proteins, columns = group labels, values = `meanAbundanceImp`

2. **Evaluate contrasts** using `get_contrast()` (tidyMS_missigness.R:112):
   - Parses contrast expressions like `"group_A - group_Ctrl"` using `rlang::parse_expr()`
   - Extracts: `estimate = group_1 - group_2`, `group_1` value, `group_2` value

3. **Compute `avgAbd`** = (group_1 + group_2) / 2

4. **Compute missingness indicator `indic`:**
   - Pivot a binary matrix: `is_missing = 1 if nrNAs == nrReplicates, else 0`
   - Apply the same contrast expression
   - `indic = is_missing_group1 - is_missing_group2`
     - `indic = 0`: both groups observed (or both missing)
     - `indic = -1`: group_1 entirely missing, group_2 observed
     - `indic = +1`: group_1 observed, group_2 entirely missing

5. **Extract nrMeasured per contrast side** for annotation

6. **Zero out invalid contrasts:**
   ```
   if indic < 0 AND estimate < 0 → estimate = 0
   if indic > 0 AND estimate > 0 → estimate = 0
   ```
   **Rationale:** When one group is entirely missing and imputed with LOD, the direction of the fold change is ambiguous. If the observed group's mean is _below_ the LOD-imputed group, the contrast is zeroed out as unreliable. Only reports differences where the observed group is _above_ the imputed group (consistent with the "missing because below LOD" assumption).

### Step 6: Statistical Testing

**Method:** `add_p_values()` (private, tidyMS_missigness_V2.R:200)

```
statistic = estimate / sdT
p.value   = 2 × pt(|statistic|, df=df, lower.tail=FALSE)    # two-tailed
conf.low  = estimate - qt(α/2, df) × sdT
conf.high = estimate + qt(α/2, df) × sdT
```

Where α = 1 - confint (default confint = 0.95).

### Step 7: P-value Adjustment

In `ContrastsMissing$get_contrasts()`:

```r
result <- self$p.adjust(result, column = "p.value", group_by_col = "contrast", newname = "FDR")
```

Default: Benjamini-Hochberg correction, applied per contrast.

### Step 8: Column Renaming

```r
result <- rename(result, diff = estimate, sigma = sd, std.error = sdT)
result <- mutate(result, modelName = "groupAverage", .before = 1)
```

---

### Output Columns

| Column | Description |
|--------|-------------|
| `modelName` | `"groupAverage"` |
| `contrast` | Contrast name |
| `protein_Id` (hierarchy keys) | Protein identifier |
| `diff` | Log-space fold change (group_1 - group_2) |
| `avgAbd` | (group_1 + group_2) / 2 |
| `meanAbundanceImp_group_1` | Imputed mean for group 1 |
| `meanAbundanceImp_group_2` | Imputed mean for group 2 |
| `indic` | Missingness indicator (-1, 0, +1) |
| `nrMeasured_group_1` | Observations in group 1 |
| `nrMeasured_group_2` | Observations in group 2 |
| `df` | Degrees of freedom |
| `sigma` | Pooled standard deviation |
| `std.error` | Standard error of the contrast (sdT) |
| `statistic` | t-statistic |
| `p.value` | Two-tailed p-value |
| `FDR` | BH-adjusted p-value |
| `conf.low`, `conf.high` | Confidence interval bounds |

---

### Assumptions and Limitations

1. **MNAR (Missing Not At Random):** Missing values are assumed to be below the detection limit. This is a strong assumption — if data is MAR, LOD imputation introduces bias.

2. **Global LOD:** One LOD value for all proteins. In practice, different proteins may have different detection thresholds.

3. **Homogeneous variance:** Pooled variance assumes equal variance across groups within a protein.

4. **Equal group sizes for sdT:** The `sdT = sqrt(pool.var × 2 / (pool.n / n.groups))` formula uses average group size, which is approximate for unbalanced designs.

5. **75th percentile fallback:** Proteins without variance estimates receive the 75th percentile of all estimated SDs — conservative but ad hoc.

6. **Zero-out logic:** Contrasts involving one fully-missing group are zeroed when the direction is "wrong" (inconsistent with missing=low assumption). This is conservative but loses power.

7. **No moderation / shrinkage:** Unlike limma or ContrastsModerated, there is no empirical Bayes moderation of variances. Each protein's pooled variance is used as-is (or replaced with the global fallback).

8. **No repeated measures:** Cannot handle blocking or paired designs.

9. **No robust fitting:** Sensitive to outliers in group means.

---

## Proposed Redesign: lm-based or limma-based ContrastsMissing

### Motivation

The current implementation manually computes group means, pooled variance, and t-statistics. This reimplements standard statistical machinery. An lm- or limma-based approach would:

1. **Reuse existing model fitting** — `build_model()` or `build_model_limma()` already handles the linear model
2. **Gain empirical Bayes moderation** (limma path) — variance shrinkage across proteins
3. **Support complex designs** — multi-factor, interactions, blocking
4. **Standardize the output** — same contrast interface as other backends

### Option A: lm-based

- Impute individual missing values (not just group means) with LOD before fitting
- Feed imputed data to `strategy_lm()` → `build_model()` → `Contrasts` / `ContrastsModerated`
- Pro: reuses existing per-protein lm infrastructure, moderated t-statistics available
- Con: imputing individual values at LOD adds artificial precision (variance underestimated for imputed observations)
- Mitigation: could use weighted lm with weight=0 or very low weight for imputed observations

### Option B: limma-based

- Impute individual missing values with LOD
- Feed imputed data to `strategy_limma()` → `build_model_limma()` → `ContrastsLimma`
- Pro: empirical Bayes variance shrinkage, handles complex designs, weights support already implemented
- Con: same individual-value imputation concern; limma's variance shrinkage may partially compensate

### Option C: Hybrid — keep group-mean approach but add moderation

- Keep the current group-mean imputation strategy
- After computing pooled variances, apply limma-style `squeezeVar()` for empirical Bayes shrinkage
- Pro: preserves the current "missing = below LOD" philosophy while gaining variance moderation
- Con: custom implementation, doesn't leverage existing model infrastructure

### Key Design Questions

1. Should we impute **individual values** (for lm/limma fitting) or **group means** (current approach)?
2. If imputing individual values, should imputed observations receive **reduced weight** in the model?
3. Should the LOD estimation change (e.g., per-protein LOD, more sophisticated detection limit models)?
4. Should the zero-out logic for one-sided missingness be preserved?
5. How to handle the `indic` column and missingness metadata in the new approach?



# New Imputation

I want an alternative lm based implementaton of ContrastImpute, I would be an extension of the contrast lm class.                                                                          

We build the models, as usual. Then for those proteins where we fail to build, we impute the NA with the lod. but also set all values below lod to lod. We then fit the model, however we do not use the convariance matrix, of the model fit. But we use a covariance matrix, which is a mean or median of the covariance matrices which we did estimate from the other data.

The reason being, since we impute with a constant value we would underestimate the variance, therefore we take the variance from the other models.

The idea of using the lod not just replace the NA but also as minimum value is to avoid situations that a group where we observed the data will get a higher estimate that a group where we failed to obtain an estimate.

Last but not least degrees of freedom, ideally only as many as indeed observed values.

So basic principle impute -> fit model -> update model to reflect missingness and observed variances. What do you think (no coding here. we discuss)



Interesting approach. Let me think through the mechanics and flag some considerations.

### What I like

**Reusing the lm infrastructure** — you get contrast matrices, multi-factor designs, and the existing `Contrasts` class machinery for free. No manual group-mean arithmetic.

**Clamping to LOD (not just imputing NAs)** — this solves the current zero-out logic problem elegantly. If group A has one observation at 8.5 and group B is all-missing (imputed at LOD=10), the current code gives a "wrong direction" fold change and zeroes it out. With clamping, group A's value becomes max(8.5, LOD) and the contrast naturally behaves.

**Borrowing the covariance matrix** — this is the key insight. The imputed fit gives you correct *point estimates* (coefficients), but the variance is artificial. Swapping in a typical covariance matrix from successful fits restores honest uncertainty.

### Questions and concerns

**1. Which covariance matrix to borrow?**

You said mean or median of covariance matrices from successful fits. But covariance matrices aren't scalars — you need to aggregate them sensibly:

- **Option A:** Take the median of `sigma^2` (residual variance) across successful fits, then scale the `(X'X)^{-1}` from the imputed fit by that median sigma. This is simple and preserves the design structure.
- **Option B:** Take the full `vcov()` matrices and compute element-wise median. This is messier and the matrices may not be positive definite after aggregation.

Option A seems cleaner — you're essentially saying "the coefficients are from the imputed fit, but the residual variance is borrowed." This is conceptually close to what limma's `squeezeVar` does.

**2. Degrees of freedom — which df?**

You want df to reflect actual observed values. Two choices:

- **Residual df of the imputed model** (n_total - p) — but this overcounts since imputed values aren't real observations
- **n_observed - p** — honest but could be very small or even zero (which is why the model failed in the first place)

If n_observed < p, you have no real df. You could use the df from the borrowed variance estimate (the successful fits), which gives you a "prior df" similar to limma's approach.

**3. When exactly do you impute?**

You said: fit normally first, then impute for failures. But the `strategy_lm()` / `build_model()` path fits per-protein. So:

- Proteins with enough data → normal fit, normal covariance → use as-is (these are "donors")
- Proteins with too many NAs → model fails → impute with LOD + clamp → refit → replace covariance with donor estimate

This creates two tiers of results. The interface should probably flag which tier each protein came from (similar to the current `indic` column).

COMMENT: 

**4. What does "model fails" mean exactly?**

- `lm()` returns NA coefficients? 
- Singular fit?
- Or a threshold like "more than X% missing"?

You might also want to impute for proteins that *technically* fit but have, say, only 1 observation in one group — the fit succeeds but the variance estimate is garbage.

**5. Relationship to existing classes**

Architecturally, this could be:

- A **post-processor on `Model`** — fit all proteins (imputing where needed), then selectively replace `vcov` on the model objects
- A **new `Contrasts*` class** that wraps `Contrasts` — takes a normal `Model`, identifies failures, refits with imputation, borrows covariance, and merges results
- A **facade** like `ContrastsRLMFacade` that orchestrates: fit → identify failures → impute → refit → patch covariance → merge

The facade pattern seems most natural given how the other "fallback" approaches work (e.g., `ContrastsMissing` is already used as a fallback in `prolfquapp`).

### Summary of the design as I understand it

```
For each protein:
  1. Try lm fit on observed data
  2. If success → keep fit, contribute vcov to "donor pool"
  3. If failure:
     a. Replace NAs with LOD, clamp all values to max(value, LOD)
     b. Fit lm on imputed data → get coefficients
     c. Replace vcov with donor estimate (median sigma² × (X'X)⁻¹)
     d. Set df = n_observed - p (or donor df if too small)
  4. Extract contrasts using standard Contrasts machinery
```

Does this match your thinking? The main open question from my side is whether option A (borrow scalar sigma², keep imputed design's `(X'X)^{-1}`) is what you have in mind for the covariance replacement.