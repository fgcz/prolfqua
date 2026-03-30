# The Limpa Pipeline: A Modular Walkthrough

## Architecture

The limpa pipeline has four modules connected by well-defined data objects. Each module can be understood, tested, and (in some cases) swapped independently.

```
┌─────────────────────┐      ┌──────────────────────┐      ┌──────────────────┐      ┌──────────────────────┐
│  Module 1: Quantify  │ ──▶  │  Module 2: Fit Model  │ ──▶  │ Module 3: Contrast│ ──▶  │ Module 4: Moderate   │
│                      │      │                       │      │                   │      │                      │
│  dpcQuant()          │      │  voomaLmFit           │      │  contrasts.fit()  │      │  eBayes()            │
│  dpcQuantByRow()     │      │  WithImputation()     │      │                   │      │  topTable()          │
│  (or MaxLFQ, etc.)   │      │                       │      │                   │      │  treat()             │
└─────────────────────┘      └──────────────────────┘      └──────────────────┘      └──────────────────────┘
        │                              │                            │                          │
   EList with $se              MArrayLM (fit)              MArrayLM (fit2)           moderated MArrayLM
```

### The Interface Objects

| After module | Object type | Key contents |
|-------------|-------------|-------------|
| **1 → 2** | `EList` | `$E` (expression matrix, no NAs), `$se` (standard errors), `$genes$Imputed` (logical matrix) |
| **2 → 3** | `MArrayLM` | `$coefficients`, `$stdev.unscaled`, `$sigma`, `$df.residual`, `$cov.coefficients` |
| **3 → 4** | `MArrayLM` | Same, but coefficients now represent contrasts rather than design columns |
| **4 → user** | `MArrayLM` | Adds `$t` (moderated t), `$p.value`, `$s2.post`, `$df.total` |

---

## Module 1: Quantify Proteins

### Purpose

Convert a precursor-level matrix (precursors × samples, with NAs) into a protein-level matrix (proteins × samples, **no NAs**) with accompanying **standard errors** that measure quantification uncertainty per observation.

### Entry Points

There are two functions, depending on the input structure:

**`dpcQuant(y, protein.id)`** — the standard path. Multiple precursors per protein are summarised via an additive model (sample effects + precursor effects) with DPC-aware posterior estimation.

**`dpcQuantByRow(y)`** — the simplified path for data that is already one row per protein (or one precursor per protein). No precursor summarisation needed; performs DPC-aware maximum posterior estimation row by row. If `dpcQuant()` detects all proteins have exactly one precursor, it delegates here automatically.

Note: `dpcImpute()` is the deprecated name for `dpcQuantByRow()`.

### What Happens Inside

#### 1a. DPC Estimation

If a DPC is not provided, the function estimates it internally:

```r
beta0 <- estimateDPCIntercept(y, dpc.slope = dpc.slope)   # default slope = 0.8
dpc <- c(beta0 = beta0, beta1 = dpc.slope)
```

The DPC models detection probability as a logistic function of intensity: $P(d=1 \mid y) = \text{logit}^{-1}(\alpha_0 + \alpha_1 y)$. It can also be pre-estimated via `dpcCN()` for the complete-normal model and passed in.

#### 1b. Hyperparameter Estimation

Empirical Bayes hyperparameters are estimated from the global data:

```r
h <- dpcQuantHyperparam(y, protein.id = protein.id, dpc.slope = dpc[2])
```

This yields three quantities governing the prior on sample effects $\gamma_i$:

| Hyperparameter | Meaning | Code |
|---------------|---------|------|
| $\mu_\gamma$ | Global mean log-expression | `h$prior.mean` |
| $\sigma_{\text{mean}}$ | SD of mean expression across proteins | `h$prior.sd` |
| $\sigma_{\text{logFC}}$ | SD of within-protein variation (i.e., log-fold-changes) | `h$prior.logFC` |

Both `prior.sd` and `prior.logFC` are floored at 1 to prevent degenerate priors.

#### 1c. Protein Summarisation

For each protein, `peptides2Proteins()` minimises:

$$Q = D_{\text{obs}} + D_{\text{mis}} + D_{\text{prior}}$$

- $D_{\text{obs}}$: squared residuals from observed precursor values
- $D_{\text{mis}}$: log-probability of being missing under the DPC (integrals via Gauss-Hermite quadrature)
- $D_{\text{prior}}$: empirical Bayes penalty pulling extreme estimates towards the global mean

The flag `standard.errors = TRUE` triggers computation of posterior standard errors — this is what makes Module 2's bivariate variance model possible.

### Output

An `EList` object:

```r
y.protein$E               # G × n protein-level log2-intensities (no NAs)
y.protein$se              # G × n posterior standard errors
y.protein$genes$Imputed   # G × n logical: TRUE where all precursors were missing
y.protein$dpc             # c(alpha0, alpha1)
```

### Alternative: Using External Quantification

Module 2 does not require DPC-Quant. If you have protein-level data from MaxLFQ, MSstats, or any other method, you can use `voomaLmFitWithImputation()` with `predictor = NULL` and it reduces to standard vooma. If you have some per-protein quality metric (e.g., PSM counts), you can pass that as the `predictor` argument instead.

---

## Module 2: Fit the Linear Model

### Purpose

Estimate per-protein regression coefficients under a user-specified experimental design, using precision weights derived from a bivariate variance model that accounts for both the mean-variance trend and quantification uncertainty.

### Function

```r
fit <- voomaLmFitWithImputation(
    y         = y.protein,           # EList from Module 1
    design    = design,              # design matrix
    predictor = log(y.protein$se),   # log-SE matrix from DPC-Quant
    imputed   = y.protein$genes$Imputed,  # imputation indicator
    block     = NULL,                # optional blocking factor
    plot      = TRUE                 # show variance trend
)
```

The convenience wrapper `dpcDE()` does exactly this — it extracts `$se` and `$genes$Imputed` from a `dpcQuant()` output and calls `voomaLmFitWithImputation()`.

### What Happens Inside

#### 2a. Preliminary OLS Fit (line 74)

```r
fit <- lmFit(y, design, weights = prior.weights)
```

Produces per-protein residual SDs (`fit$sigma`), fitted values, and average log-expression. Identical to vooma.

#### 2b. df Correction for Entirely Imputed Groups (lines 84–120)

Unique to limpa. When all precursors for a protein are missing in one experimental group, the DPC-Quant value for that group is driven entirely by the prior, not by data. The code detects this via the `imputed` matrix and the design's hat matrix:

```r
h <- 1 - hat(design)
DFImputed <- imputed %*% (1 - h)
```

For affected proteins, it re-fits with imputed values set to `NA`, obtaining corrected `sigma` and `df.residual`. This prevents artificially deflated variance estimates.

#### 2c. Bivariate Variance Model (lines 122–145)

When `predictor` is provided, six lines transform univariate vooma into bivariate limpa:

```r
sx  <- A                                         # mean log-expression per protein (length G)
sy  <- sqrt(fit$sigma)                           # (variance)^(1/4) per protein (length G)
mu  <- fitted.values                             # G × n matrix

sxc <- rowMeans(predictor, na.rm = TRUE)         # [1] mean log-uncertainty per protein
vartrend <- lm.fit(cbind(1, sx, sxc), sy)        # [2] bivariate linear model
beta <- coef(vartrend)                            # [3] (beta_0, beta_1, beta_2)
sx <- vartrend$fitted.values                      # [4] overwrite: protein-level 1D index
mu <- beta[1] + beta[2]*mu + beta[3]*predictor    # [5] overwrite: observation-level 1D index
```

**The single-index model trick:** The linear model projects from 2D (mean expression, log-uncertainty) to 1D. This index is what lowess will smooth against, avoiding a 2D surface smoother.

When `predictor` is `NULL`, this block is skipped and the function reduces to standard vooma.

#### 2d. Lowess + Weights (lines 154–176)

```r
l <- lowess(sx, sy, f = span)                    # smooth in 1D index space
f <- approxfun(l, rule = 2, ties = list("ordered", mean))
w <- 1/f(mu)^4                                   # per-observation precision weights
```

Identical machinery to vooma — but `sx` and `mu` have been swapped out by Step 2c if a predictor was provided.

#### 2e. Optional Iteration (lines 186–275)

If `block` or `sample.weights` are specified:

1. Estimate sample weights / intra-block correlation with current vooma weights
2. Re-fit preliminary model → re-estimate bivariate trend → re-compute weights
3. Re-estimate sample weights / correlation

Two iterations total (initial + one refinement).

#### 2f. Final WLS Fit (line 278)

```r
fit <- lmFit(y, design, block = block, correlation = correlation, weights = weights)
```

Per-protein WLS: $\hat{\beta}_g = (X^T W_g X)^{-1} X^T W_g y_g$

### Output

An `MArrayLM` object containing:

| Field | Shape | Description |
|-------|-------|-------------|
| `$coefficients` | $G \times p$ | $\hat{\beta}_g$ for each protein (one column per design parameter) |
| `$stdev.unscaled` | $G \times p$ | $\sqrt{\text{diag}((X^T W_g X)^{-1})}$ — the unscaled SE for each coefficient |
| `$sigma` | length $G$ | Residual SD per protein: $s_g = \sqrt{\frac{(y_g - X\hat{\beta}_g)^T W_g (y_g - X\hat{\beta}_g)}{n - p}}$ |
| `$df.residual` | length $G$ | Residual df per protein (corrected for imputed groups) |
| `$cov.coefficients` | $p \times p$ | $(X^T W_g X)^{-1}$ — **not** multiplied by $\sigma^2_g$ yet |

**Key relationship:** The variance of $\hat{\beta}_{gk}$ (coefficient $k$ for protein $g$) is:

$$\text{Var}(\hat{\beta}_{gk}) = \sigma^2_g \cdot \text{stdev.unscaled}_{gk}^2 = \sigma^2_g \cdot [(X^T W_g X)^{-1}]_{kk}$$

This factorisation into $\sigma^2_g$ (protein-specific) and $(X^T W_g X)^{-1}$ (shared structure) is what enables `eBayes()` to moderate across proteins.

---

## Module 3: Compute Contrasts

### Purpose

Re-parameterise the fitted model so that the coefficients represent the **specific comparisons** of interest, rather than the raw design matrix columns.

### When You Need This

If your design matrix directly encodes the comparisons you care about (e.g., a single treatment column in an intercept + treatment design), you can skip this module — the coefficient of interest is already a column of `fit$coefficients`.

You need `contrasts.fit()` when:

- Your design uses group means (one column per group) and you want pairwise differences
- You want to test linear combinations like (A + B)/2 − (C + D)/2
- You want to test multiple contrasts simultaneously

### How It Works

```r
# Define contrast matrix
contrast.matrix <- makeContrasts(
    TreatmentVsControl = Treatment - Control,
    DrugAVsDrugB       = DrugA - DrugB,
    levels = design
)

# Apply contrasts
fit2 <- contrasts.fit(fit, contrasts = contrast.matrix)
```

#### What `contrasts.fit()` Does Mathematically

If $C$ is the contrast matrix ($p \times k$, where $k$ is the number of contrasts), then for each protein:

$$\hat{\beta}^*_g = \hat{\beta}_g \, C$$

$$\text{Cov}(\hat{\beta}^*_g) = \sigma^2_g \, C^T (X^T W_g X)^{-1} C$$

The new unscaled standard errors are:

$$\text{stdev.unscaled}^*_{gk} = \sqrt{[C^T (X^T W_g X)^{-1} C]_{kk}}$$

**Important:** `contrasts.fit()` is a pure linear algebra operation. It does not re-estimate anything — it simply re-projects the existing coefficient estimates and their covariance structure.

### Output

An `MArrayLM` object with the same structure as Module 2's output, but:

- `$coefficients` now has $k$ columns (one per contrast) instead of $p$
- `$stdev.unscaled` and `$cov.coefficients` are updated accordingly
- `$sigma` and `$df.residual` are unchanged — they are properties of the residuals, which contrasts don't alter

---

## Module 4: Moderate and Extract Results

### Purpose

Apply empirical Bayes shrinkage of per-protein variance estimates towards a common prior, producing moderated t-statistics with better statistical properties (more stable, better calibrated) than ordinary per-protein t-tests. Then extract ranked results.

### `eBayes()`: Empirical Bayes Moderation

```r
fit2 <- eBayes(fit2, trend = FALSE)
```

Note: the paper uses `trend = TRUE` for competing limma-based pipelines (the "limma-trend" approach), but when `voomaLmFitWithImputation()` has already produced vooma-style weights, the mean-variance trend is already accounted for in the weights, so `trend = FALSE` is the default for the limpa pipeline.

#### What `eBayes()` Does

##### 4a. Fit the Prior Distribution

The per-protein residual variances $s^2_g$ are assumed to follow a scaled inverse-chi-squared distribution:

$$\frac{1}{s^2_g} \sim \frac{1}{s^2_0} \cdot \frac{\chi^2_{d_0}}{d_0}$$

where $s^2_0$ is the prior variance and $d_0$ is the prior degrees of freedom. These are estimated from the observed distribution of $s^2_g$ across all proteins by fitting a scaled F-distribution to the sample variances (via `fitFDist()`).

##### 4b. Compute Posterior (Moderated) Variances

For each protein, the moderated variance is:

$$\tilde{s}^2_g = \frac{d_0 \, s^2_0 + d_g \, s^2_g}{d_0 + d_g}$$

where $d_g$ = `fit$df.residual[g]`. This is a weighted average of the prior variance $s^2_0$ and the protein's own residual variance $s^2_g$, with weights proportional to their respective degrees of freedom.

The posterior degrees of freedom are:

$$\tilde{d}_g = d_0 + d_g$$

##### 4c. Compute Moderated t-Statistics

For contrast $k$ and protein $g$:

$$\tilde{t}_{gk} = \frac{\hat{\beta}^*_{gk}}{\tilde{s}_g \cdot \text{stdev.unscaled}^*_{gk}}$$

This replaces the ordinary t-statistic (which uses $s_g$) with one using the moderated $\tilde{s}_g$. The moderated t follows a t-distribution on $\tilde{d}_g$ degrees of freedom under the null.

**Why this helps:** Proteins with very few observations (or many imputed values, hence low df) get their variance estimates pulled towards the global trend. This prevents both false positives (from proteins with artificially low variance) and false negatives (from proteins with artificially high variance due to noise).

##### 4d. Compute p-Values and B-Statistics

P-values are computed from the moderated t-distribution:

$$p_{gk} = 2 \cdot P(t_{\tilde{d}_g} > |\tilde{t}_{gk}|)$$

The B-statistic (log-odds of differential expression) is also computed for each protein.

### Output

The `MArrayLM` object is updated with:

| Field | Shape | Description |
|-------|-------|-------------|
| `$t` | $G \times k$ | Moderated t-statistics |
| `$p.value` | $G \times k$ | P-values from moderated t |
| `$s2.post` | length $G$ | Posterior (moderated) variance $\tilde{s}^2_g$ |
| `$s2.prior` | scalar | Prior variance $s^2_0$ |
| `$df.prior` | scalar | Prior degrees of freedom $d_0$ |
| `$df.total` | length $G$ | $\tilde{d}_g = d_0 + d_g$ |
| `$lods` | $G \times k$ | B-statistic (log-odds of DE) |

### `topTable()`: Extract Ranked Results

```r
results <- topTable(fit2, coef = "TreatmentVsControl", number = Inf, sort.by = "P")
```

Returns a data frame with columns:

| Column | Source |
|--------|--------|
| `logFC` | `fit2$coefficients[, coef]` |
| `AveExpr` | `fit2$Amean` (average log2-expression) |
| `t` | `fit2$t[, coef]` (moderated t-statistic) |
| `P.Value` | `fit2$p.value[, coef]` |
| `adj.P.Val` | Benjamini-Hochberg adjusted p-values |
| `B` | `fit2$lods[, coef]` (log-odds of DE) |

### `treat()`: Testing with a Fold-Change Threshold

An alternative to `eBayes()` when you want to test whether the fold-change exceeds a minimum threshold, rather than just testing whether it's non-zero:

```r
fit2 <- treat(fit2, lfc = log2(1.5))   # test |logFC| > log2(1.5)
results <- topTreat(fit2, coef = "TreatmentVsControl", number = Inf)
```

`treat()` applies the same empirical Bayes moderation but tests $H_0: |\beta_{gk}| \leq \text{lfc}$ instead of $H_0: \beta_{gk} = 0$.

---

## Complete Modular Pipeline in Code

### Minimal: Using `dpcDE()` Wrapper

```r
library(limpa)

# Module 1: Quantify
y.protein <- dpcQuant(y_precursor, protein.id = "Protein.Group")

# Modules 2 (fit): dpcDE wraps voomaLmFitWithImputation
fit <- dpcDE(y.protein, design)

# Module 3: Contrasts (skip if coefficient of interest is a design column)
fit2 <- contrasts.fit(fit, contrasts = makeContrasts(B - A, levels = design))

# Module 4: Moderate + extract
fit2 <- eBayes(fit2)
results <- topTable(fit2, coef = 1, number = Inf)
```

### Full Manual Pipeline

```r
library(limpa)
library(limma)

# ─── Module 1: Quantify ──────────────────────────────────────────────
# Option A: DPC-Quant (multi-precursor → protein)
dpc <- dpcCN(y_precursor)
y.protein <- dpcQuant(y_precursor, protein.id = "Protein.Group", dpc = dpc)

# Option B: DPC-QuantByRow (already one row per protein)
# y.protein <- dpcQuantByRow(y_protein_level)

# Option C: External quantification (MaxLFQ, MSstats, etc.)
# y.protein <- ... (no $se available, predictor = NULL in Module 2)


# ─── Module 2: Fit Model ─────────────────────────────────────────────
predictor <- log(y.protein$se)             # log-SE matrix from DPC-Quant
imputed   <- y.protein$genes$Imputed       # imputation indicator

fit <- voomaLmFitWithImputation(
    y         = y.protein,
    design    = design,
    predictor = predictor,
    imputed   = imputed,
    block     = patient_id,                # optional: triggers iteration
    plot      = TRUE
)


# ─── Module 3: Contrasts ─────────────────────────────────────────────
# Example 1: Simple pairwise
contrast.matrix <- makeContrasts(
    TEN_vs_MPR = TEN - MPR,
    levels = design
)
fit2 <- contrasts.fit(fit, contrasts = contrast.matrix)

# Example 2: Interaction / average effects
# contrast.matrix <- makeContrasts(
#     AvgTreatment = (DrugA + DrugB)/2 - Control,
#     DrugDiff     = DrugA - DrugB,
#     levels = design
# )
# fit2 <- contrasts.fit(fit, contrasts = contrast.matrix)

# Example 3: No contrasts needed (coefficient already in design)
# fit2 <- fit   # just pass through


# ─── Module 4: Moderate + Extract ────────────────────────────────────
fit2 <- eBayes(fit2)
results <- topTable(fit2, coef = "TEN_vs_MPR", number = Inf)

# Alternative: fold-change threshold testing
# fit2 <- treat(fit2, lfc = log2(1.5))
# results <- topTreat(fit2, coef = "TEN_vs_MPR", number = Inf)

# Gene set testing (works on the same fit object)
# camera(y.protein, index = gene_sets, design = design, contrast = contrast.matrix[, 1])
# fry(y.protein, index = gene_sets, design = design, contrast = contrast.matrix[, 1])
```

---

## Data Flow Summary

```
 Precursor matrix                    Protein EList                   MArrayLM (raw)
 (precursors × n)                    (proteins × n)                 (proteins × p coefficients)
 ┌──────────────┐                   ┌──────────────┐               ┌──────────────────┐
 │ y_ji with NAs │ ── Module 1 ──▶  │ $E  (no NAs) │ ── Module 2 ─▶│ $coefficients     │
 │               │    dpcQuant()    │ $se           │   voomaLmFit │ $stdev.unscaled   │
 │               │                  │ $Imputed      │   WithImp()  │ $sigma             │
 └──────────────┘                   └──────────────┘               │ $df.residual       │
                                                                    └────────┬─────────┘
                                                                             │
                                              ┌──────────────────────────────┘
                                              │
                                              ▼
                                    MArrayLM (contrasted)            MArrayLM (moderated)
                                    (proteins × k contrasts)         (proteins × k contrasts)
                                    ┌──────────────────┐            ┌──────────────────┐
                              ──────│ $coefficients     │── Mod 4 ──│ $t (moderated)    │
                   Module 3         │ $stdev.unscaled   │  eBayes() │ $p.value          │
                   contrasts.fit()  │ $sigma (unchanged)│           │ $s2.post          │
                                    │ $df.res (same)    │           │ $df.total          │
                                    └──────────────────┘            └──────────────────┘
                                                                             │
                                                                             ▼
                                                                     topTable() / topTreat()
                                                                    ┌──────────────────┐
                                                                    │ logFC, AveExpr,   │
                                                                    │ t, P.Value,       │
                                                                    │ adj.P.Val, B      │
                                                                    └──────────────────┘
```

---

## Swappability Guide

| Module | Limpa default | Can swap for | What you lose |
|--------|--------------|--------------|---------------|
| **1** | `dpcQuant()` | MaxLFQ + knn/Perseus imputation | Standard errors → no bivariate trend, no df correction |
| **1** | `dpcQuant()` | `dpcQuantByRow()` | Multi-precursor summarisation (OK if already 1-row-per-protein) |
| **2** | `voomaLmFitWithImputation()` with predictor | Same function, `predictor = NULL` | Bivariate trend → univariate vooma only |
| **2** | `voomaLmFitWithImputation()` | `vooma()` + `lmFit()` | df correction, imputed-group handling, bivariate trend |
| **2** | `voomaLmFitWithImputation()` | Plain `lmFit()` + `eBayes(trend=TRUE)` | All variance modelling (limma-trend approach) |
| **3** | `contrasts.fit()` | Direct coefficient extraction | Ability to test arbitrary linear combinations |
| **4** | `eBayes()` | `treat()` | Tests non-zero logFC instead of exceeding a threshold |
| **4** | `eBayes()` | Nothing (mandatory) | Cannot skip — moderation is what makes the pipeline work |

---

## Side-by-Side: vooma + limma vs. limpa (Modular View)

| Module | vooma + limma | limpa |
|--------|--------------|-------|
| **1. Quantify** | External: MaxLFQ + imputation | `dpcQuant()`: additive model + DPC + EB prior |
| **1. Output** | Expression matrix only | Expression matrix + SE matrix + imputed indicator |
| **2. Trend** | 1D lowess on mean expression | `lm.fit` on (mean expr, log-SE) → 1D lowess on projected index |
| **2. df correction** | None | Detects entirely-imputed groups, corrects `sigma` and `df.residual` |
| **2. Weights** | $w_{gi} = 1/f(\hat{\mu}^*_{gi})^4$ | $w_{gi} = 1/f(\hat{\beta}_0 + \hat{\beta}_1 \hat{\mu}^*_{gi} + \hat{\beta}_2 \psi_{gi})^4$ |
| **3. Contrasts** | `contrasts.fit()` | `contrasts.fit()` (identical) |
| **4. Moderate** | `eBayes()` | `eBayes()` (identical) |
