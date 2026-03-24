# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## What is prolfqua

An R package for mass spectrometry-based label-free quantification (LFQ)
proteomics analysis. It provides a complete workflow: QC, normalization,
protein aggregation, statistical modelling, hypothesis testing, and
sample size estimation. Data is always in long (tidy) format. Branch
`Modelling2R6` is the active development branch.

## Build & Test Commands

``` bash
make test            # Run testthat suite (runs document first)
make check-fast      # R CMD check without vignettes (quick validation)
make check           # Full R CMD check (document → build → check)
make document        # Generate roxygen2 docs (NAMESPACE + man/)
make install         # Install package locally
make lint            # Run lintr static analysis
make format          # Format with air
make build-vignettes # Build vignettes into inst/doc
make site            # Build pkgdown site locally
```

Single test file:

``` bash
Rscript -e "testthat::test_file('tests/testthat/test-LFQData.R')"
```

Library setup:

``` bash
Rscript -e ".libPaths()"
```

Use the normal user / system R libraries for this workspace; `renv`
autoload is disabled.

## Code Style

- **Line length:** 120 chars, **indentation:** 2 spaces (`.lintr`)
- `object_name_linter` is disabled — the codebase uses camelCase for R6
  classes and snake_case/mixed for functions
- NAMESPACE is **auto-generated** by roxygen2 — never edit directly; run
  `make document`
- Roxygen is configured with `r6 = TRUE` for R6 class documentation

## Architecture

### Core Data Flow

    Raw Data + AnalysisConfiguration → LFQData
        ├── get_Transformer() → LFQDataTransformer  (log2, robscale, normalize)
        ├── get_Aggregator()  → LFQDataAggregator   (peptide → protein rollup)
        ├── get_Stats()       → LFQDataStats         (CV, variance per group)
        ├── get_Plotter()     → LFQDataPlotter       (heatmaps, PCA, boxplots)
        ├── get_Summariser()  → LFQDataSummariser    (missingness, hierarchy counts)
        └── get_Imputer()     → LFQDataImp           (missing value imputation)

    LFQData → build_contrast_analysis(lfqdata, modelstr, contrasts, method)
        └── Returns a Facade with uniform API:
            $get_contrasts(), $get_missing(), $get_Plotter(), $to_wide()

### Facade Pattern (`ContrastsFacades.R`, `build_contrast_analysis.R`)

[`build_contrast_analysis()`](https://wolski.github.io/prolfqua/reference/build_contrast_analysis.md)
is the recommended entry point. Each method dispatches to a Facade class
that wires strategy → model → contrasts → moderation internally.

**Aggregated input** (protein-level, `subject_Id == hierarchy_keys`):
`lm`, `rlm`, `lm_missing`, `lm_impute`, `limma`, `deqms`, `firth`

**Nested input** (peptide-level, `subject_Id` is strict subset of
`hierarchy_keys`): `lmer`, `ropeca`

### Weights & `nr_children`

`config$nr_children` names the column tracking child-feature counts
(e.g. peptides per protein). After `get_Aggregator()` rollup, this
column is populated automatically. For peptide/precursor-level data it
is typically 1.

**Protein-level input must carry `nr_children`.** DEqMS uses it for
count-based variance modelling (`ContrastsDEqMSFacade`). If the column
is missing,
[`setup_analysis()`](https://wolski.github.io/prolfqua/reference/setup_analysis.md)
adds it set to 1 with a warning — but this defeats the purpose for
aggregated protein data where the actual peptide count matters.

### Key Design Patterns

**Decorator/Composition**: LFQData factory methods (`get_Transformer()`,
`get_Plotter()`, etc.) return decorator objects that wrap the LFQData.
Decorators hold a reference in their `lfq` field.

**Method chaining**: Transformer methods return `self` for chaining,
access result via `$lfq`:

``` r
lfqdata <- lfqdata$get_Transformer()$log2()$robscale()$lfq
```

**Strategy pattern for models**:
[`strategy_lm()`](https://wolski.github.io/prolfqua/reference/strategy.md),
[`strategy_lmer()`](https://wolski.github.io/prolfqua/reference/strategy.md),
[`strategy_rlm()`](https://wolski.github.io/prolfqua/reference/strategy.md),
[`strategy_glm()`](https://wolski.github.io/prolfqua/reference/strategy.md)
return lists with `model_fun`, `contrast_fun`, `isSingular`, `anova_df`.
[`strategy_limma()`](https://wolski.github.io/prolfqua/reference/strategy_limma.md)
returns a simpler list (formula, trend, robust) consumed by
[`build_model_limma()`](https://wolski.github.io/prolfqua/reference/build_model_limma.md).

**Config immutability**: AnalysisConfiguration is always deep-cloned
when passed to new LFQData instances. Never modify config in-place on an
existing LFQData.

### R6 Classes (22 classes across R/)

| Category            | Classes                                                                                                                                                                | Files                                                   |
|---------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------|
| Core data           | `LFQData`, `AnalysisConfiguration`                                                                                                                                     | LFQData.R, AnalysisConfiguration.R                      |
| Decorators          | `LFQDataTransformer`, `LFQDataAggregator`, `LFQDataStats`, `LFQDataPlotter`, `LFQDataSummariser`, `LFQDataImp`                                                         | LFQData\*.R                                             |
| Model interfaces    | `ModelInterface`, `Model`, `ModelFirth`, `ModelLimma`                                                                                                                  | Model\*.R, ContrastsLimma.R                             |
| Contrast interfaces | `ContrastsInterface`, `Contrasts`, `ContrastsModerated`, `ContrastsLimma`, `ContrastsProDA`, `ContrastsROPECA`, `ContrastsMissing`, `ContrastsFirth`, `ContrastsTable` | Contrasts\*.R, ContrastFirth.R, ContrastsSimpleImpute.R |
| Visualization       | `ContrastsPlotter`                                                                                                                                                     | ContrastsPlotter.R                                      |
| Utilities           | `MissingHelpers`                                                                                                                                                       | tidyMS_missigness_V2.R                                  |

### AnalysisConfiguration

Flat R6 class that maps column roles in the data: - **hierarchy**:
ordered measurement levels (protein_Id → peptide_Id → precursor_Id →
fragment_Id). `hierarchyDepth` controls which level is modelled. -
**factors**: explanatory variables (group, treatment). `factorDepth`
controls interaction depth. - **workIntensity**: response column. Uses a
stack (`set_response()` / `pop_response()` / `get_response()`) for
working with multiple intensity columns. - **fileName**: sample
identifier column.

Concrete config factories:
[`create_config_MQ_peptide()`](https://wolski.github.io/prolfqua/reference/concrete_AnalysisConfiguration.md),
`create_config_Skyline()`, `create_config_Spectronaut_Peptide()`, etc.
in `tidyMS_R6_ConcreteConfigurations.R`.

### Key Functions (not in classes)

- `build_contrast_analysis(lfqdata, modelstr, contrasts, method)` — main
  entry point, returns a Facade (in build_contrast_analysis.R)
- `setup_analysis(data, config)` — prepare data for analysis (in
  AnalysisConfiguration.R)
- `build_model(data, strategy, subject_Id)` — fit per-protein models (in
  tidyMS_R6Model.R)
- `build_model_impute(lfqdata, strategy)` — fit with LOD imputation +
  borrowed covariance for missing groups (in tidyMS_R6Model.R)
- `build_model_limma(lfqdata, strategy)` — fit limma matrix model (in
  ContrastsLimma.R)
- [`strategy_lm()`](https://wolski.github.io/prolfqua/reference/strategy.md),
  [`strategy_lmer()`](https://wolski.github.io/prolfqua/reference/strategy.md),
  [`strategy_rlm()`](https://wolski.github.io/prolfqua/reference/strategy.md),
  [`strategy_glm()`](https://wolski.github.io/prolfqua/reference/strategy.md)
  — per-protein model strategies (in tidyMS_R6_Modelling.R)
- [`strategy_limma()`](https://wolski.github.io/prolfqua/reference/strategy_limma.md)
  — limma matrix model strategy (in ContrastsLimma.R)
- [`sim_lfq_data_peptide_config()`](https://wolski.github.io/prolfqua/reference/sim_lfq_data_peptide_config.md)
  — simulate test data (in simulate_LFQ_data.R)

### File Naming Convention

- `R/LFQData*.R` — Core data container and its decorator classes
- `R/Model*.R`, `R/Contrasts*.R` — Modelling and hypothesis testing
- `R/AnalysisConfiguration.R` — Configuration system
- `R/tidyMS_R6_*.R` — Model strategies, concrete configs, correlation
  analysis
- `R/tidyMS_*.R` — Utility functions (plotting, stats, aggregation,
  missingness)

## Testing

7 test files in `tests/testthat/`: - `test-LFQData.R` — Core data
container and decorators - `test-Model.R` — Model fitting and
coefficient extraction - `test-Contrasts.R` — Contrast computation -
`test-ContrastsLimma.R` — Limma backend (ModelLimma, ContrastsLimma,
merge, 2-factor) - `test-ContrastsPlotter.R` — Contrast visualization -
`test-plotting_functions.R` — Low-level plots -
`test-tidyconfig_functions.R` — Configuration and utilities

## Cross-Package Context

prolfqua is part of the prolfqua ecosystem (see `../CLAUDE.md`).
Downstream packages depend on its R6 classes and exported API: -
**prolfquapp** — CLI wrapper for core facility workflows -
**prophosqua** — Phosphoproteomics analysis - **prolfquabenchmark** —
Benchmarking vignettes

Renaming R6 methods, changing exported function signatures, or modifying
AnalysisConfiguration fields can silently break these packages.
