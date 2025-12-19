# CLAUDE.md

This file provides guidance to Claude Code when working with the prolfqua R package.

## Overview

prolfqua is an R package for label-free quantification (LFQ) proteomics data analysis. It provides tools for:
- Data preprocessing and normalization
- Statistical modeling and differential expression analysis
- Visualization and quality control

## Development Conventions

### NAMESPACE Management
- **NAMESPACE is managed by roxygen2**. Never edit it directly.
- Use `@export` tags in R code to export functions/classes.
- Run `devtools::document()` to update NAMESPACE and man pages.

### Testing
- **Test via `@examples` first**. Roxygen2 examples are the primary test mechanism.
- Temporary testthat tests can be created during development, but examples are preferred.
- Run examples with `devtools::run_examples()` or `R CMD check`.

### Documentation
- Use roxygen2 comments (`#'`) for all exported functions and classes.
- Run `devtools::document()` to generate man pages.

## Common Development Commands

```bash
# Install from source
R CMD INSTALL .

# Generate documentation
Rscript -e "devtools::document()"

# Run examples
Rscript -e "devtools::run_examples()"

# Check package
R CMD check prolfqua_*.tar.gz

# Run specific examples from a file
Rscript -e "devtools::run_examples(run_donttest = TRUE)"
```

## Key Classes

### LFQData (R6)
Main data container storing:
- `data`: Long-format data.frame
- `config`: AnalysisConfiguration object
- `prefix`: Output prefix

Key methods: `to_wide()`, `get_Transformer()`, `get_Aggregator()`, `get_Plotter()`

### LFQDataSE (R6)
SummarizedExperiment-backed alternative to LFQData. Provides:
- Same interface as LFQData
- Internal SE storage for Bioconductor interoperability
- Performance benefits for matrix operations

Use `as_LFQDataSE()` to convert from LFQData.

### Configuration Classes
- `AnalysisConfiguration`: Wraps table annotation and parameters
- `AnalysisTableAnnotation`: Column mappings (hierarchy, factors, response)
- `AnalysisParameters`: Analysis settings

## File Organization

```
R/
  LFQData.R              # Core data class
  LFQDataSE.R            # SE-backed variant
  LFQDataTransformer.R   # Data transformations
  LFQDataAggregator.R    # Peptide-to-protein aggregation
  LFQDataPlotter.R       # Visualization
  LFQDataStats.R         # Statistics
  LFQDataSummariser.R    # Data summaries
  AnalysisConfiguration.R # Configuration classes
```
