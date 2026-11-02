# prolfqua: Proteomics Label Free Quantification

prolfqua covers the steps of a mass spectrometry based label-free
quantification (LFQ) analysis: quality control, normalization,
aggregation of peptides to proteins, modelling, hypothesis testing and
sample size estimation. A single function, build_contrast_analysis(),
fits a model per protein and computes contrasts with one of many
methods: linear, robust and rank-based regression, limma and limma-voom,
DEqMS variance moderation, limpa's detection probability model, Firth
logistic and binomial models of detection, mixed models and ROPECA on
peptide-level data, a hurdle model combining intensity and detection,
and variants that impute groups without observations at the limit of
detection. All methods return their results through the same interface
and columns, so they can be compared and combined. prolfqua models
simple designs with one explanatory variable as well as experiments with
several factors, and was developed while analysing dozens of proteomics
datasets at the Functional Genomics Center Zurich.

## See also

Useful links:

- <https://github.com/fgcz/prolfqua>

- <https://fgcz.github.io/prolfqua/>

- Report bugs at <https://github.com/fgcz/prolfqua/issues>

## Author

**Maintainer**: Christian Panse <cp@fgcz.ethz.ch>
([ORCID](https://orcid.org/0000-0003-1975-3064))

Authors:

- Witold Wolski <witold.wolski@fgcz.uzh.ch>
  ([ORCID](https://orcid.org/0000-0002-6468-120X))

Other contributors:

- Functional Genomics Center Zurich \[funder\]
