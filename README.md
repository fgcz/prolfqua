# LFQService - for Label Free Quantification Services

![prolfqua](inst/Figures/imgfile.png)

The R package contains vignettes and functions for analyzing mass spec device based LFQ proteomics experiments at the [FGCZ](http://www.fgcz.ch/).


## 1. System Requirements  

## 1.1. Install R

a Windows/Linux/MacOSX x64 platform R 3.6.

```{r}
# requires install.packages(c(BiocManager, 'remotes'))
BiocManager::install('wolski/LFQService')
```

### 1.2. Required packages

```{r}
pkgs <- c('broom', 'bookdown', 'conflicted', 'corrplot', 'dplyr',
  'flextable', 'docopt',
  'ggplot2', 'ggbeeswarm', 'ggfortify', 'glue', 'GGally', 'pheatmap',
  'kableExtra', 'limma', 'lme4', 'lmerTest', 'magrittr', 'multcomp',
  'plotly', 'purrr', 'readxl', 'tidyverse', 'yaml',
  'tidyr', 'writexl')

pkgs <- pkgs[(!pkgs %in% unique(installed.packages()[,'Package']))]
if(length(pkgs) > 0){install.packages(pkgs)}
```



## 2. Running R-scripts

Generate QC report from bat file.
First add `<LFQService_path>/win` to your path variable. Then you can generate a QC report from a maxquant QC by running.


```
lfq_MQ_SampleSizeReport.bat .\data\1296877_QC.zip
```


## 3. Best of code snippets


```
Rscript ~/__checkouts/R/LFQService/inst/run_scripts/lfq_MQ_SampleSizeReport.R --help
Rscript ~/__checkouts/R/LFQService/inst/run_scripts/lfq_MQ_SampleSizeReport.R ~/Downloads/1330043.zip
```


## Motivation

The package for _pro_teomics _l_abel _f_ree _qua_ntification `prolfqua` (read : prolevka) evolved from a set of scripts and functions written in the R programming language to visualize and analyze mass spectrometric data, and some of them are still in R packages such as quantable, protViz or imsbInfer. For computing protein fold changes among treatment conditions, we first used t-test or linear models, then started to use functions implemented in the package limma to obtain moderated p-values. We did also try to use other packages such as MSStats, ROPECA or MSqRob all implemented in R, with the idea to integrate the various approaches to protein fold-change estimation. Although all these packages were written in R,  model specification, input and output formats differ widely and wildly, which made our aim to use the original implementations challenging. Therefore, and also to understand the algorithms used, we attempted to reimplement those methods, if possible. 

When developing `prolfqua` we were inspired by packages such as `sf` or `stars` which use data in the long table format, use dplyr for data transformation and ggplot2 for visualization. In `prolfqua` the data needed for analysis is represented using a single data.frame and a configuration object. The configuration annotates the table, specifies what information is in which column. The results of the statistical modelling we also store in data-frames, making filtering, transforming the results with `dplyr` and `ggplot2` easy.

The use of an annotated table makes integrating new data if provided in long formatted tables simple. Hence for Spectronaut or Skyline outputs, or MSStats formatted input, all is needed is a table annotation (see code snipped). For software, which writes the data in a wide tables, e.g. Maxquant, we implemented methods which first transform the data into a long format.  Relying on the long data table format enabled us to apply a large variety of useful visualizations and data preprocessing methods in R. 

A further design decision, which differentiates `prolfqua` is that it embraces and supports R's linear model formula interface, and R's lme4 mixed effect models formula interface. R's formula interface for regression models is flexible, widely used and documented. The linear model and linear mixed model interfaces allow specifying a wide range of essential models, including parallel designs, factorial designs, repeated measurements and many more. Since `prolfqua` uses R modelling infrastructure, we can fit all these models to proteomics data.

This is not possible with any other package dedicated to proteomics data analysis. For instance, MSStats, although it uses the same modelling infrastructure, supports only a small subset of possible models. Limma, on the other hand, supports the R formula interface but not for linear mixed models. Because ROPECA relies on `limma` it is limited to the same subset of models. MSqRob is limited to random effects model's only, and it remains unclear how to fit these models to factorial designs, and how interactions among factors can be computed and tested. MSqRob was benchmarked only on parallel-group designs.

The use of R's formula interface does not limit `prolfqua` to outputs provided by R modelling infrastructure. `prolefqa` also implements p-Value moderation, as in the limma publication or computing probabilities of differential regulation as suggested in the ROPECA publication. 
Moreover, the design decision allowed us to integrate Bayesian regression models provided by the r-package `brms` easily. Because of that, we can benchmark all those methods: linear models, mixed effect models, p-value moderation, ROPECA as well as Bayesian regression models within the same framework. This helped us to deepen our understanding of the practical application of these methods.

Last but not least `prolfqua` supports the LFQ data analysis workflow, e.g. computing coefficients of Variations (CV) for peptide and proteins, it can be used for sample size estimation, visualization and summarization of missing data and intensity distributions, multivariate analysis of the data, etc.
It also implements various protein intensity summarization and inference methods, e.g. Top 3, or Tukeys median polish. Last but not least, ANOVA analysis or model selection using the likelihood ratio test for thousand of proteins is implemented and can be performed.


To use `prolfqua` knowledge of the R regression model infrastructure is of advantage. Acknowledging, the complexity of the formula interface, we provide an  MSstats emulator, where the model formula is inferred from the annotation file structure of an MSStats analysis.


# How to get started

## 2 grp analysis
- See vignettes.

## Parallel design
What is parallel design.
- See vignettes.

## repeated measurements
What is a paired measurement.
- See vignettes.

## Contrasts
- See vignettes.


# TODO
- Make public

