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

prolfqua (prolevka) evolved from a set of scripts and functions written in the R programming languate to visualize and analyse mass spectrometric data, and some of them are still in R packages such as quantable, protViz or imsbInfer. For computing protein fold changes among treatment conditions, we first used linear models, than started to use functions implemented in the package limma to obtain moderated p-values. We did also try to use other packages such as MSStats, ROPECA or MSqRob all implemented in R, with the idea to integrate the various approaches to protein fold change estimation. However, although all these packages are implmented in R the way how the data is handled, the models are specified, and how the results are presented differce among the those packages widely, making the task of using the original implementations challenging. Because of this, and also to better understand the algorithms used we tried to reimplmement those methods, if posssible. One of the the design decisions made when developing `prolfqua` was inspired by packages such as `sf` or `stars` which are tightly integrated with dplyr and for data managment use tables in long format. In `prolfqua` the data needed for analysis is represented using a single data.frame and an configuration object which annotates the table columns. Also the results of the statistical modelling are represented in data.frames, making filtering, transforming the results with dplyr and ggplot2 easy.

Integrating new data if provided long format, which is the case for many software tools, e.g. spectronaut or skyline peptide or protein reports, is simple. For software, which writes the data in a long format, e.g. Maxquant we implemented methods which first transforms them into long format.  



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

