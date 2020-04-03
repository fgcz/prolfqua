# LFQService - for Label Free Quantification Services


The R package contains vignettes and functions for analyzing mass spec device based LFQ proteomics experiments at the [FGCZ](http://www.fgcz.ch/).




## 1. System Requirements  

## 1.1. Install R

a Windows/Linux/MacOSX x64 platform R 3.5.

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
