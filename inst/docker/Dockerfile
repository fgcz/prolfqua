FROM rocker/r-ver:3.6.1
MAINTAINER Christian Panse <Christian.Panse@gmail.com>




RUN apt-get update && apt-get install apt-utils curl vim -y

RUN Rscript -e 'install.packages(c("bookdown", \
  "broom", \
  "conflicted", \
  "corrplot", \
  "dplyr", \
  "GGally", \
  "ggbeeswarm", "ggfortify", "ggplot2", "glue"), repos="https://stat.ethz.ch/CRAN/")' \
  && Rscript -e 'install.packages(c("heatmap3", "kableExtra", "lme4"), repos="https://stat.ethz.ch/CRAN/")' \
  && Rscript -e 'install.packages(c("lmerTest", "magrittr", "multcomp", "protViz", "purrr", "readxl", "tidyr"), repos="https://stat.ethz.ch/CRAN/")' \
  && Rscript -e 'install.packages(c("writexl", "yaml", "shiny", "testthat"), repos="https://stat.ethz.ch/CRAN/")' \
  && Rscript -e 'install.packages(c("tidyverse"), repos="https://stat.ethz.ch/CRAN/")'


# wew@fgcz.ethz.ch packages
RUN Rscript -e 'install.packages(c("quantable"), repos="https://stat.ethz.ch/CRAN/")'

#  plotly \
  #limma \


LABEL description="1.0"
LABEL description="deploy WEW's sample size estimation."
RUN Rscript -e 'install.packages(c("BiocManager"), repos="https://stat.ethz.ch/CRAN/")' \
  && Rscript -e 'BiocManager::install("limma")'
RUN  Rscript -e 'BiocManager::install(c("specL", "msqc1", "NestLink", "tartare"))'
