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
# RUN Rscript -e 'install.packages(c("quantable"), repos="https://stat.ethz.ch/CRAN/")'

#  plotly \
  #limma \

LABEL description="1.0"
LABEL description="deploy WEW's sample size estimation."
RUN Rscript -e 'install.packages(c("BiocManager"), repos="https://stat.ethz.ch/CRAN/")' \
  && Rscript -e 'BiocManager::install("limma")'
RUN  Rscript -e 'BiocManager::install(c("specL", "msqc1", "NestLink", "tartare"))'



RUN apt-get install libxml2-dev libcurl4-openssl-dev libssl-dev zlib1g-dev pandoc  libfontconfig1-dev libcairo2-dev -y
RUN Rscript -e 'install.packages(c("data.table", "kableExtra", "plotly", "tidyverse", "writexl"), repos="https://stat.ethz.ch/CRAN/")'
RUN Rscript -e 'install.packages(c("gdtools", "flextable"), repos="https://stat.ethz.ch/CRAN/")'

RUN mkdir /tmp/LFQService
COPY . /tmp/LFQService/
RUN cd /tmp/ && R CMD build LFQService --no-build-vignettes && R CMD INSTALL LFQS*z

## docker run -v /srv/www/htdocs//p2370/bfabric/Proteomics/MaxQuant/2019/2019-10/2019-10-19/workunit_222404/:/scratch/workunit_222404 -a stdin -a stdout -i -t d29638c801af Rscript /usr/local/lib/R/site-library/LFQService/run_scripts/lfq_MQ_SampleSizeReport.R /scratch/workunit_222404/1361948.zip --outdir=/tmp/WU222404/

