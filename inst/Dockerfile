FROM rocker/rstudio
MAINTAINER Christian Panse <Christian.Panse@gmail.com>

RUN apt-get update && apt-get install apt-utils curl libcurl4-openssl-dev libicu-dev pandoc wget unzip texlive texlive-latex-extra texlive-fonts-extra texlive-xetex pandoc-citeproc lmodern vim -y

LABEL description="1.0"
LABEL description="deploy WEW's sample size estimation."


RUN R -e "install.packages(c('BiocManager', 'remotes', 'rmarkdown', 'roxygen2', 'testthat', 'rcmdcheck'), repos='https://cloud.r-project.org'); pkgs <- c('DT', 'GGally', 'KernSmooth', 'MASS','bookdown', 'broom', 'coda', 'conflicted', 'corrplot', 'docopt', 'dplyr', 'flextable', 'forcats', 'ggbeeswarm', 'ggfortify', 'ggplot2', 'ggpubr', 'ggrepel', 'ggridges', 'glue', 'htmlwidgets', 'knitr', 'limma', 'lme4', 'lmerTest', 'magrittr', 'multcomp', 'pROC', 'pheatmap', 'plotly', 'plyr', 'progress', 'purrr', 'readr', 'readxl', 'rlang', 'rmarkdown', 'testthat', 'tibble', 'tidyr', 'tidyverse', 'yaml', 'writexl') ; lapply(pkgs, function(i){if (!i %in% installed.packages()){BiocManager::install(i)}}); "

RUN R -e "remotes::install_gitlab('wolski/prolfquadata', host='gitlab.bfabric.org')"

RUN mkdir /tmp/prolfqua
COPY . /tmp/prolfqua/
RUN cd /tmp/ && R CMD build prolfqua && R CMD INSTALL prolfqua*z

# ADD DOCKER FILE CODE HERE
RUN R -e "'install.packages('protViz')"

# DOCKERID=9e39049d81a5 && docker run -e PASSWORD=prolfqua -p 8787:8787 -e ADD=shiny -ti -v /scratch/:/scratch/  $DOCKERID
