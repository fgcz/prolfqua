FROM bioconductor/bioconductor_docker:RELEASE_3_22

LABEL description="Bioconductor check image for prolfqua"

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

WORKDIR /tmp/prolfqua

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
      make \
      qpdf \
      texinfo \
      texlive \
      texlive-bibtex-extra \
      texlive-fonts-extra \
      texlive-latex-extra \
      texlive-science \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

COPY DESCRIPTION ./

RUN Rscript -e "install.packages('remotes', repos = 'https://cloud.r-project.org')" \
    && Rscript -e "remotes::install_deps(dependencies = TRUE, upgrade = 'never')"

COPY . .

RUN mkdir -p check-logs \
    && make check 2>&1 | tee check-logs/check.log \
    && make check-bioc 2>&1 | tee check-logs/check-bioc.log
