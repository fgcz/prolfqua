[![R-CMD-check-prolfqua](https://github.com/wolski/prolfqua/actions/workflows/r.yaml/badge.svg)](https://github.com/wolski/prolfqua/actions/workflows/r.yaml)

<img src="man/figures/imgfile.png" width="200">

# prolfqua - an R package for Proteomics Label Free Quantification Services

The R package contains functions for analyzing mass spectrometry based LFQ experiments.
This package is developed at the [FGCZ](http://www.fgcz.ch/).

# How to install prolfqua?

Requirements : A Windows/Linux/MacOSX x64 platform with R 4 or higher
To install the package please execute in R

```
install.packages('remotes')
remotes::install_github('wolski/prolfqua')
```

If you want to build the vignettes

```
install.packages('remotes')
remotes::install_gitlab("wolski/prolfquaData", host="gitlab.bfabric.org")
remotes::install_github('wolski/prolfqua', build_vignettes = TRUE)

```

Let us please know about any installation problems or errors when using the package
https://github.com/wolski/prolfqua/issues


# How to get started

See [Bioconductor 2021 Conference poster](https://fgcz-proteomics.uzh.ch/~wolski/PosterBioconductor.html). 
Watch the lightning (8 min) talk at [EuroBioc2020](https://www.youtube.com/watch?v=jOXU4X7nV9I&t) on YouTube.


Or read the pkgdown generate website https://wolski.github.io/prolfqua/index.html



Detailed documentation with R code:

- [Comparing two Conditions](https://wolski.github.io/prolfqua/articles/Comparing2Groups.html)
- [QC and sample size estimation](https://wolski.github.io/prolfqua/articles/QualityControlAndSampleSizeEstimation.html)
- [Analysing factorial designs](https://wolski.github.io/prolfqua/articles/Modelling2Factors.html)
- [Benchmarking LFQ pipeline with prolfqua](https://wolski.github.io/prolfqua/articles/BenchmarkingIonstarData.html)

Example QC and sample size report

- [QC and sample size Report](https://wolski.github.io/prolfqua/articles/QCandSampleSize.html)



# How to cite?

If you are using the package in your work please cite:
https://f1000research.com/slides/9-1476


## Motivation

The package for **pro**teomics **l**abel **f**ree **qua**ntification `prolfqua` (read : prolevka) evolved from a set of scripts and functions written in the R programming language to visualize and analyze mass spectrometric data, and some of them are still in R packages such as quantable, protViz or imsbInfer. For computing protein fold changes among treatment conditions, we first used t-test or linear models, then started to use functions implemented in the package limma to obtain moderated p-values. We did also try to use other packages such as MSStats, ROPECA or MSqRob all implemented in R, with the idea to integrate the various approaches to protein fold-change estimation. Although all these packages were written in R,  model specification, input and output formats differ widely and wildly, which made our aim to use the original implementations challenging. Therefore, and also to understand the algorithms used, we attempted to reimplement those methods, if possible. 

When developing _prolfqua_ we were inspired by packages such as _sf_ or _stars_ which use data in long table format and _dplyr_ for data transformation and ggplot2 for visualization.  In the long table format each column stores a different attribute, e.g. there is only a single column with the raw intensities. In the wide table format there might be several columns with the same attribute, e.g. for each recorded sample a raw intensity column.
In _prolfqua_ the data needed for analysis is represented using a single data-frame in long format and a configuration object. The configuration annotates the table, specifies what information is in which column. The results of the statistical modelling are stored in data frames.  Relying on the long data table format enabled us to access a large variety of useful visualizations as well as data preprocessing methods implemented in the R packages _dplyr_ and _ggplot2_.

The use of an annotated table makes integrating new data if provided in long formatted tables simple.  Hence for Spectronaut or Skyline text output, all is needed is a table annotation (see code snipped).  Since MSStats formatted input is a table in long format _prolefqa_ works with MSstats formatted files. For software, which writes the data in a wide table format, e.g. Maxquant, we implemented methods which first transform the data into a long format.  

A further design decision, which differentiates `prolfqua` is that it embraces and supports R's linear model formula interface, or R lme4 formula interface. R's formula interface for linear models is flexible, widely used and documented. The linear model and linear mixed model interfaces allow specifying a wide range of essential models, including parallel designs, factorial designs, repeated measurements and many more. Since `prolfqua` uses R modelling infrastructure directly, we can fit all these models to proteomics data.
This is not easily possible with any other package dedicated to proteomics data analysis. For instance, MSStats, although using the same modelling infrastructure, supports only a small subset of possible models. Limma, on the other hand, supports R formula interface but not for linear mixed models. Since the ROPECA package relies on _limma_ it is limited to the same subset of models. MSqRob is limited to random effects model's, and it is unclear how to fit these models to factorial designs, and how interactions among factors can be computed and tested.

The use of R's formula interface does not limit _prolfqua_ to the output provided by the R modelling infrastructure. _prolfqua_ also implements p-value moderations, as in the limma publication or computing probabilities of differential regulation, as suggested in the ROPECA publication. 
Moreover, the design decision to use the R formula interface allowed us to integrate Bayesian regression models provided by the r-package _brms_. Because of that, we can benchmark all those methods: linear models, mixed effect models, p-value moderation, ROPECA as well as Bayesian regression models within the same framework, which enabled us to evaluate the practical relevance of these methods.

Last but not least _prolfqua_ supports the LFQ data analysis workflow, e.g. computing coefficients of Variations (CV) for peptide and proteins, sample size estimation, visualization and summarization of missing data and intensity distributions, multivariate analysis of the data, etc.
It also implements various protein intensity summarization and inference methods, e.g. top 3, or Tukeys median polish etc. Last but not least, ANOVA analysis or model selection using the likelihood ratio test for thousand of proteins can be performed. 

To use `prolfqua` knowledge of the R regression model infrastructure is of advantage. Acknowledging, the complexity of the formula interface,  we provide an  MSstats emulator, where the model specification is generated based on the annotation file structure. 



## Running R-scripts

Generate QC report from bat file.
First add `<prolfqua_path>/win` to your path variable. Then you can generate a QC report from a maxquant QC by running.


```
lfq_MQ_SampleSizeReport.bat .\data\1296877_QC.zip
```


```
Rscript ~/__checkouts/R/prolfqua/inst/run_scripts/lfq_MQ_SampleSizeReport.R --help
Rscript ~/__checkouts/R/prolfqua/inst/run_scripts/lfq_MQ_SampleSizeReport.R ~/Downloads/1330043.zip
```


# Related resources

- [Triqler](https://github.com/statisticalbiotechnology/triqler)
- [MSQRob](https://github.com/statOmics/MSqRob)
- [PECA/ROPECA](http://bioconductor.org/packages/release/bioc/html/PECA.html)
- [MSstats](https://www.bioconductor.org/packages/release/bioc/html/MSstats.html)
- PEPA test - analysis including shared peptides
  - [DAPAR](https://github.com/samWieczorek/DAPAR/)
  - [DAPARData](https://github.com/samWieczorek/DAPARdata/)
  - [PEPA validation](https://github.com/ThomasBurger/pepa-validation)

#  Relevant background information

- [Contrasts in R by Rose Maier](https://rstudio-pubs-static.s3.amazonaws.com/65059_586f394d8eb84f84b1baaf56ffb6b47f.html)
- [R Companion](https://rcompanion.org/rcompanion/h_01.html)
- [Extending the Linear Model with R](http://www.maths.bath.ac.uk/~jjf23/ELM/)
- [Bayesian Data Analysis](http://www.stat.columbia.edu/~gelman/book/)
- [Bayesian essentials with R - R package](https://cran.r-project.org/web/packages/bayess/index.html)
- [Contrasts in R - an example vignette by Rose Maier](https://rstudio-pubs-static.s3.amazonaws.com/65059_586f394d8eb84f84b1baaf56ffb6b47f.html)

# R packages to compute contrasts from linear models

- [emmeans](https://CRAN.R-project.org/package=emmeans) Obtain estimated marginal means (EMMs) for many linear, generalized linear, and mixed models.
- [lmerTest](https://CRAN.R-project.org/package=lmerTest) computes contrast for [lme4](https://CRAN.R-project.org/package=lme4) models


# Future interesting topics or packages to look at

- [modelsummary](https://vincentarelbundock.github.io/modelsummary/index.html)
- [modelsummary tutorial](https://elbersb.com/public/pdf/web-7-regression-tables-graphs.pdf)
- [edgeR tutorial](https://gist.github.com/jdblischak/11384914)
- [another edgeR tutorial](https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html)

# Sample size estimation based on FDR

- [ssize](https://www.bioconductor.org/packages/release/bioc/html/ssize.html)
- [ssize.fdr](https://cran.r-project.org/web/packages/ssize.fdr/index.html)
  - related article [https://journal.r-project.org/archive/2009/RJ-2009-019/RJ-2009-019.pdf]
- [proper](https://bioconductor.org/packages/release/bioc/html/PROPER.html)

# What package name?

What name should we use?

https://twitter.com/WitoldE/status/1338799648149041156

- prolfqua - PROteomics Label Free QUAntification package (read prolewka)

- LFQService - we do proteomics LFQ services at the FGCZ.
- nalfqua - Not Another Label Free QUAntification package (read nalewka)

