[![R-CMD-check-prolfqua](https://github.com/fgcz/prolfqua/actions/workflows/r.yaml/badge.svg)](https://github.com/fgcz/prolfqua/actions/workflows/r.yaml) ![ReleseeDownloads](https://img.shields.io/github/downloads/fgcz/prolfqua/total)
[![codecov](https://codecov.io/gh/fgcz/prolfqua/branch/Modelling2R6/graph/badge.svg?token=NP7IPP323C)](https://codecov.io/gh/fgcz/prolfqua)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.1101%2F2022.06.07.494524-ligtgreen)](https://www.biorxiv.org/content/early/2022/06/09/2022.06.07.494524)

<img src="man/figures/imgfile.png" width="200">

# prolfqua - a R package for Proteomics Differential Expression Analysis

The R package contains functions for analyzing mass spectrometry based experiments.
This package is developed at the [FGCZ](http://fgcz.ch/).
The package documentation including vignettes can be accessed at https://fgcz.github.io/prolfqua/index.html

`prolfqua` makes easy things easy while remaining fully hackable.


# How to install prolfqua?

Requirements : A Windows|Linux|MacOSX platform with R (>= 4.1) installed.


We recommend to install the package using the latest [release](https://github.com/fgcz/prolfqua/releases)
Download the `prolfqua_X.Y.Z.tar.gz` from the [github release page](https://github.com/fgcz/prolfqua/releases). and then execute:

```
install.packages("prolfqua_X.Y.Z.tar.gz",repos = NULL, type="source")
```


To install the package without vignettes from github you can execute in R.

```
install.packages('remotes')
remotes::install_github('fgcz/prolfqua', dependencies = TRUE)
```


If you want to build the vignettes on your system:

```
install.packages('remotes')
remotes::install_github('fgcz/prolfqua', build_vignettes = TRUE, dependencies = TRUE)

```


Let us please know about any installation problems or errors when using the package:
https://github.com/fgcz/prolfqua/issues



# How to get started

- See [Bioconductor 2021 Conference poster](https://fgcz-proteomics.uzh.ch/~wolski/PosterBioconductor.html). 
- Watch the lightning (8 min) talk at [EuroBioc2020](https://www.youtube.com/watch?v=jOXU4X7nV9I&t) on YouTube or [slides](https://f1000research.com/slides/9-1476).
- See our article at [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.06.07.494524v1)
which describes our package.
- Read the pkgdown generate website https://fgcz.github.io/prolfqua/index.html

Or run the following example code:

```
R.version.string; packageVersion("prolfqua")

## read MQ peptide.txt and annotation table
startdata <- prolfqua::tidyMQ_Peptides(system.file(
  "samples/maxquant_txt/tiny2.zip", package = "prolfqua"))
annot <- readxl::read_xlsx(system.file(
  "samples/maxquant_txt/annotation_Ionstar2018_PXD003881.xlsx", package = "prolfqua"))
startdata <- dplyr::inner_join(annot, startdata, by = "raw.file")

## create MaxQuant configuration
config <- prolfqua::create_config_MQ_peptide()

## specify explanatory variable
config$table$factors[['dilution.']] = "sample"

## create R6 object
lfqpep <- prolfqua::LFQData$new(startdata, config, setup = TRUE) 

## remove observation with 0 intensity and filter for 2 peptides per protein
lfqpep$remove_small_intensities()$filter_proteins_by_peptide_count()

## transform intensities
lfqpep <- lfqpep$get_Transformer()$log2()$robscale()$lfq
lfqpep$rename_response("log_peptide_abundance")
agr <- lfqpep$get_Aggregator()
lfqpro <- agr$medpolish()
lfqpro$rename_response("log_protein_abundance")

## plot Figure 3 panels A-D
pl <- lfqpep$get_Plotter()
panelA <- pl$intensity_distribution_density() +
  ggplot2::labs(tag = "A") + ggplot2::theme(legend.position = "none")
panelB <- agr$plot()$plots[[54]] + ggplot2::labs(tag = "B")
panelC <- lfqpro$get_Stats()$violin() + ggplot2::labs(tag = "C")
pl <- lfqpro$get_Plotter()
panelD <- pl$boxplots()$boxplot[[54]] + ggplot2::labs(tag = "D")
ggpubr::ggarrange(panelA, panelB, panelC, panelD)

```


Detailed documentation with R code:

- [Comparing two Conditions](https://fgcz.github.io/prolfqua/articles/Comparing2Groups.html)
- [QC and protein wise sample size estimation](https://fgcz.github.io/prolfqua/articles/QualityControlAndSampleSizeEstimation.html)
- [Analysing factorial designs](https://fgcz.github.io/prolfqua/articles/Modelling2Factors.html)

Example QC and sample size report

- [QC and sample size Report](https://fgcz.github.io/prolfqua/articles/QCandSampleSize.html)

# How to cite?

Please do reference the [prolfqua article at biorxiv.org](https://www.biorxiv.org/content/10.1101/2022.06.07.494524v1)

```
@article {Wolski2022.06.07.494524,
	author = {Wolski, Witold Eryk and Nanni, Paolo and Grossmann, Jonas and d{\textquoteright}Errico, Maria and Schlapbach, Ralph and Panse, Christian},
	title = {prolfqua: A Comprehensive R-package for Proteomics Differential Expression Analysis},
	elocation-id = {2022.06.07.494524},
	year = {2022},
	doi = {10.1101/2022.06.07.494524},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2022/06/09/2022.06.07.494524},
	eprint = {https://www.biorxiv.org/content/early/2022/06/09/2022.06.07.494524.full.pdf},
	journal = {bioRxiv}
}
```

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

- [proDA](https://www.bioconductor.org/packages/release/bioc/html/proDA.html)
- [MSstats](https://www.bioconductor.org/packages/release/bioc/html/MSstats.html)
- [MSQRob](https://github.com/statOmics/MSqRob)
- [Triqler](https://github.com/statisticalbiotechnology/triqler)
- [DAPAR](https://github.com/samWieczorek/DAPAR/)
- [DAPARData](https://github.com/samWieczorek/DAPARdata/)
- [PECA/ROPECA](http://bioconductor.org/packages/release/bioc/html/PECA.html)

#  Relevant background information

- [R Companion](https://rcompanion.org/rcompanion/h_01.html)
- [Extending the Linear Model with R](http://www.maths.bath.ac.uk/~jjf23/ELM/)
- [Bayesian Data Analysis](http://www.stat.columbia.edu/~gelman/book/)
- [Bayesian essentials with R - R package](https://CRAN.R-project.org/package=bayess)
- [Contrasts in R - an example vignette by Rose Maier](https://rstudio-pubs-static.s3.amazonaws.com/65059_586f394d8eb84f84b1baaf56ffb6b47f.html)

# R packages to compute contrasts from linear models

- [emmeans](https://CRAN.R-project.org/package=emmeans) Obtain estimated marginal means (EMMs) for many linear, generalized linear, and mixed models.
- [lmerTest](https://CRAN.R-project.org/package=lmerTest) computes contrast for [lme4](https://CRAN.R-project.org/package=lme4) models
- [multcomp](https://CRAN.R-project.org/package=multcomp) computes contrast for linear models and adjusts p-values (multiple comparison)

# Future interesting topics or packages to look at

- [modelsummary](https://vincentarelbundock.github.io/modelsummary/index.html)
- [modelsummary tutorial](https://elbersb.com/public/pdf/web-7-regression-tables-graphs.pdf)
- [edgeR tutorial](https://gist.github.com/jdblischak/11384914)
- [another edgeR tutorial](https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html)

- https://fromthebottomoftheheap.net/2021/02/02/random-effects-in-gams/

# Sample size estimation based on FDR

- [ssize](https://www.bioconductor.org/packages/release/bioc/html/ssize.html)
- [ssize.fdr](https://CRAN.R-project.org/package=ssize.fdr)
  - related article [https://journal.r-project.org/archive/2009/RJ-2009-019/RJ-2009-019.pdf]
- [proper](https://bioconductor.org/packages/release/bioc/html/PROPER.html)

# What package name?

What name should we use?

https://twitter.com/WitoldE/status/1338799648149041156

- prolfqua - PROteomics Label Free QUAntification package (read prolewka)
- LFQService - we do proteomics LFQ services at the FGCZ.
- nalfqua - Not Another Label Free QUAntification package (read nalewka)
- prodea - proteomics differential expression analysis ?

